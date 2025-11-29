"""
Main pipeline for building the cross-species antibody-antigen structure dataset.

This module orchestrates:
1. SAbDab data download and filtering
2. UniProt/Ensembl ortholog mapping
3. PDB structure retrieval
4. Structure processing and alignment
5. Output generation
"""

import os
import json
import time
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime

import gemmi
import pandas as pd
from tqdm import tqdm

from .sabdab import (
    download_sabdab_summary,
    parse_sabdab_summary,
    filter_human_antigen_complexes,
    get_unique_complexes
)
from .mapping import (
    get_uniprot_from_pdb_chain,
    get_uniprot_info,
    find_mouse_ortholog,
    find_pdb_structures_for_uniprot,
    get_pdb_chain_for_entity
)
from .structure import (
    download_structure,
    parse_structure,
    save_structure,
    save_structure_cif,
    save_cif_with_metadata,
    align_structures_pymol,
    align_structures_biopython,
    merge_structures,
    ChainSelect,
    PYMOL_AVAILABLE
)


@dataclass
class DataPoint:
    """Represents a single data point in the dataset."""
    id: str
    pdb_id_human: str
    pdb_id_mouse: str
    human_antigen_chain: str
    mouse_antigen_chain: str
    heavy_chain: str
    light_chain: str
    human_uniprot: str
    mouse_uniprot: str
    gene_name: str
    protein_name: str
    sequence_identity: float
    resolution_human: Optional[float] = None
    resolution_mouse: Optional[float] = None
    alignment_rmsd: Optional[float] = None
    status: str = "pending"
    error_message: str = ""
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    human_antigen_chains: List[str] = field(default_factory=list)


class CrossSpeciesDatasetPipeline:
    """
    Pipeline for building cross-species antibody-antigen structure dataset.
    """

    def __init__(
        self,
        data_dir: str,
        output_dir: str,
        resolution_threshold: float = 2.5,
        sequence_identity_threshold: float = 50.0,
        use_pymol: bool = True
    ):
        """
        Initialize the pipeline.

        Args:
            data_dir: Directory for intermediate data (downloads, cache)
            output_dir: Directory for final output
            resolution_threshold: Maximum resolution in Angstroms
            sequence_identity_threshold: Minimum sequence identity percentage
            use_pymol: Use PyMOL for alignment (falls back to Biopython if unavailable)
        """
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.resolution_threshold = resolution_threshold
        self.sequence_identity_threshold = sequence_identity_threshold
        self.use_pymol = use_pymol and PYMOL_AVAILABLE

        # Create directories
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(data_dir, "pdb_cache"), exist_ok=True)

        # State tracking
        self.data_points: List[DataPoint] = []
        self.processing_log: List[Dict] = []

    def log(self, message: str, level: str = "INFO"):
        """Log a message."""
        timestamp = datetime.now().isoformat()
        log_entry = {"timestamp": timestamp, "level": level, "message": message}
        self.processing_log.append(log_entry)
        print(f"[{timestamp}] {level}: {message}")

    def run(self, limit: Optional[int] = None, dry_run: bool = False) -> pd.DataFrame:
        """
        Run the full pipeline.

        Args:
            limit: Maximum number of data points to process (None = all)
            dry_run: If True, only analyze without downloading structures

        Returns:
            DataFrame with processing results
        """
        self.log("Starting cross-species dataset pipeline")

        # Step 1: Download and parse SAbDab data
        self.log("Step 1: Downloading and parsing SAbDab data")
        sabdab_df = self._load_sabdab_data()

        # Step 2: Filter for suitable complexes
        self.log("Step 2: Filtering for human antigen complexes")
        filtered_df = filter_human_antigen_complexes(
            sabdab_df,
            resolution_threshold=self.resolution_threshold
        )
        unique_df = get_unique_complexes(filtered_df)

        if limit:
            unique_df = unique_df.head(limit)
            self.log(f"Limited to {limit} entries for processing")

        # Step 3: Process each complex
        self.log(f"Step 3: Processing {len(unique_df)} complexes")

        successful = 0
        failed = 0

        for idx, row in tqdm(unique_df.iterrows(), total=len(unique_df), desc="Processing"):
            try:
                data_point = self._process_complex(row, dry_run=dry_run)
                if data_point:
                    self.data_points.append(data_point)
                    if data_point.status == "success":
                        successful += 1
                    else:
                        failed += 1
            except Exception as e:
                self.log(f"Error processing {row['pdb']}: {e}", level="ERROR")
                failed += 1

            # Rate limiting
            time.sleep(0.5)

        self.log(f"Processing complete: {successful} successful, {failed} failed")

        # Step 4: Generate summary
        self.log("Step 4: Generating summary CSV")
        summary_df = self._generate_summary()

        return summary_df

    def _load_sabdab_data(self) -> pd.DataFrame:
        """Load SAbDab data, downloading if necessary."""
        summary_path = download_sabdab_summary(self.data_dir)
        return parse_sabdab_summary(summary_path)

    def _process_complex(self, row: pd.Series, dry_run: bool = False) -> Optional[DataPoint]:
        """
        Process a single antibody-antigen complex.

        Args:
            row: SAbDab entry row
            dry_run: If True, only analyze without downloading structures

        Returns:
            DataPoint or None
        """
        pdb_id = row['pdb']
        antigen_chain_raw = str(row['antigen_chain'])
        antigen_chains = [c.strip() for c in antigen_chain_raw.split('|') if c.strip()]
        antigen_chain = antigen_chains[0] if antigen_chains else antigen_chain_raw.strip()
        heavy_chain = row['Hchain']
        light_chain = row['Lchain']

        # Step 1: Get UniProt ID for human antigen
        human_uniprot = get_uniprot_from_pdb_chain(pdb_id, antigen_chain)
        if not human_uniprot:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse="",
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain="",
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot="",
                mouse_uniprot="",
                gene_name="",
                protein_name="",
                sequence_identity=0.0,
                status="failed",
                error_message="Could not map antigen chain to UniProt"
            )

        # Step 2: Get human protein info
        human_info = get_uniprot_info(human_uniprot)
        if not human_info:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse="",
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain="",
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot="",
                gene_name="",
                protein_name="",
                sequence_identity=0.0,
                status="failed",
                error_message="Could not retrieve UniProt info"
            )

        # Step 3: Find mouse ortholog
        mouse_info = find_mouse_ortholog(human_uniprot)
        if not mouse_info:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse="",
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain="",
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot="",
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=0.0,
                status="failed",
                error_message="No mouse ortholog found"
            )

        # Check sequence identity threshold
        seq_identity = mouse_info.get('sequence_identity', 0.0)
        if seq_identity < self.sequence_identity_threshold:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse="",
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain="",
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_info.get('accession', ''),
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                status="failed",
                error_message=f"Sequence identity {seq_identity:.1f}% below threshold {self.sequence_identity_threshold}%"
            )

        # Step 4: Find PDB structure for mouse ortholog
        mouse_structures = find_pdb_structures_for_uniprot(mouse_info['accession'])
        if not mouse_structures:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse="",
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain="",
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_info.get('accession', ''),
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                status="failed",
                error_message="No PDB structure found for mouse ortholog"
            )

        # Select best mouse structure (first one for now, could be improved)
        mouse_structure_info = mouse_structures[0]
        mouse_pdb_id = mouse_structure_info['pdb_id']

        # Get chain ID for mouse structure
        mouse_chains = get_pdb_chain_for_entity(
            mouse_pdb_id,
            mouse_structure_info['entity_id']
        )
        mouse_chain = mouse_chains[0] if mouse_chains else 'A'

        if dry_run:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse=mouse_pdb_id,
                human_antigen_chain=antigen_chain,
                human_antigen_chains=antigen_chains,
                mouse_antigen_chain=mouse_chain,
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_info.get('accession', ''),
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                status="dry_run",
                error_message=""
            )

        # Step 5: Download and process structures
        try:
            rmsd = self._process_structures(
                pdb_id=pdb_id,
                antigen_chain=antigen_chain,
                antigen_chains=antigen_chains,
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                mouse_pdb_id=mouse_pdb_id,
                mouse_chain=mouse_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_info.get('accession', ''),
                gene_name=human_info.get('gene_name', '')
            )

            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse=mouse_pdb_id,
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain=mouse_chain,
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_info.get('accession', ''),
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                alignment_rmsd=rmsd,
                status="success",
                error_message=""
            )

        except Exception as e:
            return DataPoint(
                id=f"DP_{pdb_id}_{antigen_chain}",
                pdb_id_human=pdb_id,
                pdb_id_mouse=mouse_pdb_id,
                human_antigen_chain=antigen_chain,
                mouse_antigen_chain=mouse_chain,
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_info.get('accession', ''),
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                status="failed",
                error_message=f"Structure processing error: {str(e)}"
            )

    def _process_structures(
        self,
        pdb_id: str,
        antigen_chain: str,
        antigen_chains: List[str],
        heavy_chain: str,
        light_chain: str,
        mouse_pdb_id: str,
        mouse_chain: str,
        human_uniprot: str,
        mouse_uniprot: str,
        gene_name: str
    ) -> float:
        """
        Download, process and align structures.

        Returns:
            Alignment RMSD
        """
        cache_dir = os.path.join(self.data_dir, "pdb_cache")

        # Create output folder for this data point
        dp_id = f"DP_{pdb_id}_{antigen_chain}"
        dp_dir = os.path.join(self.output_dir, dp_id)
        os.makedirs(dp_dir, exist_ok=True)

        # Download human complex structure
        human_cif, human_pdb = download_structure(pdb_id, cache_dir)
        if not human_cif and not human_pdb:
            raise ValueError(f"Could not download structure for {pdb_id}")

        # Download mouse antigen structure
        mouse_cif, mouse_pdb = download_structure(mouse_pdb_id, cache_dir)
        if not mouse_cif and not mouse_pdb:
            raise ValueError(f"Could not download structure for {mouse_pdb_id}")

        # Parse structures
        human_file = human_cif or human_pdb
        mouse_file = mouse_cif or mouse_pdb

        human_structure = parse_structure(human_file, pdb_id)
        mouse_structure = parse_structure(mouse_file, mouse_pdb_id)

        if not human_structure or not mouse_structure:
            raise ValueError("Could not parse structures")

        # Extract and save antibody (H + L chains)
        antibody_chains = []
        if heavy_chain and heavy_chain != 'NA':
            antibody_chains.append(heavy_chain)
        if light_chain and light_chain != 'NA':
            antibody_chains.append(light_chain)

        # Save antibody
        antibody_pdb = os.path.join(dp_dir, f"{dp_id}_antibody.pdb")
        antibody_cif = os.path.join(dp_dir, f"{dp_id}_antibody.cif")
        save_structure(human_structure, antibody_pdb, chain_ids=antibody_chains)
        if human_cif:
            save_cif_with_metadata(human_cif, antibody_cif, chain_ids=antibody_chains)
        else:
            save_structure_cif(human_structure, antibody_cif, chain_ids=antibody_chains)

        # Save human antigen
        human_ag_pdb = os.path.join(dp_dir, f"{dp_id}_human_ag.pdb")
        human_ag_cif = os.path.join(dp_dir, f"{dp_id}_human_ag.cif")
        human_ag_full_pdb = None
        human_ag_full_cif = None
        save_structure(human_structure, human_ag_pdb, chain_ids=[antigen_chain])
        if human_cif:
            save_cif_with_metadata(human_cif, human_ag_cif, chain_ids=[antigen_chain])
        else:
            save_structure_cif(human_structure, human_ag_cif, chain_ids=[antigen_chain])

        # Save full human antigen (all listed antigen chains)
        if antigen_chains:
            human_ag_full_pdb = os.path.join(dp_dir, f"{dp_id}_human_ag_full.pdb")
            human_ag_full_cif = os.path.join(dp_dir, f"{dp_id}_human_ag_full.cif")
            save_structure(human_structure, human_ag_full_pdb, chain_ids=antigen_chains)
            if human_cif:
                save_cif_with_metadata(human_cif, human_ag_full_cif, chain_ids=antigen_chains)
            else:
                save_structure_cif(human_structure, human_ag_full_cif, chain_ids=antigen_chains)

        # Align mouse antigen to human antigen position
        mouse_ag_pdb = os.path.join(dp_dir, f"{dp_id}_mouse_ag.pdb")
        mouse_ag_cif = os.path.join(dp_dir, f"{dp_id}_mouse_ag.cif")
        mouse_ag_full_pdb = os.path.join(dp_dir, f"{dp_id}_mouse_ag_full.pdb")
        mouse_ag_full_cif = os.path.join(dp_dir, f"{dp_id}_mouse_ag_full.cif")
        human_complex_pdb = os.path.join(dp_dir, f"{dp_id}_human_complex.pdb")
        human_complex_cif = os.path.join(dp_dir, f"{dp_id}_human_complex.cif")
        mouse_complex_pdb = os.path.join(dp_dir, f"{dp_id}_mouse_complex.pdb")
        mouse_complex_cif = os.path.join(dp_dir, f"{dp_id}_mouse_complex.cif")

        # First save unaligned mouse structure
        save_structure(mouse_structure, mouse_ag_pdb, chain_ids=[mouse_chain])

        # Perform alignment
        rotran = None
        aligned_structure = mouse_structure
        if self.use_pymol:
            # Use PyMOL for alignment
            rmsd = align_structures_pymol(
                mobile_file=mouse_ag_pdb,
                reference_file=human_ag_pdb,
                output_file=mouse_ag_pdb,
                mobile_chain=mouse_chain,
                reference_chain=antigen_chain
            )
            aligned_structure = parse_structure(mouse_ag_pdb, mouse_pdb_id) or mouse_structure
        else:
            # Use Biopython for alignment
            rmsd, aligned_structure, rotran = align_structures_biopython(
                mobile_structure=mouse_structure,
                reference_structure=human_structure,
                mobile_chain=mouse_chain,
                reference_chain=antigen_chain
            )
            save_structure(aligned_structure, mouse_ag_pdb, chain_ids=[mouse_chain])

        # Save aligned mouse structure as CIF too
        aligned_mouse = parse_structure(mouse_ag_pdb, "aligned_mouse")
        if aligned_mouse:
            if mouse_cif and rotran:
                # Apply the alignment transform to original mmCIF so metadata is preserved
                save_cif_with_metadata(
                    mouse_cif,
                    mouse_ag_cif,
                    chain_ids=[mouse_chain],
                    rotran=rotran
                )
            elif mouse_cif:
                # No transform available (e.g., PyMOL path); fall back to aligned structure output
                save_structure_cif(aligned_mouse, mouse_ag_cif)
            else:
                save_structure_cif(aligned_mouse, mouse_ag_cif)

        # Save aligned mouse full structure (all chains)
        if rotran and mouse_cif:
            save_cif_with_metadata(
                mouse_cif,
                mouse_ag_full_cif,
                chain_ids=None,
                rotran=rotran
            )
        elif mouse_cif:
            save_structure_cif(aligned_structure, mouse_ag_full_cif)
        elif aligned_mouse:
            save_structure_cif(aligned_mouse, mouse_ag_full_cif)

        save_structure(aligned_structure, mouse_ag_full_pdb)

        # Build complexes: antibody + antigen (human) and antibody + aligned mouse antigen
        complex_chain_ids_human = antibody_chains + (antigen_chains if antigen_chains else [antigen_chain])
        # Human complex PDB/CIF (metadata from human cif)
        save_structure(human_structure, human_complex_pdb, chain_ids=complex_chain_ids_human)
        if human_cif:
            save_cif_with_metadata(human_cif, human_complex_cif, chain_ids=complex_chain_ids_human)
        else:
            save_structure_cif(human_structure, human_complex_cif, chain_ids=complex_chain_ids_human)

        # Mouse complex: merge antibody (human structure) + aligned mouse antigen structure
        merged_struct = merge_structures(
            [
                (human_structure, antibody_chains),
                (aligned_structure, None)
            ],
            new_id=f"{dp_id}_mouse_complex"
        )
        save_structure(merged_struct, mouse_complex_pdb)
        # For CIF, best-effort: combine using antibody metadata; antigen metadata retained if available
        if human_cif and mouse_cif and rotran:
            # Start from human CIF (for antibody metadata) and append transformed mouse atom rows
            try:
                doc_h = gemmi.cif.read_file(human_cif)
                block_h = doc_h.sole_block()
                table = block_h.find_mmcif_category("_atom_site")
                if not table:
                    raise ValueError("No atom_site in human CIF")
                tags = list(table.tags)
                idx_label = tags.index("_atom_site.label_asym_id") if "_atom_site.label_asym_id" in tags else -1
                idx_auth = tags.index("_atom_site.auth_asym_id") if "_atom_site.auth_asym_id" in tags else -1

                # Keep antibody rows only
                to_remove = []
                for i, row in enumerate(table):
                    label_chain = row[idx_label] if idx_label >= 0 else ""
                    auth_chain = row[idx_auth] if idx_auth >= 0 else ""
                    if label_chain not in antibody_chains and auth_chain not in antibody_chains:
                        to_remove.append(i)
                for offset, i in enumerate(to_remove):
                    table.remove_row(i - offset)

                # Append mouse rows (all chains) with transform
                doc_m = gemmi.cif.read_file(mouse_cif)
                block_m = doc_m.sole_block()
                table_m = block_m.find_mmcif_category("_atom_site")
                if not table_m:
                    raise ValueError("No atom_site in mouse CIF")
                tags_m = list(table_m.tags)
                idx_x = tags.index("_atom_site.Cartn_x") if "_atom_site.Cartn_x" in tags else -1
                idx_y = tags.index("_atom_site.Cartn_y") if "_atom_site.Cartn_y" in tags else -1
                idx_z = tags.index("_atom_site.Cartn_z") if "_atom_site.Cartn_z" in tags else -1
                idx_x_m = tags_m.index("_atom_site.Cartn_x") if "_atom_site.Cartn_x" in tags_m else -1
                idx_y_m = tags_m.index("_atom_site.Cartn_y") if "_atom_site.Cartn_y" in tags_m else -1
                idx_z_m = tags_m.index("_atom_site.Cartn_z") if "_atom_site.Cartn_z" in tags_m else -1
                rot, tran = rotran
                for row in table_m:
                    new_row = ["" for _ in range(len(tags))]
                    for j, tag in enumerate(tags_m):
                        try:
                            tgt_idx = tags.index(tag)
                        except ValueError:
                            continue
                        new_row[tgt_idx] = row[j]
                    if idx_x >= 0 and idx_y >= 0 and idx_z >= 0 and idx_x_m >= 0 and idx_y_m >= 0 and idx_z_m >= 0:
                        try:
                            x = float(row[idx_x_m])
                            y = float(row[idx_y_m])
                            z = float(row[idx_z_m])
                            new_x = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z + tran[0]
                            new_y = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z + tran[1]
                            new_z = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z + tran[2]
                            new_row[idx_x] = f"{new_x:.3f}"
                            new_row[idx_y] = f"{new_y:.3f}"
                            new_row[idx_z] = f"{new_z:.3f}"
                        except Exception:
                            pass
                    table.append_row(new_row)

                block_h.check_empty_loops("_atom_site")
                doc_h.write_file(mouse_complex_cif)
            except Exception as e:
                self.log(f"Mouse complex CIF fallback for {dp_id}: {e}", level="ERROR")
                save_structure_cif(merged_struct, mouse_complex_cif)
        else:
            save_structure_cif(merged_struct, mouse_complex_cif)

        # Save metadata
        metadata = {
            'id': dp_id,
            'pdb_id_human': pdb_id,
            'pdb_id_mouse': mouse_pdb_id,
            'human_antigen_chain': antigen_chain,
            'human_antigen_chains': antigen_chains,
            'mouse_antigen_chain': mouse_chain,
            'antibody_chains': antibody_chains,
            'human_uniprot': human_uniprot,
            'mouse_uniprot': mouse_uniprot,
            'gene_name': gene_name,
            'alignment_rmsd': rmsd,
            'files': {
                'antibody_pdb': os.path.basename(antibody_pdb),
                'antibody_cif': os.path.basename(antibody_cif),
                'human_ag_pdb': os.path.basename(human_ag_pdb),
                'human_ag_cif': os.path.basename(human_ag_cif),
                'human_ag_full_pdb': os.path.basename(human_ag_full_pdb) if antigen_chains else "",
                'human_ag_full_cif': os.path.basename(human_ag_full_cif) if antigen_chains else "",
                'mouse_ag_pdb': os.path.basename(mouse_ag_pdb),
                'mouse_ag_cif': os.path.basename(mouse_ag_cif),
                'mouse_ag_full_pdb': os.path.basename(mouse_ag_full_pdb),
                'mouse_ag_full_cif': os.path.basename(mouse_ag_full_cif),
                'human_complex_pdb': os.path.basename(human_complex_pdb),
                'human_complex_cif': os.path.basename(human_complex_cif),
                'mouse_complex_pdb': os.path.basename(mouse_complex_pdb),
                'mouse_complex_cif': os.path.basename(mouse_complex_cif)
            }
        }

        with open(os.path.join(dp_dir, 'metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2)

        return rmsd

    def _generate_summary(self) -> pd.DataFrame:
        """Generate summary CSV file."""
        if not self.data_points:
            return pd.DataFrame()

        # Convert data points to DataFrame
        data = [asdict(dp) for dp in self.data_points]
        df = pd.DataFrame(data)

        # Save to CSV
        csv_path = os.path.join(self.output_dir, 'dataset_summary.csv')
        df.to_csv(csv_path, index=False)
        self.log(f"Summary saved to {csv_path}")

        # Print statistics
        total = len(df)
        successful = len(df[df['status'] == 'success'])
        failed = len(df[df['status'] == 'failed'])

        self.log(f"Dataset statistics:")
        self.log(f"  Total entries processed: {total}")
        self.log(f"  Successful: {successful}")
        self.log(f"  Failed: {failed}")

        if successful > 0:
            avg_identity = df[df['status'] == 'success']['sequence_identity'].mean()
            avg_rmsd = df[df['status'] == 'success']['alignment_rmsd'].mean()
            self.log(f"  Average sequence identity: {avg_identity:.1f}%")
            if pd.notna(avg_rmsd):
                self.log(f"  Average alignment RMSD: {avg_rmsd:.2f} A")

        # Save processing log
        log_path = os.path.join(self.output_dir, 'processing_log.json')
        with open(log_path, 'w') as f:
            json.dump(self.processing_log, f, indent=2)

        return df


if __name__ == "__main__":
    # Test run
    import tempfile

    data_dir = tempfile.mkdtemp()
    output_dir = tempfile.mkdtemp()

    pipeline = CrossSpeciesDatasetPipeline(
        data_dir=data_dir,
        output_dir=output_dir
    )

    # Dry run with limit
    result = pipeline.run(limit=5, dry_run=True)
    print(result)
