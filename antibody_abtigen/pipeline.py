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
    get_unique_complexes,
    save_filtered_summary,
    find_mouse_antigen_in_sabdab
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
    clean_structure,
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
    # New fields for extension modules
    epitope_identity: Optional[float] = None      # Epitope sequence identity (%)
    epitope_rmsd: Optional[float] = None          # Epitope structure RMSD (Å)
    epitope_consistent: Optional[bool] = None     # Whether passes consistency check
    mouse_structure_source: str = 'RCSB'          # 'SAbDab' | 'RCSB'
    mouse_gene_name: str = ""                     # Mouse gene name for grouping


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
        use_pymol: bool = True,
        epitope_rmsd_threshold: float = 1.5,
        epitope_identity_threshold: float = 80.0,
        prefer_sabdab_mouse: bool = True
    ):
        """
        Initialize the pipeline.

        Args:
            data_dir: Directory for intermediate data (downloads, cache)
            output_dir: Directory for final output
            resolution_threshold: Maximum resolution in Angstroms
            sequence_identity_threshold: Minimum sequence identity percentage
            use_pymol: Use PyMOL for alignment (falls back to Biopython if unavailable)
            epitope_rmsd_threshold: Epitope RMSD threshold for consistency (Å)
            epitope_identity_threshold: Epitope sequence identity threshold (%)
            prefer_sabdab_mouse: Prefer mouse structures from SAbDab over RCSB
        """
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.resolution_threshold = resolution_threshold
        self.sequence_identity_threshold = sequence_identity_threshold
        self.use_pymol = use_pymol and PYMOL_AVAILABLE
        self.epitope_rmsd_threshold = epitope_rmsd_threshold
        self.epitope_identity_threshold = epitope_identity_threshold
        self.prefer_sabdab_mouse = prefer_sabdab_mouse

        # New directory structure (following document specification)
        # All data goes under data_dir (./data/antigen_antibody/)
        # SAbDab directories (Phase 1-3)
        self.sabdab_raw_dir = os.path.join(data_dir, "SAbDab", "raw")
        self.sabdab_cleaned_dir = os.path.join(data_dir, "SAbDab", "cleaned")
        self.sabdab_antigen_dir = os.path.join(data_dir, "SAbDab", "antigen")
        self.sabdab_antibody_dir = os.path.join(data_dir, "SAbDab", "antibody")

        # Phase 4-5 directories (also under data_dir per document spec)
        self.human_mouse_dir = os.path.join(data_dir, "HumanMouse")
        self.mouse_antigen_dir = os.path.join(data_dir, "MouseAntigen")

        # Create all directories
        for d in [
            self.sabdab_raw_dir, self.sabdab_cleaned_dir,
            self.sabdab_antigen_dir, self.sabdab_antibody_dir,
            self.human_mouse_dir, self.mouse_antigen_dir
        ]:
            os.makedirs(d, exist_ok=True)

        # State tracking
        self.data_points: List[DataPoint] = []
        self.processing_log: List[Dict] = []
        self.sabdab_full_df: Optional[pd.DataFrame] = None  # Store full SAbDab data for mouse search

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

        # Step 4: Build mouse antigen library (Phase 5)
        if not dry_run and successful > 0:
            self.log("Step 4: Building mouse antigen library")
            self._build_mouse_antigen_library()
            self._update_pairs_csv()

        # Step 5: Generate summary
        self.log("Step 5: Generating summary CSV")
        summary_df = self._generate_summary()

        return summary_df

    def _load_sabdab_data(self) -> pd.DataFrame:
        """Load SAbDab data, downloading if necessary."""
        # Download to SAbDab/raw/ directory
        summary_path = download_sabdab_summary(self.sabdab_raw_dir)
        full_df = parse_sabdab_summary(summary_path)
        # Store full dataframe for mouse antigen search
        self.sabdab_full_df = full_df
        return full_df

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

        # Step 3: Try to find mouse structure (prioritize SAbDab, fallback to RCSB)
        mouse_pdb_id = None
        mouse_chain = None
        mouse_uniprot = ""
        mouse_gene_name = ""
        mouse_structure_source = "RCSB"
        seq_identity = 0.0

        # Step 3a: Try SAbDab first if enabled
        antigen_name = row.get('antigen_name', '')
        if self.prefer_sabdab_mouse and self.sabdab_full_df is not None:
            sabdab_mouse = find_mouse_antigen_in_sabdab(
                self.sabdab_full_df,
                antigen_name,
                human_uniprot,
                resolution_threshold=self.resolution_threshold
            )
            if sabdab_mouse:
                mouse_pdb_id = sabdab_mouse['pdb']
                # Parse multi-chain format like "A | B" - take first chain
                raw_chain = sabdab_mouse['antigen_chain']
                if ' | ' in str(raw_chain):
                    mouse_chain = raw_chain.split(' | ')[0].strip()
                else:
                    mouse_chain = str(raw_chain).strip()
                mouse_gene_name = sabdab_mouse.get('antigen_name', '')
                mouse_structure_source = "SAbDab"
                seq_identity = 100.0  # Assume high identity for same-name antigen
                self.log(f"Found mouse antigen in SAbDab: {mouse_pdb_id} (chain {mouse_chain})")

        # Step 3b: Fallback to UniProt/Ensembl ortholog mapping
        if not mouse_pdb_id:
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
            mouse_uniprot = mouse_info.get('accession', '')
            mouse_gene_name = mouse_info.get('gene_name', '')
            mouse_structures = find_pdb_structures_for_uniprot(mouse_uniprot)
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
                    mouse_uniprot=mouse_uniprot,
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
            mouse_chains_list = get_pdb_chain_for_entity(
                mouse_pdb_id,
                mouse_structure_info['entity_id']
            )
            mouse_chain = mouse_chains_list[0] if mouse_chains_list else 'A'

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
                mouse_uniprot=mouse_uniprot,
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                mouse_structure_source=mouse_structure_source,
                mouse_gene_name=mouse_gene_name,
                status="dry_run",
                error_message=""
            )

        # Step 5: Download and process structures
        try:
            rmsd, pair_metadata = self._process_structures(
                pdb_id=pdb_id,
                antigen_chain=antigen_chain,
                antigen_chains=antigen_chains,
                heavy_chain=heavy_chain,
                light_chain=light_chain,
                mouse_pdb_id=mouse_pdb_id,
                mouse_chain=mouse_chain,
                human_uniprot=human_uniprot,
                mouse_uniprot=mouse_uniprot,
                gene_name=human_info.get('gene_name', ''),
                mouse_gene_name=mouse_gene_name,
                mouse_structure_source=mouse_structure_source
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
                mouse_uniprot=mouse_uniprot,
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                alignment_rmsd=rmsd,
                mouse_structure_source=mouse_structure_source,
                mouse_gene_name=mouse_gene_name,
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
                mouse_uniprot=mouse_uniprot,
                gene_name=human_info.get('gene_name', ''),
                protein_name=human_info.get('protein_name', ''),
                sequence_identity=seq_identity,
                mouse_structure_source=mouse_structure_source,
                mouse_gene_name=mouse_gene_name,
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
        gene_name: str,
        mouse_gene_name: str = "",
        mouse_structure_source: str = "RCSB"
    ) -> Tuple[float, Dict]:
        """
        Download, process and align structures following the document's 5-phase structure.

        Phase 1: Download raw structures to SAbDab/raw/
        Phase 2: Clean structures to SAbDab/cleaned/
        Phase 3: Split antigen/antibody to SAbDab/antigen/ and SAbDab/antibody/
        Phase 4: Align human-mouse pairs to HumanMouse/{pair}/
        Phase 5: Build mouse antigen library (done in separate method)

        Returns:
            Tuple of (RMSD, metadata dict)
        """
        # === Phase 1: Download raw structures ===
        # Download to SAbDab/raw/
        human_cif, human_pdb = download_structure(pdb_id, self.sabdab_raw_dir)
        if not human_cif and not human_pdb:
            raise ValueError(f"Could not download structure for {pdb_id}")

        mouse_cif, mouse_pdb = download_structure(mouse_pdb_id, self.sabdab_raw_dir)
        if not mouse_cif and not mouse_pdb:
            raise ValueError(f"Could not download structure for {mouse_pdb_id}")

        # Parse structures
        human_file = human_cif or human_pdb
        mouse_file = mouse_cif or mouse_pdb

        human_structure = parse_structure(human_file, pdb_id)
        mouse_structure = parse_structure(mouse_file, mouse_pdb_id)

        if not human_structure or not mouse_structure:
            raise ValueError("Could not parse structures")

        # === Phase 2: Clean structures ===
        # Determine chains to keep
        antibody_chains = []
        if heavy_chain and heavy_chain != 'NA':
            antibody_chains.append(heavy_chain)
        if light_chain and light_chain != 'NA':
            antibody_chains.append(light_chain)

        all_human_chains = antibody_chains + (antigen_chains if antigen_chains else [antigen_chain])

        # Clean human complex (keep antigen + antibody, remove water/hetatm)
        cleaned_human_cif = os.path.join(self.sabdab_cleaned_dir, f"{pdb_id}.cif")
        if human_cif and not os.path.exists(cleaned_human_cif):
            clean_structure(human_cif, cleaned_human_cif, all_human_chains)

        # === Phase 3: Split antigen and antibody ===
        # Save antibody to SAbDab/antibody/
        antibody_cif = os.path.join(self.sabdab_antibody_dir, f"{pdb_id}_antibody.cif")
        if human_cif and not os.path.exists(antibody_cif):
            save_cif_with_metadata(human_cif, antibody_cif, chain_ids=antibody_chains)

        # Save human antigen to SAbDab/antigen/
        antigen_cif = os.path.join(self.sabdab_antigen_dir, f"{pdb_id}_antigen.cif")
        if human_cif and not os.path.exists(antigen_cif):
            ag_chains = antigen_chains if antigen_chains else [antigen_chain]
            save_cif_with_metadata(human_cif, antigen_cif, chain_ids=ag_chains)

        # === Phase 4: Human-Mouse alignment ===
        # Create pair directory: HumanMouse/{human_gene}_{mouse_gene}/
        # Sanitize directory name - remove special characters
        safe_human_gene = "".join(c if c.isalnum() or c in '-_' else '_' for c in (gene_name or 'Human'))
        safe_human_gene = safe_human_gene[:50]
        safe_mouse_gene = "".join(c if c.isalnum() or c in '-_' else '_' for c in (mouse_gene_name or 'Mouse'))
        safe_mouse_gene = safe_mouse_gene[:50]
        pair_name = f"{safe_human_gene}_{safe_mouse_gene}"
        pair_dir = os.path.join(self.human_mouse_dir, pair_name)
        os.makedirs(pair_dir, exist_ok=True)

        # Save human antigen (ALL chains) to pair directory
        human_ag_cif = os.path.join(pair_dir, "human_ag.cif")
        all_antigen_chains = antigen_chains if antigen_chains else [antigen_chain]
        if human_cif:
            save_cif_with_metadata(human_cif, human_ag_cif, chain_ids=all_antigen_chains)
        else:
            save_structure_cif(human_structure, human_ag_cif, chain_ids=all_antigen_chains)

        # Prepare for alignment (use first chain for alignment reference)
        human_ag_pdb = os.path.join(pair_dir, "human_ag_ref.pdb")
        save_structure(human_structure, human_ag_pdb, chain_ids=[antigen_chain])

        mouse_ag_pdb = os.path.join(pair_dir, "mouse_ag_unaligned.pdb")
        save_structure(mouse_structure, mouse_ag_pdb, chain_ids=[mouse_chain])

        # Perform alignment
        rotran = None
        aligned_structure = mouse_structure
        if self.use_pymol:
            rmsd = align_structures_pymol(
                mobile_file=mouse_ag_pdb,
                reference_file=human_ag_pdb,
                output_file=mouse_ag_pdb,
                mobile_chain=mouse_chain,
                reference_chain=antigen_chain
            )
            aligned_structure = parse_structure(mouse_ag_pdb, mouse_pdb_id) or mouse_structure
        else:
            rmsd, aligned_structure, rotran = align_structures_biopython(
                mobile_structure=mouse_structure,
                reference_structure=human_structure,
                mobile_chain=mouse_chain,
                reference_chain=antigen_chain
            )

        # Save aligned mouse structure
        mouse_aligned_cif = os.path.join(pair_dir, "mouse_ag_aligned.cif")
        mouse_aligned_pdb = os.path.join(pair_dir, "mouse_ag_aligned.pdb")

        if rotran and mouse_cif:
            save_cif_with_metadata(mouse_cif, mouse_aligned_cif, chain_ids=[mouse_chain], rotran=rotran)
        else:
            save_structure_cif(aligned_structure, mouse_aligned_cif, chain_ids=[mouse_chain])

        save_structure(aligned_structure, mouse_aligned_pdb, chain_ids=[mouse_chain])

        # Clean up temporary files
        if os.path.exists(mouse_ag_pdb):
            os.remove(mouse_ag_pdb)
        if os.path.exists(human_ag_pdb):
            os.remove(human_ag_pdb)

        # Save metadata for this pair
        pair_metadata = {
            'pair_name': pair_name,
            'pdb_id_human': pdb_id,
            'pdb_id_mouse': mouse_pdb_id,
            'human_antigen_chains': all_antigen_chains,
            'mouse_antigen_chain': mouse_chain,
            'antibody_heavy_chain': heavy_chain,
            'antibody_light_chain': light_chain,
            'human_uniprot': human_uniprot,
            'mouse_uniprot': mouse_uniprot,
            'human_gene_name': gene_name,
            'mouse_gene_name': mouse_gene_name,
            'alignment_rmsd': rmsd,
            'mouse_structure_source': mouse_structure_source,
            'files': {
                'human_ag': 'human_ag.cif',
                'mouse_ag_aligned': 'mouse_ag_aligned.cif'
            }
        }

        with open(os.path.join(pair_dir, 'metadata.json'), 'w') as f:
            json.dump(pair_metadata, f, indent=2)

        return rmsd, pair_metadata

    def _build_mouse_antigen_library(self) -> Dict[str, str]:
        """
        Phase 5: Build mouse antigen library by grouping aligned mouse antigens by gene name.

        Groups all aligned mouse antigen chains by mouse_gene_name and saves them
        to MouseAntigen/{mouse_gene_name}.cif

        Returns:
            Dictionary mapping mouse_gene_name to output file path
        """
        self.log("Phase 5: Building mouse antigen library")

        # Group successful data points by mouse gene name
        gene_to_datapoints: Dict[str, List[DataPoint]] = {}
        for dp in self.data_points:
            if dp.status == "success" and dp.mouse_gene_name:
                if dp.mouse_gene_name not in gene_to_datapoints:
                    gene_to_datapoints[dp.mouse_gene_name] = []
                gene_to_datapoints[dp.mouse_gene_name].append(dp)

        self.log(f"Found {len(gene_to_datapoints)} unique mouse genes to process")

        library_files: Dict[str, str] = {}

        import shutil

        for gene_name, datapoints in gene_to_datapoints.items():
            try:
                # Find the best structure for this gene (lowest RMSD)
                best_dp = min(datapoints, key=lambda dp: dp.alignment_rmsd or float('inf'))

                # Get the aligned mouse structure from HumanMouse directory
                # Use same sanitization as in _process_structures
                safe_human_gene = "".join(c if c.isalnum() or c in '-_' else '_' for c in (best_dp.gene_name or 'Human'))[:50]
                safe_mouse_gene = "".join(c if c.isalnum() or c in '-_' else '_' for c in (gene_name or 'Mouse'))[:50]
                pair_name = f"{safe_human_gene}_{safe_mouse_gene}"
                pair_dir = os.path.join(self.human_mouse_dir, pair_name)
                aligned_mouse_cif = os.path.join(pair_dir, "mouse_ag_aligned.cif")

                if not os.path.exists(aligned_mouse_cif):
                    self.log(f"Warning: Aligned mouse structure not found for {gene_name} at {aligned_mouse_cif}", level="WARNING")
                    continue

                # Copy to MouseAntigen directory
                output_path = os.path.join(self.mouse_antigen_dir, f"{safe_mouse_gene}.cif")
                shutil.copy2(aligned_mouse_cif, output_path)

                library_files[gene_name] = output_path
                self.log(f"Added {gene_name} to mouse antigen library (from {best_dp.pdb_id_mouse})")

            except Exception as e:
                self.log(f"Error processing {gene_name}: {e}", level="ERROR")

        self.log(f"Mouse antigen library built with {len(library_files)} entries")
        return library_files

    def _update_pairs_csv(self):
        """
        Generate human_mouse_pairs.csv in data directory.

        This CSV provides a complete mapping of:
        - Human antigen (PDB, chains, gene, UniProt)
        - Mouse antigen (PDB, chain, gene, UniProt)
        - Antibody (PDB, heavy chain, light chain)
        """
        pairs_data = []
        for dp in self.data_points:
            if dp.status == "success":
                # Sanitize gene names for pair_name
                safe_human_gene = "".join(c if c.isalnum() or c in '-_' else '_' for c in (dp.gene_name or 'Human'))[:50]
                safe_mouse_gene = "".join(c if c.isalnum() or c in '-_' else '_' for c in (dp.mouse_gene_name or 'Mouse'))[:50]
                pair_name = f"{safe_human_gene}_{safe_mouse_gene}"

                # Format antigen chains as comma-separated string
                antigen_chains_str = ','.join(dp.human_antigen_chains) if dp.human_antigen_chains else dp.human_antigen_chain

                pairs_data.append({
                    # Pair identification
                    'pair_name': pair_name,
                    'pair_dir': f"HumanMouse/{pair_name}",

                    # Human antigen info
                    'human_antigen_pdb': dp.pdb_id_human,
                    'human_antigen_chains': antigen_chains_str,
                    'human_antigen_gene': dp.gene_name,
                    'human_antigen_uniprot': dp.human_uniprot,
                    'human_antigen_name': dp.protein_name,
                    'human_antigen_file': f"SAbDab/antigen/{dp.pdb_id_human}_antigen.cif",

                    # Mouse antigen info
                    'mouse_antigen_pdb': dp.pdb_id_mouse,
                    'mouse_antigen_chain': dp.mouse_antigen_chain,
                    'mouse_antigen_gene': dp.mouse_gene_name,
                    'mouse_antigen_uniprot': dp.mouse_uniprot,
                    'mouse_antigen_source': dp.mouse_structure_source,
                    'mouse_antigen_file': f"MouseAntigen/{safe_mouse_gene}.cif",

                    # Antibody info
                    'antibody_pdb': dp.pdb_id_human,
                    'antibody_heavy_chain': dp.heavy_chain,
                    'antibody_light_chain': dp.light_chain,
                    'antibody_file': f"SAbDab/antibody/{dp.pdb_id_human}_antibody.cif",

                    # Alignment metrics
                    'sequence_identity': dp.sequence_identity,
                    'alignment_rmsd': dp.alignment_rmsd,

                    # Epitope validation
                    'epitope_identity': dp.epitope_identity,
                    'epitope_rmsd': dp.epitope_rmsd,
                    'epitope_consistent': dp.epitope_consistent
                })

        if pairs_data:
            pairs_df = pd.DataFrame(pairs_data)
            csv_path = os.path.join(self.data_dir, 'human_mouse_pairs.csv')
            pairs_df.to_csv(csv_path, index=False)
            self.log(f"Saved {len(pairs_data)} pairs to {csv_path}")

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
