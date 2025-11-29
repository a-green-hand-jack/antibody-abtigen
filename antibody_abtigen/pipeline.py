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
    align_structures_pymol,
    align_structures_biopython,
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
        antigen_chain = str(row['antigen_chain']).split('|')[0].strip()  # Take first chain if multiple
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
        save_structure_cif(human_structure, antibody_cif, chain_ids=antibody_chains)

        # Save human antigen
        human_ag_pdb = os.path.join(dp_dir, f"{dp_id}_human_ag.pdb")
        human_ag_cif = os.path.join(dp_dir, f"{dp_id}_human_ag.cif")
        save_structure(human_structure, human_ag_pdb, chain_ids=[antigen_chain])
        save_structure_cif(human_structure, human_ag_cif, chain_ids=[antigen_chain])

        # Align mouse antigen to human antigen position
        mouse_ag_pdb = os.path.join(dp_dir, f"{dp_id}_mouse_ag.pdb")
        mouse_ag_cif = os.path.join(dp_dir, f"{dp_id}_mouse_ag.cif")

        # First save unaligned mouse structure
        save_structure(mouse_structure, mouse_ag_pdb, chain_ids=[mouse_chain])

        # Perform alignment
        if self.use_pymol:
            # Use PyMOL for alignment
            rmsd = align_structures_pymol(
                mobile_file=mouse_ag_pdb,
                reference_file=human_ag_pdb,
                output_file=mouse_ag_pdb,
                mobile_chain=mouse_chain,
                reference_chain=antigen_chain
            )
        else:
            # Use Biopython for alignment
            rmsd, aligned_structure = align_structures_biopython(
                mobile_structure=mouse_structure,
                reference_structure=human_structure,
                mobile_chain=mouse_chain,
                reference_chain=antigen_chain
            )
            save_structure(aligned_structure, mouse_ag_pdb, chain_ids=[mouse_chain])

        # Save aligned mouse structure as CIF too
        aligned_mouse = parse_structure(mouse_ag_pdb, "aligned_mouse")
        if aligned_mouse:
            save_structure_cif(aligned_mouse, mouse_ag_cif)

        # Save metadata
        metadata = {
            'id': dp_id,
            'pdb_id_human': pdb_id,
            'pdb_id_mouse': mouse_pdb_id,
            'human_antigen_chain': antigen_chain,
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
                'mouse_ag_pdb': os.path.basename(mouse_ag_pdb),
                'mouse_ag_cif': os.path.basename(mouse_ag_cif)
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
