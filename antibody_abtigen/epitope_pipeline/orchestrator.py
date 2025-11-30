"""
Pipeline orchestrator for the epitope-centric pipeline.

This module implements the main orchestrator that coordinates all pipeline stages:
1. Structure Cleaning (cleaner.py)
2. Epitope Extraction (extractor.py)
3. ESM-2 Embedding (encoder.py)
4. Similarity Grouping (grouper.py)
5. Structure Alignment (aligner.py)

The orchestrator provides:
- Full pipeline execution with a single command
- Checkpointing for resume from failures
- Progress tracking and logging
- Summary statistics generation
"""

import json
import logging
import time
from pathlib import Path
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime

from .core import (
    PipelineOrchestrator,
    PipelineConfig,
    CleanedStructure,
    EpitopeResidues,
    EpitopePipelineError,
)
from .cleaner import GemmiStructureCleaner, FilterResult, save_filter_log
from .extractor import GeometricEpitopeExtractor
from .encoder import ESM2EpitopeEncoder, EncoderOutput
from .storage import HDF5EmbeddingStore
from .grouper import (
    NumpyEpitopeGrouper,
    GroupResult,
    GroupingOutput,
    save_groups_json,
    save_grouping_stats_csv,
)
from .aligner import (
    PyMOLStructureAligner,
    GroupAlignmentOutput,
    save_alignment_summary_csv,
)
from .epitope_log import (
    save_epitope_residues_csv,
    save_epitope_summary_csv,
    save_embedding_stats_csv,
)

logger = logging.getLogger(__name__)


@dataclass
class PipelineCheckpoint:
    """Checkpoint for resuming pipeline execution."""
    stage: str  # 'clean', 'embed', 'group', 'align', 'done'
    processed_pdb_ids: List[str]
    timestamp: str
    config_hash: str
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict) -> 'PipelineCheckpoint':
        return cls(**data)


@dataclass
class PipelineResult:
    """Result of running the full pipeline."""
    success: bool
    stages_completed: List[str]
    total_structures: int
    cleaned_structures: int
    embedded_structures: int
    groups_found: int
    groups_aligned: int
    total_time_seconds: float
    output_dir: Path
    errors: List[str] = field(default_factory=list)
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        result = asdict(self)
        result['output_dir'] = str(self.output_dir)
        return result


class EpitopePipeline(PipelineOrchestrator):
    """
    Main orchestrator for the epitope-centric pipeline.

    Coordinates all stages of the pipeline:
    1. Clean: Filter and clean CIF structures
    2. Embed: Generate ESM-2 embeddings for epitopes
    3. Group: Cluster epitopes by embedding similarity
    4. Align: Align structures within groups

    Usage:
        >>> config = create_default_config(data_dir="./data")
        >>> pipeline = EpitopePipeline(config)
        >>> result = pipeline.run_full(input_dir, output_dir)
    """

    def __init__(
        self,
        config: Optional[PipelineConfig] = None,
        sabdab_summary_path: Optional[Path] = None,
        verbose: bool = True
    ):
        """
        Initialize pipeline.

        Args:
            config: Pipeline configuration (optional, uses defaults if None)
            sabdab_summary_path: Path to SAbDab summary TSV
            verbose: Print progress messages
        """
        self.config = config
        self.sabdab_summary_path = sabdab_summary_path
        self.verbose = verbose

        # Initialize components lazily
        self._cleaner = None
        self._extractor = None
        self._encoder = None
        self._store = None
        self._grouper = None
        self._aligner = None

    def _log(self, message: str):
        """Log message if verbose."""
        if self.verbose:
            print(message)
        logger.info(message)

    def _init_cleaner(self) -> GemmiStructureCleaner:
        """Initialize structure cleaner."""
        if self._cleaner is None:
            self._cleaner = GemmiStructureCleaner(
                sabdab_summary_path=self.sabdab_summary_path
            )
        return self._cleaner

    def _init_extractor(self) -> GeometricEpitopeExtractor:
        """Initialize epitope extractor."""
        if self._extractor is None:
            distance = 5.0
            if self.config:
                distance = self.config.contact_distance_threshold
            self._extractor = GeometricEpitopeExtractor(distance_threshold=distance)
        return self._extractor

    def _init_encoder(self) -> ESM2EpitopeEncoder:
        """Initialize ESM-2 encoder."""
        if self._encoder is None:
            device = "cuda"
            use_fp16 = True
            cache_dir = None
            if self.config:
                device = self.config.device
                use_fp16 = self.config.use_fp16
                cache_dir = self.config.esm_cache_dir
            self._encoder = ESM2EpitopeEncoder(
                device=device,
                use_fp16=use_fp16,
                cache_dir=cache_dir
            )
        return self._encoder

    def _init_store(self) -> HDF5EmbeddingStore:
        """Initialize embedding store."""
        if self._store is None:
            self._store = HDF5EmbeddingStore()
        return self._store

    def _init_grouper(self) -> NumpyEpitopeGrouper:
        """Initialize epitope grouper."""
        if self._grouper is None:
            threshold = 0.85
            if self.config:
                threshold = self.config.similarity_threshold
            self._grouper = NumpyEpitopeGrouper(
                similarity_threshold=threshold,
                min_group_size=2,
                exclude_same_pdb=True
            )
        return self._grouper

    def _init_aligner(self) -> PyMOLStructureAligner:
        """Initialize structure aligner."""
        if self._aligner is None:
            use_super = True
            if self.config:
                use_super = self.config.use_pymol_super
            self._aligner = PyMOLStructureAligner(use_super=use_super)
        return self._aligner

    def run(
        self,
        input_cif_paths: List[Path],
        output_json: Path,
        checkpoint_every: int = 100
    ) -> Dict:
        """
        Run the full pipeline (interface method).

        Args:
            input_cif_paths: List of CIF files to process
            output_json: Where to save final results
            checkpoint_every: Checkpoint frequency

        Returns:
            Summary dictionary
        """
        output_dir = output_json.parent
        result = self.run_full(
            input_dir=None,
            output_dir=output_dir,
            cif_files=input_cif_paths,
            checkpoint_every=checkpoint_every
        )
        return result.to_dict()

    def run_full(
        self,
        input_dir: Optional[Path],
        output_dir: Path,
        cif_files: Optional[List[Path]] = None,
        limit: Optional[int] = None,
        checkpoint_every: int = 100,
        skip_stages: Optional[List[str]] = None,
        resume_from: Optional[Path] = None
    ) -> PipelineResult:
        """
        Run the full 4-stage pipeline.

        Args:
            input_dir: Directory containing raw CIF files
            output_dir: Output directory for all results
            cif_files: Optional explicit list of CIF files (overrides input_dir)
            limit: Maximum structures to process
            checkpoint_every: Save checkpoint after N structures
            skip_stages: List of stages to skip ('clean', 'embed', 'group', 'align')
            resume_from: Path to checkpoint file to resume from

        Returns:
            PipelineResult with statistics
        """
        start_time = time.time()
        skip_stages = skip_stages or []
        errors = []

        # Setup output directories
        output_dir = Path(output_dir)
        cleaned_dir = output_dir / "cleaned"
        embeddings_dir = output_dir / "embeddings"
        grouping_dir = output_dir / "grouping"
        aligned_dir = output_dir / "aligned"

        for d in [cleaned_dir, embeddings_dir, grouping_dir, aligned_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Get CIF files
        if cif_files:
            all_cif_files = [Path(f) for f in cif_files]
        elif input_dir:
            input_dir = Path(input_dir)
            all_cif_files = sorted(input_dir.glob("*.cif"))
        else:
            raise ValueError("Either input_dir or cif_files must be provided")

        if limit:
            all_cif_files = all_cif_files[:limit]

        self._log(f"\n{'='*70}")
        self._log("Epitope-Centric Pipeline")
        self._log(f"{'='*70}")
        self._log(f"Input: {len(all_cif_files)} CIF files")
        self._log(f"Output: {output_dir}")
        self._log(f"Stages: clean → embed → group → align")
        self._log("")

        stages_completed = []
        structures: Dict[str, CleanedStructure] = {}
        epitopes: Dict[str, EpitopeResidues] = {}
        encoder_outputs: List[EncoderOutput] = []
        grouping_output: Optional[GroupingOutput] = None
        alignment_outputs: List[GroupAlignmentOutput] = []

        # Stage 1: Clean
        if 'clean' not in skip_stages:
            self._log("Stage 1/4: Cleaning structures...")
            try:
                structures, filter_results = self._run_clean_stage(
                    all_cif_files, cleaned_dir
                )
                stages_completed.append('clean')
                self._log(f"  Cleaned: {len(structures)}/{len(all_cif_files)} structures")
            except Exception as e:
                errors.append(f"Clean stage failed: {e}")
                logger.exception("Clean stage failed")
        else:
            # Load existing cleaned structures - only load the ones in all_cif_files
            self._log("Stage 1/4: Loading pre-cleaned structures...")
            structures = self._load_cleaned_structures_from_list(all_cif_files)
            self._log(f"  Loaded: {len(structures)} structures")

        if not structures:
            return PipelineResult(
                success=False,
                stages_completed=stages_completed,
                total_structures=len(all_cif_files),
                cleaned_structures=0,
                embedded_structures=0,
                groups_found=0,
                groups_aligned=0,
                total_time_seconds=time.time() - start_time,
                output_dir=output_dir,
                errors=errors + ["No structures to process"]
            )

        # Stage 2: Embed
        if 'embed' not in skip_stages:
            self._log("\nStage 2/4: Generating embeddings...")
            try:
                epitopes, encoder_outputs = self._run_embed_stage(
                    structures, embeddings_dir
                )
                stages_completed.append('embed')
                self._log(f"  Embedded: {len(encoder_outputs)} structures")
            except Exception as e:
                errors.append(f"Embed stage failed: {e}")
                logger.exception("Embed stage failed")
        else:
            # Load existing embeddings
            self._log("\nStage 2/4: Loading pre-computed embeddings...")
            encoder_outputs = self._load_embeddings(embeddings_dir)
            # Re-extract epitopes for alignment
            epitopes = self._extract_epitopes(structures)
            self._log(f"  Loaded: {len(encoder_outputs)} embeddings")

        if not encoder_outputs:
            return PipelineResult(
                success=False,
                stages_completed=stages_completed,
                total_structures=len(all_cif_files),
                cleaned_structures=len(structures),
                embedded_structures=0,
                groups_found=0,
                groups_aligned=0,
                total_time_seconds=time.time() - start_time,
                output_dir=output_dir,
                errors=errors + ["No embeddings generated"]
            )

        # Stage 3: Group
        if 'group' not in skip_stages:
            self._log("\nStage 3/4: Grouping epitopes...")
            try:
                grouping_output = self._run_group_stage(
                    encoder_outputs, grouping_dir
                )
                stages_completed.append('group')
                self._log(f"  Groups: {len(grouping_output.groups)}")
                self._log(f"  Grouped: {grouping_output.grouped_epitopes}/{grouping_output.total_epitopes}")
            except Exception as e:
                errors.append(f"Group stage failed: {e}")
                logger.exception("Group stage failed")
        else:
            # Load existing groups
            self._log("\nStage 3/4: Loading pre-computed groups...")
            grouping_output = self._load_groups(grouping_dir)
            if grouping_output:
                self._log(f"  Loaded: {len(grouping_output.groups)} groups")

        if not grouping_output or not grouping_output.groups:
            return PipelineResult(
                success=len(stages_completed) > 0,
                stages_completed=stages_completed,
                total_structures=len(all_cif_files),
                cleaned_structures=len(structures),
                embedded_structures=len(encoder_outputs),
                groups_found=0,
                groups_aligned=0,
                total_time_seconds=time.time() - start_time,
                output_dir=output_dir,
                errors=errors + ["No groups found"]
            )

        # Stage 4: Align
        if 'align' not in skip_stages:
            self._log("\nStage 4/4: Aligning structures...")
            try:
                alignment_outputs = self._run_align_stage(
                    grouping_output.groups, epitopes, structures, aligned_dir
                )
                stages_completed.append('align')
                self._log(f"  Aligned: {len(alignment_outputs)} groups")
            except Exception as e:
                errors.append(f"Align stage failed: {e}")
                logger.exception("Align stage failed")

        # Save final summary
        total_time = time.time() - start_time
        result = PipelineResult(
            success=len(errors) == 0,
            stages_completed=stages_completed,
            total_structures=len(all_cif_files),
            cleaned_structures=len(structures),
            embedded_structures=len(encoder_outputs),
            groups_found=len(grouping_output.groups) if grouping_output else 0,
            groups_aligned=len(alignment_outputs),
            total_time_seconds=total_time,
            output_dir=output_dir,
            errors=errors,
            metadata={
                'timestamp': datetime.now().isoformat(),
                'config': asdict(self.config) if self.config else None,
            }
        )

        # Save summary JSON
        summary_path = output_dir / "pipeline_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(result.to_dict(), f, indent=2, default=str)

        self._log(f"\n{'='*70}")
        self._log("Pipeline Complete!")
        self._log(f"{'='*70}")
        self._log(f"Total time: {total_time:.1f}s")
        self._log(f"Structures: {len(all_cif_files)} → {len(structures)} cleaned → {len(encoder_outputs)} embedded")
        self._log(f"Groups: {len(grouping_output.groups) if grouping_output else 0} found, {len(alignment_outputs)} aligned")
        self._log(f"Output: {output_dir}")

        if errors:
            self._log(f"\nWarnings: {len(errors)}")
            for err in errors[:5]:
                self._log(f"  - {err}")

        return result

    def _run_clean_stage(
        self,
        cif_files: List[Path],
        output_dir: Path
    ) -> tuple:
        """Run structure cleaning stage."""
        cleaner = self._init_cleaner()
        cleaner.output_dir = output_dir

        structures = {}
        filter_results = []

        for cif_path in cif_files:
            try:
                cleaned, filter_result = cleaner.clean_structure(cif_path)
                filter_results.append(filter_result)
                if cleaned:
                    structures[cleaned.pdb_id] = cleaned
            except Exception as e:
                logger.warning(f"Failed to clean {cif_path.name}: {e}")

        # Save filter log
        save_filter_log(filter_results, output_dir / "filtering_log.csv")

        return structures, filter_results

    def _run_embed_stage(
        self,
        structures: Dict[str, CleanedStructure],
        output_dir: Path
    ) -> tuple:
        """Run embedding generation stage."""
        extractor = self._init_extractor()
        encoder = self._init_encoder()
        store = self._init_store()

        epitopes = {}
        encoder_outputs = []

        for pdb_id, structure in structures.items():
            try:
                # Extract epitope
                epitope = extractor.extract_epitope(structure)
                epitopes[epitope.epitope_id] = epitope

                # Encode
                output = encoder.encode_full(epitope, structure)
                encoder_outputs.append(output)
            except Exception as e:
                logger.warning(f"Failed to embed {pdb_id}: {e}")

        # Save embeddings
        if encoder_outputs:
            h5_path = output_dir / "embeddings.h5"
            store.save_encoder_outputs(encoder_outputs, h5_path)

            # Save CSV logs
            save_epitope_residues_csv(
                list(epitopes.values()),
                structures,
                output_dir / "epitope_residues.csv"
            )
            save_epitope_summary_csv(
                list(epitopes.values()),
                structures,
                output_dir / "epitope_summary.csv",
                encoder_outputs
            )
            save_embedding_stats_csv(encoder_outputs, output_dir / "embedding_stats.csv")

        return epitopes, encoder_outputs

    def _run_group_stage(
        self,
        encoder_outputs: List[EncoderOutput],
        output_dir: Path
    ) -> GroupingOutput:
        """Run epitope grouping stage."""
        grouper = self._init_grouper()

        # Build index
        grouper.build_index_from_encoder_outputs(encoder_outputs, use_epitope_embedding=True)

        # Compute similarity and find groups
        grouper.compute_similarity_matrix()
        grouping_output = grouper.find_groups_detailed()

        # Save outputs
        save_groups_json(grouping_output, output_dir / "groups.json")
        save_grouping_stats_csv(grouping_output, output_dir / "grouping_stats.csv")
        grouper.save_sparse_similarity(output_dir / "similarity_sparse.h5")

        return grouping_output

    def _run_align_stage(
        self,
        groups: List[GroupResult],
        epitopes: Dict[str, EpitopeResidues],
        structures: Dict[str, CleanedStructure],
        output_dir: Path
    ) -> List[GroupAlignmentOutput]:
        """Run structure alignment stage."""
        aligner = self._init_aligner()

        alignment_outputs = []

        for group in groups:
            try:
                output = aligner.align_group_detailed(
                    group, epitopes, structures, output_dir
                )
                alignment_outputs.append(output)
            except Exception as e:
                logger.warning(f"Failed to align {group.group_id}: {e}")

        # Save summary
        if alignment_outputs:
            save_alignment_summary_csv(alignment_outputs, output_dir / "alignment_summary.csv")

        return alignment_outputs

    def _load_cleaned_structures_from_list(
        self,
        cif_files: List[Path]
    ) -> Dict[str, CleanedStructure]:
        """Load pre-cleaned structures from a list of CIF files."""
        cleaner = self._init_cleaner()
        structures = {}

        for cif_path in cif_files:
            try:
                cleaned, _ = cleaner.clean_structure(cif_path, skip_filtering=True)
                if cleaned:
                    structures[cleaned.pdb_id] = cleaned
            except Exception as e:
                logger.warning(f"Failed to load {cif_path.name}: {e}")

        return structures

    def _load_cleaned_structures(self, cleaned_dir: Path) -> Dict[str, CleanedStructure]:
        """Load pre-cleaned structures from directory."""
        cleaner = self._init_cleaner()
        structures = {}

        for cif_path in cleaned_dir.glob("*_cleaned.cif"):
            try:
                cleaned, _ = cleaner.clean_structure(cif_path, skip_filtering=True)
                if cleaned:
                    structures[cleaned.pdb_id] = cleaned
            except Exception as e:
                logger.warning(f"Failed to load {cif_path.name}: {e}")

        return structures

    def _load_embeddings(self, embeddings_dir: Path) -> List[EncoderOutput]:
        """Load pre-computed embeddings."""
        store = self._init_store()
        h5_path = embeddings_dir / "embeddings.h5"

        if h5_path.exists():
            return store.load_all_encoder_outputs(h5_path)
        return []

    def _extract_epitopes(
        self,
        structures: Dict[str, CleanedStructure]
    ) -> Dict[str, EpitopeResidues]:
        """Extract epitopes from structures."""
        extractor = self._init_extractor()
        epitopes = {}

        for pdb_id, structure in structures.items():
            try:
                epitope = extractor.extract_epitope(structure)
                epitopes[epitope.epitope_id] = epitope
            except Exception as e:
                logger.warning(f"Failed to extract epitope from {pdb_id}: {e}")

        return epitopes

    def _load_groups(self, grouping_dir: Path) -> Optional[GroupingOutput]:
        """Load pre-computed groups."""
        groups_path = grouping_dir / "groups.json"

        if not groups_path.exists():
            return None

        with open(groups_path, 'r') as f:
            data = json.load(f)

        from .grouper import GroupMember

        groups = []
        for g in data.get('groups', []):
            members = [
                GroupMember(
                    epitope_id=m['epitope_id'],
                    pdb_id=m['pdb_id'],
                    is_reference=m['is_reference'],
                    antigen_chains=m.get('antigen_chains', []),
                    epitope_residues=m.get('epitope_residues', {}),
                    total_residues=m.get('total_residues', 0),
                    similarity_to_ref=m.get('similarity_to_ref')
                )
                for m in g.get('members', [])
            ]
            groups.append(GroupResult(
                group_id=g['group_id'],
                reference_epitope_id=g['reference_epitope_id'],
                member_count=g['member_count'],
                avg_similarity=g['avg_similarity'],
                min_similarity=g.get('min_similarity', 0),
                max_similarity=g.get('max_similarity', 1),
                members=members
            ))

        metadata = data.get('metadata', {})
        return GroupingOutput(
            groups=groups,
            total_epitopes=metadata.get('total_epitopes', 0),
            grouped_epitopes=metadata.get('grouped_epitopes', 0),
            singleton_count=metadata.get('singleton_count', 0),
            similarity_threshold=metadata.get('similarity_threshold', 0.85)
        )


def run_pipeline(
    input_dir: Path,
    output_dir: Path,
    sabdab_summary: Optional[Path] = None,
    limit: Optional[int] = None,
    device: str = "cuda",
    similarity_threshold: float = 0.85,
    distance_threshold: float = 5.0,
    verbose: bool = True
) -> PipelineResult:
    """
    Convenience function to run the full pipeline.

    Args:
        input_dir: Directory containing raw CIF files
        output_dir: Output directory
        sabdab_summary: Path to SAbDab summary TSV
        limit: Maximum structures to process
        device: Device for ESM-2 ('cuda' or 'cpu')
        similarity_threshold: Grouping threshold
        distance_threshold: Epitope extraction distance
        verbose: Print progress

    Returns:
        PipelineResult
    """
    from .core import create_default_config

    config = create_default_config(
        data_dir=output_dir,
        device=device,
        similarity_threshold=similarity_threshold
    )
    config.contact_distance_threshold = distance_threshold

    pipeline = EpitopePipeline(
        config=config,
        sabdab_summary_path=sabdab_summary,
        verbose=verbose
    )

    return pipeline.run_full(
        input_dir=input_dir,
        output_dir=output_dir,
        limit=limit
    )
