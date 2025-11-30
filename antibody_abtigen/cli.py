"""
Command-line interface for antibody-abtigen.
"""

import os
import sys
import click

from .pipeline import CrossSpeciesDatasetPipeline
from .filtering import filter_dataset
from .epitope_pipeline import (
    GemmiStructureCleaner,
    GeometricEpitopeExtractor,
    FilterResult,
    save_filter_log,
    generate_filter_summary,
)


@click.group()
@click.version_option()
def cli():
    """Cross-Species Antibody-Antigen Structure Dataset Builder."""
    pass


@cli.command('build')
@click.option(
    '--output', '-o',
    default='./output',
    help='Output directory for the dataset',
    type=click.Path()
)
@click.option(
    '--data-dir', '-d',
    default='./data',
    help='Directory for intermediate data (downloads, cache)',
    type=click.Path()
)
@click.option(
    '--limit', '-n',
    default=None,
    type=int,
    help='Maximum number of entries to process (None = all)'
)
@click.option(
    '--resolution', '-r',
    default=2.5,
    type=float,
    help='Maximum resolution threshold in Angstroms (default: 2.5)'
)
@click.option(
    '--identity', '-i',
    default=50.0,
    type=float,
    help='Minimum sequence identity threshold in percent (default: 50.0)'
)
@click.option(
    '--dry-run',
    is_flag=True,
    default=False,
    help='Only analyze data without downloading structures'
)
@click.option(
    '--no-pymol',
    is_flag=True,
    default=False,
    help='Use Biopython instead of PyMOL for alignment'
)
@click.option(
    '--force-download',
    is_flag=True,
    default=False,
    help='Force re-download of SAbDab summary file'
)
@click.option(
    '--epitope-rmsd',
    default=1.5,
    type=float,
    help='Maximum epitope RMSD threshold in Angstroms (default: 1.5)'
)
@click.option(
    '--epitope-identity',
    default=80.0,
    type=float,
    help='Minimum epitope sequence identity threshold in percent (default: 80.0)'
)
@click.option(
    '--no-sabdab-mouse',
    is_flag=True,
    default=False,
    help='Do not prefer SAbDab mouse structures over RCSB'
)
def build_command(
    output: str,
    data_dir: str,
    limit: int,
    resolution: float,
    identity: float,
    dry_run: bool,
    no_pymol: bool,
    force_download: bool,
    epitope_rmsd: float,
    epitope_identity: float,
    no_sabdab_mouse: bool
):
    """
    Build cross-species antibody-antigen structure dataset.

    Downloads antibody-antigen complex data from SAbDab,
    finds mouse orthologs for human antigens, and generates aligned
    structure triplets for machine learning and computational biology.
    """
    click.echo("=" * 60)
    click.echo("Cross-Species Antibody-Antigen Structure Dataset Builder")
    click.echo("=" * 60)
    click.echo()

    click.echo("Configuration:")
    click.echo(f"  Output directory: {output}")
    click.echo(f"  Data directory: {data_dir}")
    click.echo(f"  Entry limit: {limit if limit else 'All'}")
    click.echo(f"  Resolution threshold: {resolution} A")
    click.echo(f"  Sequence identity threshold: {identity}%")
    click.echo(f"  Epitope RMSD threshold: {epitope_rmsd} A")
    click.echo(f"  Epitope identity threshold: {epitope_identity}%")
    click.echo(f"  Prefer SAbDab mouse: {not no_sabdab_mouse}")
    click.echo(f"  Dry run: {dry_run}")
    click.echo(f"  Use PyMOL: {not no_pymol}")
    click.echo()

    output = os.path.abspath(output)
    data_dir = os.path.abspath(data_dir)

    os.makedirs(output, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    pipeline = CrossSpeciesDatasetPipeline(
        data_dir=data_dir,
        output_dir=output,
        resolution_threshold=resolution,
        sequence_identity_threshold=identity,
        use_pymol=not no_pymol,
        epitope_rmsd_threshold=epitope_rmsd,
        epitope_identity_threshold=epitope_identity,
        prefer_sabdab_mouse=not no_sabdab_mouse
    )

    if force_download:
        summary_file = os.path.join(data_dir, "sabdab_summary.tsv")
        if os.path.exists(summary_file):
            os.remove(summary_file)
            click.echo("Removed existing SAbDab summary file")

    try:
        result_df = pipeline.run(limit=limit, dry_run=dry_run)

        click.echo()
        click.echo("=" * 60)
        click.echo("Pipeline Complete!")
        click.echo("=" * 60)

        if not result_df.empty:
            total = len(result_df)
            successful = len(result_df[result_df['status'] == 'success'])
            dry_run_count = len(result_df[result_df['status'] == 'dry_run'])
            failed = len(result_df[result_df['status'] == 'failed'])

            click.echo(f"\nResults:")
            click.echo(f"  Total processed: {total}")
            if dry_run:
                click.echo(f"  Candidates found: {dry_run_count}")
            else:
                click.echo(f"  Successful: {successful}")
            click.echo(f"  Failed: {failed}")

            click.echo(f"\nOutput files:")
            click.echo(f"  Summary CSV: {os.path.join(output, 'dataset_summary.csv')}")
            click.echo(f"  Processing log: {os.path.join(output, 'processing_log.json')}")

            if not dry_run and successful > 0:
                click.echo(f"  Data point folders: {successful} folders in {output}/")

    except KeyboardInterrupt:
        click.echo("\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        click.echo(f"\nError: {e}", err=True)
        raise


@cli.command('to-yaml')
@click.option(
    '--input', '-i', 'input_dir',
    required=True,
    help='Input directory containing DP_* folders with CIF files',
    type=click.Path(exists=True)
)
@click.option(
    '--output', '-o',
    required=True,
    help='Output directory for YAML files',
    type=click.Path()
)
@click.option(
    '--no-bonds',
    is_flag=True,
    default=False,
    help='Do not include disulfide bonds in YAML'
)
def to_yaml_command(input_dir: str, output: str, no_bonds: bool):
    """
    Convert CIF files to Boltz-1 YAML format.

    Takes the raw data directory containing DP_* folders and converts
    each CIF file to a corresponding YAML configuration file.

    \b
    Input structure:
        input_dir/
        └── DP_XXXX_Y/
            ├── DP_XXXX_Y_antibody.cif
            ├── DP_XXXX_Y_human_ag.cif
            └── DP_XXXX_Y_mouse_ag.cif

    \b
    Output structure:
        output/
        ├── DP_XXXX_Y_antibody.yaml
        ├── DP_XXXX_Y_human_ag.yaml
        ├── DP_XXXX_Y_mouse_ag.yaml
        └── conversion_summary.json
    """
    from .yaml_converter import batch_convert_to_yamls

    click.echo("=" * 60)
    click.echo("CIF to Boltz YAML Converter")
    click.echo("=" * 60)
    click.echo()

    click.echo("Configuration:")
    click.echo(f"  Input directory: {input_dir}")
    click.echo(f"  Output directory: {output}")
    click.echo(f"  Include bonds: {not no_bonds}")
    click.echo()

    try:
        results = batch_convert_to_yamls(
            raw_data_dir=input_dir,
            output_dir=output,
            include_bonds=not no_bonds,
            verbose=True
        )

        click.echo()
        click.echo("=" * 60)
        click.echo("Conversion Complete!")
        click.echo("=" * 60)

    except ImportError as e:
        click.echo(f"\nError: {e}", err=True)
        click.echo("Install missing dependencies: pip install gemmi pyyaml", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"\nError: {e}", err=True)
        raise


@cli.command('triplet')
@click.option(
    '--data-dir', '-d',
    default='./data',
    help='Data directory containing SAbDab/, MouseAntigen/, human_mouse_pairs.csv',
    type=click.Path(exists=True)
)
@click.option(
    '--output', '-o',
    default=None,
    help='Output directory (default: data_dir/HumanMouseAntigenAntibody)',
    type=click.Path()
)
@click.option(
    '--dry-run',
    is_flag=True,
    default=False,
    help='Only print what would be done, do not create files'
)
def triplet_command(data_dir: str, output: str, dry_run: bool):
    """
    Generate aligned triplets (human_antigen, mouse_antigen, antibody).

    Reads human_mouse_pairs.csv and creates aligned structure triplets
    using PyMOL. Each triplet contains:
    - human_antigen.cif (reference position)
    - mouse_antigen.cif (aligned to human)
    - antibody.cif (keeps relative position to human antigen)

    \b
    Output structure:
        output/
        └── {human_gene}_{mouse_gene}_{pdb_id}/
            ├── human_antigen.cif
            ├── mouse_antigen.cif
            └── antibody.cif
    """
    from .triplet import generate_triplets

    click.echo("=" * 60)
    click.echo("Triplet Alignment Generator")
    click.echo("=" * 60)
    click.echo()

    click.echo("Configuration:")
    click.echo(f"  Data directory: {data_dir}")
    click.echo(f"  Output directory: {output or 'data_dir/HumanMouseAntigenAntibody'}")
    click.echo(f"  Dry run: {dry_run}")
    click.echo()

    try:
        count = generate_triplets(
            data_dir=data_dir,
            output_dir=output,
            dry_run=dry_run
        )

        click.echo()
        click.echo("=" * 60)
        click.echo("Triplet Generation Complete!")
        click.echo("=" * 60)
        click.echo(f"Generated {count} triplets")

    except FileNotFoundError as e:
        click.echo(f"\nError: {e}", err=True)
        click.echo("Make sure you have run 'build' first to generate human_mouse_pairs.csv", err=True)
        sys.exit(1)
    except RuntimeError as e:
        click.echo(f"\nError: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"\nError: {e}", err=True)
        raise


@cli.command('clean')
@click.option(
    '--input', '-i', 'input_dir',
    required=True,
    help='Input directory containing raw CIF files from SAbDab',
    type=click.Path(exists=True)
)
@click.option(
    '--output', '-o', 'output_dir',
    required=True,
    help='Output directory for cleaned CIF files',
    type=click.Path()
)
@click.option(
    '--sabdab-summary', '-s',
    default=None,
    help='Path to SAbDab summary TSV file (default: input_dir/../meta/sabdab_summary_all.tsv)',
    type=click.Path(exists=True)
)
@click.option(
    '--skip-filtering',
    is_flag=True,
    default=False,
    help='Skip antigen type filtering (process all structures)'
)
@click.option(
    '--extract-epitopes',
    is_flag=True,
    default=False,
    help='Also extract and log epitope statistics for accepted structures'
)
@click.option(
    '--distance-threshold',
    default=5.0,
    show_default=True,
    type=float,
    help='Distance threshold (Å) for epitope extraction'
)
def clean_command(
    input_dir: str,
    output_dir: str,
    sabdab_summary: str,
    skip_filtering: bool,
    extract_epitopes: bool,
    distance_threshold: float
):
    """
    Clean and filter CIF structures for the epitope pipeline.

    This command processes raw CIF files from SAbDab and:
    1. Filters structures by antigen type (protein only by default)
    2. Removes water and HETATM records
    3. Preserves only antigen and antibody chains
    4. Generates a filtering log (CSV)

    \b
    Input structure:
        input_dir/
        ├── 7K8T.cif
        ├── 8ABC.cif
        └── ...

    \b
    Output structure:
        output_dir/
        ├── 7K8T_cleaned.cif
        ├── 8ABC_cleaned.cif
        ├── filtering_log.csv
        └── epitope_stats.csv (if --extract-epitopes)
    """
    from pathlib import Path
    import csv
    import statistics

    click.echo("=" * 70)
    click.echo("Structure Cleaning for Epitope Pipeline")
    click.echo("=" * 70)
    click.echo()

    input_path = Path(input_dir).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    # Determine SAbDab summary path
    if sabdab_summary:
        sabdab_path = Path(sabdab_summary)
    else:
        # Default: look for meta/sabdab_summary_all.tsv relative to input
        sabdab_path = input_path.parent / "meta" / "sabdab_summary_all.tsv"
        if not sabdab_path.exists():
            sabdab_path = None

    click.echo("Configuration:")
    click.echo(f"  Input directory: {input_path}")
    click.echo(f"  Output directory: {output_path}")
    click.echo(f"  SAbDab summary: {sabdab_path if sabdab_path else 'Not found (using heuristics)'}")
    click.echo(f"  Skip filtering: {skip_filtering}")
    click.echo(f"  Extract epitopes: {extract_epitopes}")
    if extract_epitopes:
        click.echo(f"  Distance threshold: {distance_threshold} Å")
    click.echo()

    # Get CIF files
    cif_files = sorted(input_path.glob("*.cif"))
    if not cif_files:
        click.echo("Error: No CIF files found in input directory", err=True)
        sys.exit(1)

    click.echo(f"Found {len(cif_files)} CIF files")
    click.echo()

    # Initialize cleaner
    cleaner = GemmiStructureCleaner(
        output_dir=output_path,
        sabdab_summary_path=sabdab_path
    )

    # Initialize extractor if needed
    extractor = None
    if extract_epitopes:
        extractor = GeometricEpitopeExtractor(distance_threshold=distance_threshold)

    # Process structures
    filter_results = []
    cleaned_structures = []
    epitope_stats = []
    errors = []

    with click.progressbar(cif_files, label='Processing structures') as bar:
        for cif_path in bar:
            try:
                cleaned, filter_result = cleaner.clean_structure(
                    cif_path, skip_filtering=skip_filtering
                )
                filter_results.append(filter_result)

                if cleaned:
                    cleaned_structures.append(cleaned)

                    # Extract epitope if requested
                    if extractor:
                        try:
                            epitope = extractor.extract_epitope(cleaned)
                            epitope_stats.append({
                                'pdb_id': cleaned.pdb_id,
                                'epitope_residues': epitope.total_residue_count(),
                                'num_contacts': epitope.num_contacts,
                                'antigen_chains': len(epitope.antigen_chains)
                            })
                        except Exception as e:
                            errors.append(f"Epitope extraction failed for {cleaned.pdb_id}: {e}")

            except Exception as e:
                errors.append(f"Failed to process {cif_path.name}: {e}")

    click.echo()

    # Save filtering log
    log_path = output_path / "filtering_log.csv"
    save_filter_log(filter_results, log_path)
    click.echo(f"Filtering log saved to: {log_path}")

    # Save epitope stats if extracted
    if epitope_stats:
        epitope_log_path = output_path / "epitope_stats.csv"
        with open(epitope_log_path, 'w', newline='') as f:
            writer = csv.DictWriter(
                f, fieldnames=['pdb_id', 'epitope_residues', 'num_contacts', 'antigen_chains']
            )
            writer.writeheader()
            writer.writerows(epitope_stats)
        click.echo(f"Epitope stats saved to: {epitope_log_path}")

    # Print summary
    click.echo()
    click.echo(generate_filter_summary(filter_results))

    # Print epitope stats if available
    if epitope_stats:
        epitope_counts = [s['epitope_residues'] for s in epitope_stats]
        click.echo()
        click.echo("=" * 70)
        click.echo("Epitope Extraction Statistics")
        click.echo("=" * 70)
        click.echo(f"Structures with epitopes: {len(epitope_stats)}")
        click.echo(f"Average epitope size: {statistics.mean(epitope_counts):.1f} residues")
        click.echo(f"Median epitope size: {statistics.median(epitope_counts):.1f} residues")
        click.echo(f"Range: {min(epitope_counts)} - {max(epitope_counts)} residues")
        click.echo("=" * 70)

    # Print errors if any
    if errors:
        click.echo()
        click.echo(f"Warnings ({len(errors)} issues):", err=True)
        for err in errors[:10]:  # Show first 10
            click.echo(f"  - {err}", err=True)
        if len(errors) > 10:
            click.echo(f"  ... and {len(errors) - 10} more", err=True)

    # Final summary
    click.echo()
    click.echo("=" * 70)
    click.echo("Cleaning Complete!")
    click.echo("=" * 70)
    click.echo(f"Cleaned structures: {len(cleaned_structures)}/{len(cif_files)}")
    click.echo(f"Output directory: {output_path}")


@cli.command('embed')
@click.option(
    '--input', '-i', 'input_dir',
    required=True,
    help='Input directory containing cleaned CIF files',
    type=click.Path(exists=True)
)
@click.option(
    '--output', '-o', 'output_dir',
    required=True,
    help='Output directory for embeddings and CSV logs',
    type=click.Path()
)
@click.option(
    '--sabdab-summary', '-s',
    default=None,
    help='Path to SAbDab summary TSV file',
    type=click.Path(exists=True)
)
@click.option(
    '--device',
    default='cuda',
    show_default=True,
    type=click.Choice(['cuda', 'cpu']),
    help='Device to run ESM-2 model on'
)
@click.option(
    '--no-fp16',
    is_flag=True,
    default=False,
    help='Disable FP16 mixed precision (use FP32)'
)
@click.option(
    '--distance-threshold',
    default=5.0,
    show_default=True,
    type=float,
    help='Distance threshold (Å) for epitope extraction'
)
@click.option(
    '--cache-dir',
    default=None,
    help='Directory to cache ESM-2 model weights',
    type=click.Path()
)
@click.option(
    '--limit', '-n',
    default=None,
    type=int,
    help='Maximum number of structures to process (None = all)'
)
def embed_command(
    input_dir: str,
    output_dir: str,
    sabdab_summary: str,
    device: str,
    no_fp16: bool,
    distance_threshold: float,
    cache_dir: str,
    limit: int
):
    """
    Generate ESM-2 embeddings for epitopes in cleaned CIF structures.

    This command:
    1. Extracts epitope residues from cleaned structures
    2. Encodes epitopes using ESM-2 3B model
    3. Saves embeddings to HDF5 and CSV logs

    \b
    Input structure:
        input_dir/
        ├── 1A14_cleaned.cif
        ├── 1A2Y_cleaned.cif
        └── ...

    \b
    Output structure:
        output_dir/
        ├── embeddings.h5              # HDF5 with full + epitope embeddings
        ├── epitope_residues.csv       # Per-residue epitope info
        ├── epitope_summary.csv        # Per-structure summary
        └── embedding_stats.csv        # Embedding statistics
    """
    from pathlib import Path

    from .epitope_pipeline import (
        GemmiStructureCleaner,
        GeometricEpitopeExtractor,
        ESM2EpitopeEncoder,
        HDF5EmbeddingStore,
        save_epitope_residues_csv,
        save_epitope_summary_csv,
        save_embedding_stats_csv,
        generate_epitope_report,
    )

    click.echo("=" * 70)
    click.echo("ESM-2 Epitope Embedding Generator")
    click.echo("=" * 70)
    click.echo()

    input_path = Path(input_dir).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    # Determine SAbDab summary path
    if sabdab_summary:
        sabdab_path = Path(sabdab_summary)
    else:
        sabdab_path = input_path.parent / "meta" / "sabdab_summary_all.tsv"
        if not sabdab_path.exists():
            sabdab_path = None

    # Cache directory
    cache_path = Path(cache_dir) if cache_dir else None

    click.echo("Configuration:")
    click.echo(f"  Input directory: {input_path}")
    click.echo(f"  Output directory: {output_path}")
    click.echo(f"  SAbDab summary: {sabdab_path if sabdab_path else 'Not found'}")
    click.echo(f"  Device: {device}")
    click.echo(f"  FP16: {not no_fp16}")
    click.echo(f"  Distance threshold: {distance_threshold} Å")
    click.echo(f"  Cache directory: {cache_path if cache_path else 'Default'}")
    click.echo(f"  Limit: {limit if limit else 'All'}")
    click.echo()

    # Get CIF files
    cif_files = [f for f in sorted(input_path.glob("*_cleaned.cif"))
                 if "_CLEANED_" not in f.name]
    if limit:
        cif_files = cif_files[:limit]

    if not cif_files:
        click.echo("Error: No cleaned CIF files found in input directory", err=True)
        sys.exit(1)

    click.echo(f"Found {len(cif_files)} cleaned CIF files")
    click.echo()

    # Initialize components
    click.echo("Loading ESM-2 model (this may take a moment)...")
    cleaner = GemmiStructureCleaner(sabdab_summary_path=sabdab_path)
    extractor = GeometricEpitopeExtractor(distance_threshold=distance_threshold)
    encoder = ESM2EpitopeEncoder(
        device=device,
        use_fp16=not no_fp16,
        cache_dir=cache_path
    )
    store = HDF5EmbeddingStore()

    # Process structures
    structures = {}
    epitopes = []
    encoder_outputs = []
    errors = []

    with click.progressbar(cif_files, label='Processing structures') as bar:
        for cif_path in bar:
            try:
                # Clean structure (skip filtering since already cleaned)
                cleaned, _ = cleaner.clean_structure(cif_path, skip_filtering=True)
                if cleaned is None:
                    continue

                pdb_id = cleaned.pdb_id
                structures[pdb_id] = cleaned

                # Extract epitope
                epitope = extractor.extract_epitope(cleaned)
                epitopes.append(epitope)

                # Encode
                output = encoder.encode_full(epitope, cleaned)
                encoder_outputs.append(output)

            except Exception as e:
                errors.append(f"{cif_path.name}: {e}")

    click.echo()

    if not encoder_outputs:
        click.echo("Error: No structures encoded successfully!", err=True)
        if errors:
            for err in errors[:5]:
                click.echo(f"  - {err}", err=True)
        sys.exit(1)

    # Save embeddings to HDF5
    h5_path = output_path / "embeddings.h5"
    store.save_encoder_outputs(encoder_outputs, h5_path)
    click.echo(f"Embeddings saved to: {h5_path}")

    # Save CSV logs
    save_epitope_residues_csv(epitopes, structures, output_path / "epitope_residues.csv")
    save_epitope_summary_csv(epitopes, structures, output_path / "epitope_summary.csv", encoder_outputs)
    save_embedding_stats_csv(encoder_outputs, output_path / "embedding_stats.csv")
    click.echo(f"CSV logs saved to: {output_path}")

    # Print report
    click.echo()
    click.echo(generate_epitope_report(epitopes, encoder_outputs))

    # Print errors if any
    if errors:
        click.echo()
        click.echo(f"Warnings ({len(errors)} issues):", err=True)
        for err in errors[:10]:
            click.echo(f"  - {err}", err=True)
        if len(errors) > 10:
            click.echo(f"  ... and {len(errors) - 10} more", err=True)

    # Final summary
    click.echo()
    click.echo("=" * 70)
    click.echo("Embedding Complete!")
    click.echo("=" * 70)
    click.echo(f"Structures encoded: {len(encoder_outputs)}/{len(cif_files)}")
    click.echo(f"Output directory: {output_path}")


@cli.command('filter-interactions')
@click.option(
    '--input', '-i', 'input_dir',
    required=True,
    help='Input directory containing DP_* folders (from build step)',
    type=click.Path(exists=True)
)
@click.option(
    '--output', '-o',
    required=True,
    help='Output directory for filtered DP_* folders and summary CSV',
    type=click.Path()
)
@click.option(
    '--distance-threshold',
    default=5.0,
    show_default=True,
    type=float,
    help='Maximum atom-atom distance (Å) to count as a contact'
)
@click.option(
    '--min-contacts',
    default=10,
    show_default=True,
    type=int,
    help='Minimum contact pairs required to keep a data point'
)
@click.option(
    '--dry-run',
    is_flag=True,
    default=False,
    help='Do not copy data; only produce summary'
)
def filter_interactions_command(
    input_dir: str,
    output: str,
    distance_threshold: float,
    min_contacts: int,
    dry_run: bool,
):
    """
    Filter DP_* folders to keep only antibody-antigen pairs with contacts.

    Uses atom-atom distance to confirm interaction and writes filter_summary.csv.
    """
    if os.path.abspath(input_dir) == os.path.abspath(output) and not dry_run:
        click.echo("Error: output directory must differ from input (or use --dry-run).", err=True)
        raise click.Abort()

    click.echo("=" * 60)
    click.echo("Antibody-Antigen Interaction Filter")
    click.echo("=" * 60)
    click.echo(f"Input: {input_dir}")
    click.echo(f"Output: {output}")
    click.echo(f"Distance threshold: {distance_threshold} Å")
    click.echo(f"Min contacts: {min_contacts}")
    click.echo(f"Dry run: {dry_run}")
    click.echo()

    results = filter_dataset(
        input_dir=input_dir,
        output_dir=output,
        distance_threshold=distance_threshold,
        min_contacts=min_contacts,
        dry_run=dry_run,
    )

    kept = sum(1 for r in results if r.passed)
    total = len(results)
    click.echo(f"Checked {total} data points; kept {kept}.")
    summary_path = os.path.join(output, "filter_summary.csv")
    if dry_run:
        click.echo("Dry run complete (no data copied).")
    click.echo(f"Summary: {summary_path}")


# Main entry point - use cli group
def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()
