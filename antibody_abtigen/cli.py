"""
Command-line interface for antibody-abtigen.
"""

import os
import sys
import click

from .pipeline import CrossSpeciesDatasetPipeline
from .filtering import filter_dataset


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
