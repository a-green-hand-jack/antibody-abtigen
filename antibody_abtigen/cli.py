"""
Command-line interface for antibody-abtigen.
"""

import os
import sys
import click

from .pipeline import CrossSpeciesDatasetPipeline


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
def build_command(
    output: str,
    data_dir: str,
    limit: int,
    resolution: float,
    identity: float,
    dry_run: bool,
    no_pymol: bool,
    force_download: bool
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
        use_pymol=not no_pymol
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


# Main entry point - use cli group
def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()
