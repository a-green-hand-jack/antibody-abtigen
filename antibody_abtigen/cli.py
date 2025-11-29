"""
Command-line interface for antibody-abtigen.
"""

import os
import sys
import click

from .pipeline import CrossSpeciesDatasetPipeline


@click.command()
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
@click.version_option()
def main(
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

    This tool downloads antibody-antigen complex data from SAbDab,
    finds mouse orthologs for human antigens, and generates aligned
    structure triplets for machine learning and computational biology.

    Output structure:

    \b
        output/
        ├── dataset_summary.csv
        ├── processing_log.json
        └── DP_XXXX_Y/
            ├── *_antibody.pdb/cif
            ├── *_human_ag.pdb/cif
            ├── *_mouse_ag.pdb/cif (aligned)
            └── metadata.json
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


if __name__ == '__main__':
    main()
