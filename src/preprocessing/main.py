import os
import click
import sys
from sites_table_from_bed import generate_all_sites_4VEP
from postprocessing_annotation import vep2summarizedannotation
from compute_mutabilities import compute_mutabilities_wrapper

@click.command()
@click.option('--preprocessing-mode', type=click.Choice(['run_vep', 'postprocess_vep', 'compute_mutabilities']), help='Preprocessing mode')
@click.option('--depths-file', type=click.Path(exists=True), help='Path to depths file')
@click.option('--mutations-file', type=click.Path(exists=True), help='Path to mutations file')
@click.option('--bed-regions-file', type=click.Path(exists=True), help='Path to BED regions file')
@click.option('--vep-input-generated', type=click.Path(), help='Path to output file generated for VEP to annotate it')

@click.option('--vep-output-file', type=click.Path(exists=True), help='Path to VEP output file')
@click.option('--vep-postprocessed-file', type=click.Path(exists=True), help='Path to postprocessed VEP output file')

@click.option('--input-vep-postprocessed-file', type=click.Path(exists=True), help='Path to postprocessed VEP file')

@click.option('--table-observed-muts', type=click.Path(), help='Path to table of observed mutations file. We recommend: mutability_per_sample_gene_context.tsv')
@click.option('--mutabilities-table', type=click.Path(), help='Path to mutabilities table file. We recommend: mutations_per_sample_gene_impact_context.count.tsv')
def main(preprocessing_mode,
            bed_regions_file, vep_input_generated,
            vep_output_file, vep_postprocessed_file,
            input_vep_postprocessed_file,
            depths_file, mutations_file,
            table_observed_muts, mutabilities_table):
    
    run_vep = False
    postprocess_vep = False
    compute_mutabilities = False
    # Adjust based on preprocessing mode
    if preprocessing_mode == 'run_vep':
        run_vep = True
        required_inputs = [bed_regions_file, vep_input_generated]

    elif preprocessing_mode == 'postprocess_vep':
        postprocess_vep = True
        required_inputs = [vep_output_file, vep_postprocessed_file]

    elif preprocessing_mode == 'compute_mutabilities':
        compute_mutabilities = True
        required_inputs = [depths_file, mutations_file, input_vep_postprocessed_file]

    if not all(required_inputs):
        click.echo("One or more input files do not exist. Please check the paths.")
        sys.exit(1)

    if run_vep:
        generate_all_sites_4VEP(bed_regions_file, vep_input_generated)

    if not vep_output_file and not input_vep_postprocessed_file:
        click.echo("Use this file to run Ensembl VEP:")
        click.echo("./regions_sites.4VEP.tsv")
        click.echo("You can run it here: https://www.ensembl.org/Homo_sapiens/Tools/VEP, and then download the results.")
        vep_output_file = click.prompt("Introduce the path to the EnsemblVEP annotated file")

    if postprocess_vep:
        click.echo(f"The postprocessed VEP output will be stored at: {vep_postprocessed_file}")
        vep2summarizedannotation(vep_output_file, vep_postprocessed_file)

    if compute_mutabilities:
        compute_mutabilities_wrapper(input_vep_postprocessed_file, depths_file, mutations_file, table_observed_muts, mutabilities_table )

if __name__ == '__main__':
    main()
