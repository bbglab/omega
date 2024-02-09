import os

import sys
from omega.src.preprocessing.sites_table_from_bed import generate_all_sites_4VEP
from omega.src.preprocessing.postprocessing_annotation import vep2summarizedannotation
from omega.src.preprocessing.compute_mutabilities import compute_mutabilities_wrapper

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

    # TODO revise how to handle this case
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
