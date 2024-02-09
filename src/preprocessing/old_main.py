import os
import json
import sys
import click

from sites_table_from_bed import generate_all_sites_4VEP
from postprocessing_annotation import vep2summarizedannotation
from compute_mutabilities import compute_mutabilities_wrapper


@click.command()
@click.option('--depths-file', type=click.Path(exists=True), help='Path to depths file')
@click.option('--mutations-file', type=click.Path(exists=True), help='Path to mutations file')
@click.option('--vep-output-file', type=click.Path(exists=True), help='Path to VEP output file')
@click.option('--vep-postprocessed-file', type=click.Path(exists=True), help='Path to postprocessed VEP output file')
@click.option('--bed-regions-file', type=click.Path(exists=True), help='Path to BED regions file')
@click.option('--table-observed-muts', type=click.Path(), help='Path to table of observed mutations file')
@click.option('--mutabilities-table', type=click.Path(), help='Path to mutabilities table file')
def main(depths_file, mutations_file, vep_output_file, vep_postprocessed_file, bed_regions_file, table_observed_muts, mutabilities_table):
    # Adjust based on input files
    run_vep = not vep_output_file and not vep_postprocessed_file and bed_regions_file
    postprocess_vep = not vep_postprocessed_file
    compute_mutabilities = True

    if not (depths_file and mutations_file):
        click.echo("depths-file and mutations-file are required. Revise your input.")
        sys.exit(1)

    if not all([os.path.exists(file) for file in [depths_file, mutations_file, vep_output_file, vep_postprocessed_file, bed_regions_file]]):
        click.echo("One or more input files do not exist. Please check the paths.")
        sys.exit(1)

    if run_vep:
        generate_all_sites_4VEP(bed_regions_file, './regions_sites.4VEP.tsv')

    if not vep_output_file:
        click.echo("Use this file to run Ensembl VEP:")
        click.echo("./regions_sites.4VEP.tsv")
        click.echo("You can run it here: https://www.ensembl.org/Homo_sapiens/Tools/VEP, and then download the results.")
        vep_output_file = click.prompt("Introduce the path to the EnsemblVEP annotated file")

    if postprocess_vep:
        click.echo(f"The postprocessed VEP output will be stored at: {vep_postprocessed_file}")
        vep2summarizedannotation(vep_output_file, vep_postprocessed_file)

    if compute_mutabilities:
        compute_mutabilities_wrapper(vep_postprocessed_file, depths_file, mutations_file, table_observed_muts or "./mutations_per_sample_gene_impact_context.count.tsv", mutabilities_table or "./mutability_per_sample_gene_context.tsv")

if __name__ == '__main__':
    main()








def main():

    # Read JSON data from the file
    with open(sys.argv[1], 'r') as json_file:
        input_data = json.load(json_file)


    ####
    # define INPUT files and steps to run
    ####

    for info in ['mutations_file', 'depths_file']:
        if info not in input_data.keys():
            print(f"{info} is required. Revise your input file.")
            sys.exit(1)

    depths_file = input_data['depths_file']
    mutations_file = input_data['mutations_file']

    # adjust based on input files
    run_vep = True
    postprocess_vep = True
    compute_mutabilities = True


    if 'vep_postprocessed_file' in input_data.keys() :
        if os.path.exists(input_data['vep_postprocessed_file']):
            run_vep = False
            postprocess_vep = False

            print(f"Using the precomputed and postprocessed VEP output available at:\n{input_data['vep_postprocessed_file']}")

            vep_postprocessed_file = input_data['vep_postprocessed_file']
        else:
            print(f"The provided precomputed and postprocessed VEP output does NOT EXIST at:\n{input_data['vep_postprocessed_file']}")
            sys.exit(1)

    elif 'vep_output_file' in input_data.keys():
        if os.path.exists(input_data['vep_output_file']):
            run_vep = False
            print(f"Using the raw VEP output available at:\n{input_data['vep_output_file']}")

            vep_output_file = input_data['vep_output_file']
            vep_postprocessed_file = "./regions_sites.annotation_summary.tsv"

        else:
            print(f"The provided raw VEP output does NOT EXIST at:\n{input_data['vep_output_file']}")
            sys.exit(1)

    elif 'bed_regions_file' in input_data.keys():
        if os.path.exists(input_data['bed_regions_file']):
            print(f"Using the BED file available at:\n{input_data['bed_regions_file']}")

            bed_regions_file = input_data['bed_regions_file']
            vep_input_file = "./regions_sites.4VEP.tsv"
            vep_postprocessed_file = "./regions_sites.annotation_summary.tsv"

        else:
            print(f"The provided BED file does NOT EXIST at:\n{input_data['bed_regions_file']}")
            sys.exit(1)



    ####
    # define OUTPUT files
    ####
    # TODO: decide what to do if the tables are already available
    both_avail = 0
    if 'table_observed_muts' in input_data.keys():
        if os.path.exists(input_data['table_observed_muts']):
            print(f"The table with observed mutations has already been computed in the past.")
            both_avail += 1
        table_observed_muts = input_data['table_observed_muts']
    else:
        table_observed_muts = "./mutations_per_sample_gene_impact_context.count.tsv"

    if 'mutabilities_table' in input_data.keys():
        if os.path.exists(input_data['mutabilities_table']):
            print(f"The table with the mutabilities per sample, gene and context has already been computed in the past.")
            both_avail += 1
        mutabilities_table = input_data['mutabilities_table']
    else:
        mutabilities_table = "./mutability_per_sample_gene_context.tsv"

    if both_avail == 2:
        print("Both tables are already available, preprocessing not run.")
        sys.exit(0) # exit without error




    # TODO we could check for the presence of the chr prefix in the files we read just to make sure this is not causing problems
    if run_vep:
        generate_all_sites_4VEP(bed_regions_file, vep_input_file)


    ## run Ensembl VEP from their web server
    # TODO add description on how to run and download the EnsemblVEP output
    ### this would imply using "vep_input_file" to generate "vep_output_file"
    ### the user should introduce the path to vep_output_file

        print(f"Use this file to run Ensembl VEP:\n{vep_input_file}")
        print("You can run it here: https://www.ensembl.org/Homo_sapiens/Tools/VEP, and then download the results." )
        vep_output_file = input("Introduce the path to the EnsemblVEP annotated file:\n")
        # TODO check if the file exists
        vep_output_file = vep_output_file.strip()
        if not os.path.exists(vep_output_file):
            print(f"The provided raw VEP output does NOT EXIST at:\n{input_data['vep_output_file']}")
            sys.exit(1)

    if postprocess_vep:
        print(f"The postprocessed VEP output will be stored at:\n{vep_postprocessed_file}")
        vep2summarizedannotation(vep_output_file, vep_postprocessed_file)
    # else check if vep_postprocessed_file exists

    if compute_mutabilities:
        compute_mutabilities_wrapper(vep_postprocessed_file, depths_file, mutations_file, table_observed_muts, mutabilities_table)

if __name__ == '__main__':

    main()