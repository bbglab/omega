#!/usr/bin/env python

""" 
The module includes the the main function to run an HotMAPs-inspired 
method that uses the the canonical predicted structure stored in 
AlphaFold db to perform 3D-clustering of mutations using simulations 
and rank based comparison.

# =============
# EXAMPLE USAGE
# =============

# Build datasets

cd path/to/oncodrive3D
oncodrive3D build-datasets

# Run 
            
oncodrive3D run \
    -i test/maf/TCGA_WXS_ACC.in.maf  \
        -p test/mut_profile/TCGA_WXS_ACC.mutrate.json \
            -o test/results
                  
oncodrive3D run -i /workspace/projects/clustering_3d/clustering_3d/datasets_normal/kidneydata/all_mutations.all_samples.tsv -d /workspace/projects/clustering_3d/clustering_3d/datasets_normal -m /workspace/projects/clustering_3d/clustering_3d/test/normal_tests/kidneydata/mutability_kidney.json -o /workspace/projects/clustering_3d/dev_testing/result/o3d/test_normal -C kidney_mutability_unif -v
oncodrive3D run -i /workspace/projects/clustering_3d/clustering_3d/datasets_normal/kidneydata/all_mutations.all_samples.tsv -d /workspace/projects/clustering_3d/clustering_3d/datasets_normal -p /workspace/projects/clustering_3d/clustering_3d/test/normal_tests/kidneydata/mut_profiles/all_samples.192.json           
"""


# =============================================================================
# TODO: allow procesing without tumor sample info
# TODO: handle inf of the score (e.g., decimal python package)
# TODO: fix bug in requirement.txt (bgreference must be installed after setup.py)
# TODO: add script to generate conf and mutability file
# TODO: test run on normal tissue
# TODO: add filter (and logs) for mutated genes without exons coordinate when mutability is provided

# TODO: change progressbar to tqdm in run scripts
# TODO: change output names?
# TODO: change repo name
# TODO: fix doc
# TODO: update doc for normal tissue application
# TODO: suppress verbosity of multi-threading download of structures
# =============================================================================


import json
import os

import click
# import daiquiri
import numpy as np
import pandas as pd

from omega.src import __logger_name__, __version__

from omega.src.preprocessing.main import main as preprocessing_main
from omega.src.estimator.main import run_click as estimator_main

# from scripts.globals import DATE, setup_logging_decorator, startup_message

# logger = daiquiri.getLogger(__logger_name__)


@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(__version__)
def omega():
    """Omega: software for the computation of dNdS."""
    pass


@omega.command(context_settings=dict(help_option_names=['-h', '--help']),
                help="Build datasets - Required once after installation.") 
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
def preprocessing(preprocessing_mode,
                    bed_regions_file, vep_input_generated,
                    vep_output_file, vep_postprocessed_file,
                    input_vep_postprocessed_file,
                    depths_file, mutations_file,
                    table_observed_muts, mutabilities_table
            ):
    """"Build datasets necessary to run Omega."""
    # startup_message(__version__, "Initializing preprocessing...")

    # logger.info("")
    preprocessing_main(preprocessing_mode,
                        bed_regions_file, vep_input_generated,
                        vep_output_file, vep_postprocessed_file,
                        input_vep_postprocessed_file,
                        depths_file, mutations_file,
                        table_observed_muts, mutabilities_table)


@omega.command(context_settings=dict(help_option_names=['-h', '--help']),
                        help="Run dNdS analysis.")
@click.option('--observed-mutations-file', type=click.Path(exists=True), help='Path to observed mutations file')
@click.option('--mutability-file', type=click.Path(exists=True), help='Path to mutability file')
@click.option('--depths-file', type=click.Path(exists=True), help='Path to depths file')
@click.option('--vep-annotation-file', type=click.Path(exists=True), help='Path to VEP annotation file')
@click.option('--grouping-folder', type=click.Path(exists=True), help='Path to grouping folder')
@click.option('--output-fn', type=str, help='Output filename')
@click.option('--option', type=click.Choice(['bayes', 'mle']), default='bayes', help='Option type (default: bayes)')
@click.option('--cores', type=int, default=4, help='Number of cores (default: 4)')
def estimator(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, option, cores):
    # startup_message(__version__, "Running estimator...")

    # logger.info("")

    print("Received arguments:")
    print(f"observed_mutations_file: {observed_mutations_file}")
    print(f"mutability_file: {mutability_file}")
    print(f"depths_file: {depths_file}")
    print(f"vep_annotation_file: {vep_annotation_file}")
    print(f"grouping_folder: {grouping_folder}")
    print(f"output_fn: {output_fn}")
    print(f"option: {option}")
    print(f"cores: {cores}")
    estimator_main(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, option, cores)


if __name__ == "__main__":
    omega() 