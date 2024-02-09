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
# @click.option("-c", "--cores", type=click.IntRange(min=1, max=len(os.sched_getaffinity(0)), clamp=False), default=len(os.sched_getaffinity(0)),
#                 help="Number of cores to use in the computation")
# @click.option("-v", "--verbose", help="Verbose", is_flag=True)
# @setup_logging_decorator
def preprocessing(preprocessing_mode,
                    bed_regions_file, vep_input_generated,
                    vep_output_file, vep_postprocessed_file,
                    input_vep_postprocessed_file,
                    depths_file, mutations_file,
                    table_observed_muts, mutabilities_table
                #    cores, af_version,
                #    keep_pdb_files, 
                #    yes,
                #    verbose
            ):
    """"Build datasets necessary to run Omega."""
    
    # startup_message(__version__, "Initializing building datasets...")
    
    # logger.info(f"Current working directory: {os.getcwd()}")
    # logger.info(f"Build folder path: {output_dir}")
    # logger.info(f"Organism: {organism}")
    # logger.info(f"Distance threshold: {distance_threshold}Å")
    # logger.info(f"Custom IDs mapping: {uniprot_to_hugo}")
    # logger.info(f"CPU cores: {cores}")
    # logger.info(f"AlphaFold version: {af_version}")
    # logger.info(f"Keep PDB files: {keep_pdb_files}")
    # logger.info(f"Verbose: {verbose}")
    # logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    # logger.info("")

    preprocessing_main(preprocessing_mode,
                        bed_regions_file, vep_input_generated,
                        vep_output_file, vep_postprocessed_file,
                        input_vep_postprocessed_file,
                        depths_file, mutations_file,
                        table_observed_muts, mutabilities_table)

    # preprocessing_main(output_dir, 
    #                     organism, 
    #                     distance_threshold,
    #                     uniprot_to_hugo, 
    #                     cores, 
    #                     af_version, 
    #                     keep_pdb_files)


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
def estimator(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn,
                    option,
                    cores):


# @click.option("-i", "--input_maf_path", type=click.Path(exists=True), required=True, help="Path of the MAF file used as input")
# @click.option("-p", "--mut_profile_path", type=click.Path(exists=True), help="Path of the mutation profile (192 trinucleotide contexts) used as optional input")
# @click.option("-m", "--mutability_config_path", type=click.Path(exists=True), help="Path of the config file with information on mutability")
# @click.option("-o", "--output_dir", help="Path to output directory", type=str, default='results')
# @click.option("-d", "--data_dir", help="Path to datasets", type=click.Path(exists=True), default = os.path.join('datasets'))
# @click.option("-n", "--n_iterations", help="Number of densities to be simulated", type=int, default=10000)
# @click.option("-a", "--alpha", help="Significant threshold for the p-value of res and gene", type=float, default=0.01)
# @click.option("-P", "--cmap_prob_thr", type=float, default=0.5,
#               help="Threshold to define AAs contacts based on distance on predicted structure and PAE")
# @click.option("-f", "--no_fragments", help="Disable processing of fragmented (AF-F) proteins", is_flag=True)
# @click.option("-x", "--only_processed", help="Include only processed genes in the output", is_flag=True)
# @click.option("-y", "--thr_not_in_structure", type=float, default=0.1,
#               help="Threshold to filter out genes based on the ratio of mutations outside of the structure")
# @click.option("-c", "--cores", type=click.IntRange(min=1, max=len(os.sched_getaffinity(0)), clamp=False), default=len(os.sched_getaffinity(0)),
#               help="Set the number of cores to use in the computation")
# @click.option("-s", "--seed", help="Set seed to ensure reproducible results", type=int)
# @click.option("-v", "--verbose", help="Verbose", is_flag=True)
# @click.option("-t", "--cancer_type", help="Cancer type", type=str)
# @click.option("-C", "--cohort", help="Name of the cohort", type=str)
# @setup_logging_decorator

    estimator_main(observed_mutations_file, mutability_file,
                    depths_file,
                    vep_annotation_file,
                    grouping_folder, 
                    output_fn,
                    option,
                    cores)


# def run(input_maf_path, 
#         mut_profile_path,
#         mutability_config_path,
#         output_dir,
#         data_dir,
#         n_iterations,
#         alpha,
#         cmap_prob_thr,
#         no_fragments,
#         only_processed,
#         thr_not_in_structure,
#         cores,
#         seed,
#         verbose,
#         cancer_type,
#         cohort):
    """Run Oncodrive3D."""

    # ## Initialize
    # plddt_path = os.path.join(data_dir, "confidence.csv")
    # cmap_path = os.path.join(data_dir, "prob_cmaps")  
    # seq_df_path = os.path.join(data_dir, "seq_for_mut_prob.csv")                              
    # pae_path = os.path.join(data_dir, "pae")
    # cancer_type = cancer_type if cancer_type else np.nan
    # cohort = cohort if cohort else f"cohort_{DATE}"
    # path_prob = mut_profile_path if mut_profile_path else "Not provided, mutabilities will be used" if mutability_config_path else "Not provided, uniform distribution will be used"
    # path_mutability_config = mutability_config_path if mutability_config_path else "Not provided, mutabilities will not be used"

    # # Log
    # startup_message(__version__, "Initializing analysis...")

    # logger.info(f"Input MAF: {input_maf_path}")
    # logger.info(f"Input mut profile: {path_prob}")
    # logger.info(f"Input mutability config: {path_mutability_config}")
    # logger.info(f"Build directory: {data_dir}")
    # logger.info(f"Output directory: {output_dir}")
    # logger.debug(f"Path to CMAPs: {cmap_path}")
    # logger.debug(f"Path to DNA sequences: {seq_df_path}")
    # logger.debug(f"Path to PAE: {pae_path}")
    # logger.debug(f"Path to pLDDT scores: {plddt_path}")
    # logger.info(f"CPU cores: {cores}")
    # logger.info(f"Iterations: {n_iterations}")
    # logger.info(f"Significant level: {alpha}")
    # logger.info(f"Probability threshold for CMAPs: {cmap_prob_thr}")
    # logger.info(f"Disable fragments: {bool(no_fragments)}")
    # logger.info(f"Output only processed genes: {bool(only_processed)}")
    # logger.info(f"Ratio threshold mutations out of structure: {thr_not_in_structure}")
    # logger.info(f"Cohort: {cohort}")
    # logger.info(f"Cancer type: {cancer_type}")
    # logger.info(f"Verbose: {bool(verbose)}")
    # logger.info(f"Seed: {seed}")
    # logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    # logger.info("")


if __name__ == "__main__":
    omega() 