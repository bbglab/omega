#!/usr/bin/env python

""" 

# =============
# EXAMPLE USAGE
# =============

"""



import click
import daiquiri

from omega import __logger_name__, __version__

from omega.src.preprocessing.main import main as preprocessing_main
from omega.src.estimator.main import run_click as estimator_main

from omega.src.globals import DATE, setup_logging_decorator, startup_message

logger = daiquiri.getLogger(__logger_name__)


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
@click.option('--mutational-profile', type=click.Path(), default = None, help='Path to table of mutational profile.')
@click.option('--genome-assembly', type=click.Choice(['hg38', 'hg19', 'mm10']), default = 'hg38', help='Genome assembly')
@click.option('--single-sample', type=click.STRING, default = None, help='Name of the single sample. It also serves for activating the single sample mode.')
@click.option('--absent-synonymous', type=click.Choice(['ignore', 'infer_global_custom', 'infer_covariates']), default = 'ignore', help='Omega mode for genes without synonymous mutations.')
@click.option('--relative-synonymous-muts-file', type=click.Path(), default = None, help='Path to table of synonymous mutations per gene.')
@setup_logging_decorator
def preprocessing(preprocessing_mode,
                    bed_regions_file, vep_input_generated,
                    vep_output_file, vep_postprocessed_file,
                    input_vep_postprocessed_file,
                    depths_file, mutations_file,
                    table_observed_muts, mutabilities_table,
                    mutational_profile,
                    genome_assembly,
                    single_sample,
                    absent_synonymous,
                    relative_synonymous_muts_file
            ):
    """"Build datasets necessary to run Omega."""
    startup_message(__version__, "Initializing preprocessing...")

    logger.info("Received arguments:")
    logger.info(f"preprocessing_mode: {preprocessing_mode}")
    logger.info(f"bed_regions_file: {bed_regions_file}")
    logger.info(f"vep_input_generated: {vep_input_generated}")
    logger.info(f"vep_output_file: {vep_output_file}")
    logger.info(f"vep_postprocessed_file: {vep_postprocessed_file}")
    logger.info(f"input_vep_postprocessed_file: {input_vep_postprocessed_file}")
    logger.info(f"depths_file: {depths_file}")
    logger.info(f"mutations_file: {mutations_file}")
    logger.info(f"table_observed_muts: {table_observed_muts}")
    logger.info(f"mutabilities_table: {mutabilities_table}")
    logger.info(f"mutational_profile: {mutational_profile}")
    logger.info(f"genome_assembly: {genome_assembly}")
    logger.info(f"single_sample: {single_sample}")
    logger.info(f"absent_synonymous: {absent_synonymous}")
    logger.info(f"relative_synonymous_muts_file: {relative_synonymous_muts_file}")

    preprocessing_main(preprocessing_mode,
                        bed_regions_file, vep_input_generated,
                        vep_output_file, vep_postprocessed_file,
                        input_vep_postprocessed_file,
                        depths_file, mutations_file,
                        table_observed_muts, mutabilities_table,
                        mutational_profile,
                        genome_assembly,
                        single_sample,
                        absent_synonymous,
                        relative_synonymous_muts_file)


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
@setup_logging_decorator
def estimator(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, option, cores):
    startup_message(__version__, "Running estimator...")

    logger.info("Received arguments:")
    logger.info(f"observed_mutations_file: {observed_mutations_file}")
    logger.info(f"mutability_file: {mutability_file}")
    logger.info(f"depths_file: {depths_file}")
    logger.info(f"vep_annotation_file: {vep_annotation_file}")
    logger.info(f"grouping_folder: {grouping_folder}")
    logger.info(f"output_fn: {output_fn}")
    logger.info(f"option: {option}")
    logger.info(f"cores: {cores}")
    estimator_main(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, option, cores)


if __name__ == "__main__":
    omega() 