#!/usr/bin/env python

"""

# =============
# EXAMPLE USAGE
# =============

"""
import click
import daiquiri

from omega import __logger_name__, __version__
from omega.src.helpers import display_title_and_params
from omega.src.preprocessing.main   import main as preprocessing_main
from omega.src.estimator.main       import main as estimator_main
from omega.src.mutabilities.main    import main as mutabilities_main

from omega.src.globals              import setup_logging_decorator, startup_message

logger = daiquiri.getLogger(__logger_name__)


@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(__version__)
def omega():
    """Omega: software for the computation of dNdS."""
    pass


@omega.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True),
                help="Build input tables - Required once per cohort.", )
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
@click.option('--synonymous-muts-table', type=click.Path(), help='Path to table of observed synonymous mutations per gene file. We recommend: syn_muts_per_gene.tsv')
@click.option('--mutational-profile-file', type=click.Path(exists=True), help='Path to table of mutational profile.')
@click.option('--mutational-profile-global-file', type=click.Path(), default = None, help='Path to table of mutational profile for the "sample" in which global mutation rates were computed.')
@click.option('--genome-assembly', type=click.Choice(['hg38', 'hg19', 'mm10']), default = 'hg38', help='Genome assembly')
@click.option('--single-sample', type=click.STRING, default = None, help='Name of the single sample. It also serves for activating the single sample mode.')
@click.option('--absent-synonymous', type=click.Choice(['ignore', 'infer_global_custom', 'infer_covariates']), default = 'ignore', help='Omega mode for genes without synonymous mutations.')
@click.option('--synonymous-mutrates-file', type=click.Path(), default = None, help='Path to table of synonymous mutation rates per gene.')
@click.option('--verbose', is_flag=True, help='Enable verbose logging')
@setup_logging_decorator
def preprocessing(preprocessing_mode,
                    bed_regions_file, vep_input_generated,
                    vep_output_file, vep_postprocessed_file,
                    input_vep_postprocessed_file,
                    depths_file, mutations_file,
                    table_observed_muts, mutabilities_table,
                    synonymous_muts_table,
                    mutational_profile_file,
                    mutational_profile_global_file,
                    genome_assembly,
                    single_sample,
                    absent_synonymous,
                    synonymous_mutrates_file,
                    verbose
            ):
    """"Build tables necessary to run Omega."""
    startup_message(__version__, "mode: PREPROCESSING")

    display_title_and_params(title="Initializing preprocessing...")

    preprocessing_main(preprocessing_mode,
                        bed_regions_file, vep_input_generated,
                        vep_output_file, vep_postprocessed_file,
                        input_vep_postprocessed_file,
                        depths_file, mutations_file,
                        table_observed_muts, mutabilities_table,
                        synonymous_muts_table,
                        mutational_profile_file,
                        mutational_profile_global_file,
                        genome_assembly,
                        single_sample,
                        absent_synonymous,
                        synonymous_mutrates_file)


@omega.command(context_settings=dict(help_option_names=['-h', '--help']),
                        help="Run dNdS analysis.")
@click.option('--observed-mutations-file', type=click.Path(exists=True), help='Path to observed mutations file')
@click.option('--mutability-file', type=click.Path(exists=True), help='Path to mutability file')
@click.option('--depths-file', type=click.Path(exists=True), help='Path to depths file')
@click.option('--vep-annotation-file', type=click.Path(exists=True), help='Path to VEP annotation file')
@click.option('--grouping-folder', type=click.Path(exists=True), help='Path to grouping folder')
@click.option('--output-fn', type=str, help='Output filename')
@click.option('--dispersion', type=float, default=0.1, help='Choose dispersion value(default: 0.1)')
@click.option('--option', type=click.Choice(['bayes', 'mle']), default='bayes', help='Option type (default: bayes)')
@click.option('--cores', type=int, default=4, help='Number of cores (default: 4)')
@click.option('--verbose', is_flag=True, help='Enable verbose logging')
@setup_logging_decorator
def estimator(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, dispersion, option, cores, verbose):
    startup_message(__version__, "mode: ESTIMATOR")

    display_title_and_params(title="Running estimator...")

    estimator_main(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, dispersion, option, cores)


@omega.command(context_settings=dict(help_option_names=['-h', '--help']),
                        help="Run dNdS analysis.")
@click.option('--mutability-file', type=click.Path(exists=True), help='Path to mutability file')
@click.option('--depths-file', type=click.Path(exists=True), help='Path to depths file')
@click.option('--vep-annotation-file', type=click.Path(exists=True), help='Path to VEP annotation file')
@click.option('--grouping-folder', type=click.Path(exists=True), help='Path to grouping folder')
@click.option('--output-fn', type=str, help='Output filename')
@click.option('--cores', type=int, default=4, help='Number of cores (default: 4)')
@click.option('--verbose', is_flag=True, help='Enable verbose logging')
@setup_logging_decorator
def mutabilities( mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, cores, verbose):
    """Compute mutabilities per site."""
    startup_message(__version__, "mode: MUTABILITIES")

    display_title_and_params(title="Computing mutabilities per site...")

    mutabilities_main(mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, cores)


if __name__ == "__main__":
    omega()