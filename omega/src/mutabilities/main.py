import os
import daiquiri
import warnings

import pandas as pd


from omega import __logger_name__, __version__
from omega.src.assemble import Grouping, Assembler

logger = daiquiri.getLogger(__logger_name__ + '.mutabilities')

warnings.filterwarnings(module='tensorflow*', action='ignore')


def mutabilities_per_site(d: dict) -> pd.DataFrame:
    """
    This function computes the mutabilities per site using the provided data files.

    It reads the mutability, depths, and VEP annotation files, and then creates a data grid based on the grouping of
    samples, genes, and impacts. The resulting data grid is returned as a DataFrame.
    """

    mutability = pd.read_csv(d['mutability_file'], sep='\t')
    depths = pd.read_csv(d['depths_file'], sep='\t')
    vep = pd.read_csv(d['vep_annotation_file'], sep='\t')

    # group collects the grouping of samples, genes and impacts
    # the group instance will be passed to the Assembler so that 
    # it can create a data grid in accordance with the grouping
    group = Grouping()
    group.add_group('samples', os.path.join(d['grouping_folder'], 'group_samples.json'))
    group.add_group('genes', os.path.join(d['grouping_folder'], 'group_genes.json'))

    ground_control = Assembler(depths, vep, mutability, group, mode='mutabilities')
    return ground_control.mutabilities_per_site


def main( mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, cores):

    d = {'mutability_file': mutability_file,
            'depths_file': depths_file,
            'vep_annotation_file': vep_annotation_file,
            'grouping_folder' : grouping_folder
            }

    df_mutabilities = mutabilities_per_site(d)
    df_mutabilities.to_csv(output_fn, header = True, index = False, sep = '\t')
    logger.info("Mutabilities per site stored")



