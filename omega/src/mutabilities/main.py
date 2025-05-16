import os
import daiquiri
import warnings

import pandas as pd


from omega import __logger_name__, __version__
from omega.src.mutabilities.assemble import Grouping, Assembler

logger = daiquiri.getLogger(__logger_name__ + '.mutabilities')

warnings.filterwarnings(module='tensorflow*', action='ignore')


def mutabilities_per_site(d):

    mutability = pd.read_csv(d['mutability_file'], sep='\t')
    depths = pd.read_csv(d['depths_file'], sep='\t')
    vep = pd.read_csv(d['vep_annotation_file'], sep='\t')

    # group collects the grouping of samples, genes and impacts
    # the group instance will be passed to the Assembler so that 
    # it can create a data grid in accordance with the grouping
    group = Grouping()
    group.add_group('samples', os.path.join(d['grouping_folder'], 'group_samples.json'))
    group.add_group('genes', os.path.join(d['grouping_folder'], 'group_genes.json'))

    ground_control = Assembler(depths, vep, mutability, group)
    return ground_control.mutabilities_per_site

def get_mutabilities(input_json: str, output_fn: str,  cores=4):
    df_mutabilities = mutabilities_per_site(input_json)
    df_mutabilities.to_csv(output_fn, header = True, index = False, sep = '\t')
    logger.info("Mutabilities per site stored")

def main( mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn, cores):
    
    d = {'mutability_file': mutability_file,
            'depths_file': depths_file,
            'vep_annotation_file': vep_annotation_file,
            'grouping_folder' : grouping_folder
            }

    get_mutabilities(d, output_fn, cores= int(cores))



