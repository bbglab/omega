# bayes run cli
# --input: table with columns chr | pos | ref | mut | context | gene | impact | samples (reads) | samples (mutabilities)
# --impact: grouping of impacts json dict
# --genes: grouping of genes json dict
# --samples: grouping of samples json dict

# mle run cli
# --input: table with columns chr | pos | ref | mut | context | gene | impact | samples (reads) | samples (mutabilities)
# --impact: grouping of impacts json dict
# --genes: grouping of genes json dict
# --samples: grouping of samples json dict

# notes: 
# use typer for cli

# input mutations should be given in some canonical order, e.g.
# chr, pos, alt

# mutability col in the input
# is calculated using the no. syn mutations and depth per position

import os
import typer
import tqdm
from multiprocessing import Pool

import pandas as pd
import numpy as np

from assemble import Configurator, Assembler
from utils import dict_append
from omega import dNdS

import warnings
warnings.filterwarnings(module='tensorflow*', action='ignore')

app = typer.Typer()


def bayes_infer(args):
    
    gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n = args

    res = {}
    
    dnds_calculator = dNdS(l, n)

    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]

    try:
        chain = dnds_calculator.bayes_run()
        res['mean_dnds'] = [np.mean(chain)]
        res['perc_25_dnds'] = [np.percentile(chain, 25)]
        res['perc_75_dnds'] = [np.percentile(chain, 75)]
    except:
        res['mean_dnds'] = [None]
        res['perc_25_dnds'] = [None]
        res['perc_75_dnds'] = [None]
    
    return res


def mle_infer(args):

    gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n = args

    res = {}
    
    dnds_calculator = dNdS(l, n)

    res['gene'] = [gene_term]
    res['sample'] = [sample_term]
    res['impact'] = [impact_term]

    try:
        omega_hat, pvalue = dnds_calculator.mle_run()
        res['dnds'] = [omega_hat]
        res['pvalue'] = [pvalue]
    except:
        res['dnds'] = [None]
        res['pvalue'] = [None]
        
    return res


@app.command()
def bayes(conf_folder: str, data_fn: str, mut_counts_fn: str, output_fn: str, cores=4):

    conf = Configurator()

    # conf collects the groupings of samples, genes and impacts
    # the conf object will be passed to the Assembler so that it can create a data grid in accordance with the groupings

    conf.add_conf('samples', os.path.join(conf_folder, 'group_samples.json'))
    conf.add_conf('genes', os.path.join(conf_folder, 'group_genes.json'))
    conf.add_conf('impacts', os.path.join(conf_folder, 'group_impacts.json'))
    
    # load data: depth per site per sample
    # for testing data_fn = '../test/wide_input.tsv.gz'

    data = pd.read_csv(data_fn, sep='\t')

    # load data: mutation counts
    # for testing mut_counts_fn = '/workspace/datasets/transfer/ferran_to_ferriol/all_mutations_per_gene_impact_context.tsv'

    mut_counts = pd.read_csv(mut_counts_fn, sep='\t')

    # prepare input grid
    # input_grid is a generator that spits tuples
    # (gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n)
    # covering all the analysis cases specified by the groupings

    input_grid = list(Assembler(data, mut_counts, conf).input_generator())

    res = {}
    with Pool(4) as p:
        for d in tqdm.tqdm(p.imap(bayes_infer, input_grid), total=len(input_grid)):
            res = dict_append(res, d)

    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)


@app.command()
def mle(conf_folder: str, data_fn: str, mut_counts_fn: str, output_fn: str, cores=4):

    conf = Configurator()

    # conf collects the groupings of samples, genes and impacts
    # the conf object will be passed to the Assembler so that it can create a data grid in accordance with the groupings

    conf.add_conf('samples', os.path.join(conf_folder, 'group_samples.json'))
    conf.add_conf('genes', os.path.join(conf_folder, 'group_genes.json'))
    conf.add_conf('impacts', os.path.join(conf_folder, 'group_impacts.json'))
    
    # load data: depth per site per sample
    # for testing data_fn = '../test/wide_input.tsv.gz'

    data = pd.read_csv(data_fn, sep='\t')

    # load data: mutation counts
    # for testing mut_counts_fn = '/workspace/datasets/transfer/ferran_to_ferriol/all_mutations_per_gene_impact_context.tsv'

    mut_counts = pd.read_csv(mut_counts_fn, sep='\t')

    # prepare input grid
    # input_grid is a generator that spits tuples
    # (gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n)
    # covering all the analysis cases specified by the groupings

    input_grid = list(Assembler(data, mut_counts, conf).input_generator())

    res = {}

    for args in tqdm.tqdm(input_grid):
        
        d = mle_infer(args)
        res = dict_append(res, d)

    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)


if __name__ == "__main__":

    """
    python src/main.py bayes test/ \
    test/wide_input.tsv.gz \
    /workspace/datasets/transfer/ferran_to_ferriol/all_mutations_per_gene_impact_context.tsv \
    test/output/results_bayes.tsv \
    --cores 2
    """

    """
    python src/main.py mle test/ \
    test/wide_input.tsv.gz \
    /workspace/datasets/transfer/ferran_to_ferriol/all_mutations_per_gene_impact_context.tsv \
    test/output/results_mle.tsv \
    --cores 2
    """
    
    app()