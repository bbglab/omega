import os
import json
import typer
import tqdm
from enum import Enum
from multiprocessing import Pool

import pandas as pd
import numpy as np

from assemble import Grouping, Assembler
from omega import bayes_infer, mle_infer

import warnings
warnings.filterwarnings(module='tensorflow*', action='ignore')


app = typer.Typer()


def prepare_data(input_fn):
    
    with open(input_fn, 'rt') as f:
        d = json.load(f)

    mut_counts = pd.read_csv(d['observed_mutations_file'], sep='\t')
    mutability = pd.read_csv(d['mutability_file'], sep='\t')
    depths = pd.read_csv(d['depths_file'], sep='\t')
    regions = pd.read_csv(d['bed_regions_file'], sep='\t')
    vep = pd.read_csv(d['vep_annotation_file'], sep='\t')

    # group collects the grouping of samples, genes and impacts
    # the group instance will be passed to the Assembler so that 
    # it can create a data grid in accordance with the grouping
    group = Grouping()
    group.add_group('samples', os.path.join(d['grouping_folder'], 'group_samples.json'))
    group.add_group('genes', os.path.join(d['grouping_folder'], 'group_genes.json'))
    group.add_group('impacts', os.path.join(d['grouping_folder'], 'group_impacts.json'))

    ground_control = Assembler(depths, regions, vep, mut_counts, mutability, group)
    return list(ground_control.input_generator())
    

def bayes(input_json, output_fn, cores=4):

    input_grid = prepare_data(input_json)
    res = {}
    with Pool(4) as p:
        for d in tqdm.tqdm(p.imap(bayes_infer, input_grid), total=len(input_grid)):
            res = {k: res.get(k, []) + d.get(k, []) for k in d}
    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)


def mle(input_json: str, output_fn: str, cores=4):

    input_grid = prepare_data(input_json)
    res = {}
    for args in tqdm.tqdm(input_grid):
        d = mle_infer(args)
        res = {k: res.get(k, []) + d.get(k, []) for k in d}
    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)


class ModelType(str, Enum):
    bayes = "bayes"
    mle = "mle"


@app.command()
def run(input_json: str, output_fn: str, option: ModelType=ModelType.bayes, cores=4):

    if option == 'bayes':
        bayes(input_json, output_fn, cores=cores)
    if option == 'mle':
        mle(input_json, output_fn, cores=cores)


if __name__ == "__main__":

    """
    python src/main.py --option bayes --cores 2 test/input_estimation.json test/output_estimation_bayes.tsv
    python src/main.py --option mle --cores 2 test/input_estimation.json test/output_estimation_mle.tsv
    """
    
    app()