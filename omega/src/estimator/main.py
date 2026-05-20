import multiprocessing
import os
import warnings
from functools import partial

import daiquiri
import pandas as pd
import tqdm

from omega import __logger_name__
from omega.src.assemble import Assembler, Grouping
from omega.src.estimator.omega import mle_infer, bayes_infer
# Minimum non-zero p-value when clamping to float32 limits.
MIN_NORMALIZED_FLOAT32_PVALUE = 1.17e-38

LOG = daiquiri.getLogger(__logger_name__ + '.estimator')

warnings.filterwarnings(module='tensorflow*', action='ignore')

def _init_worker(gpu_id: int | None = None) -> None:
    """
    Worker initializer. Runs once in each child process.
    Set environment for devices and optionally configure TF threading.
    Import TensorFlow here if needed (after device selection).
    """
    if gpu_id is not None:
        # Make only this GPU visible to this process
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    # Import TF inside worker to avoid inheriting parent's TF runtime.
    try:
        import tensorflow as tf  # type: ignore

        tf.config.threading.set_intra_op_parallelism_threads(1)
        tf.config.threading.set_inter_op_parallelism_threads(1)
    except Exception:
        # If TF is not required by some jobs, ignore import errors here.
        pass


def prepare_data(d):
    mut_counts = pd.read_csv(d['observed_mutations_file'], sep='\t')
    mutability = pd.read_csv(d['mutability_file'], sep='\t')
    depths = pd.read_csv(d['depths_file'], sep='\t')
    vep = pd.read_csv(d['vep_annotation_file'], sep='\t')

    # group collects the grouping of samples, genes and impacts
    # the group instance will be passed to the Assembler so that
    # it can create a data grid in accordance with the grouping
    group = Grouping()
    group.add_group('samples', os.path.join(d['grouping_folder'], 'group_samples.json'))
    group.add_group('genes', os.path.join(d['grouping_folder'], 'group_genes.json'))
    group.add_group('impacts', os.path.join(d['grouping_folder'], 'group_impacts.json'))

    ground_control = Assembler(depths, vep, mutability, group, mut_counts, mode='estimator')
    return list(ground_control.input_generator())

def bayes(input_json, output_fn, cores=4):
    LOG.info("Running in bayes mode")

    input_grid = prepare_data(input_json)
    LOG.info("Data prepared")
    res = {}
    ctx = multiprocessing.get_context('spawn')
    with ctx.Pool(cores) as p:
        for d in tqdm.tqdm(p.imap(bayes_infer, input_grid), total=len(input_grid)):
            res = {k: res.get(k, []) + d.get(k, []) for k in d}
    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)

def mle(input_json: dict[str, str], output_fn: str, dispersion: float, cores=4):
    LOG.info('Running in mle mode')

    input_grid = prepare_data(input_json)
    LOG.info('Data prepared')

    res = {}
    if cores == 1:
        for args in tqdm.tqdm(input_grid):
            d = mle_infer(args, dispersion)
            res = {k: res.get(k, []) + d.get(k, []) for k in d}

    else:
        worker_fn = partial(mle_infer, dispersion=dispersion)
        ctx = multiprocessing.get_context('spawn')
        with ctx.Pool(processes=cores, initializer=_init_worker) as p:
            for d in tqdm.tqdm(p.imap(worker_fn, input_grid), total=len(input_grid)):
                res = {k: res.get(k, []) + d.get(k, []) for k in d}

    LOG.info('Writing results to  %s ', output_fn)
    df = pd.DataFrame(res)

    LOG.info('Replacing 0 p-values with %g in  %s ', MIN_NORMALIZED_FLOAT32_PVALUE, output_fn)
    df['pvalue'] = df['pvalue'].apply(lambda x: x if x > 0 else MIN_NORMALIZED_FLOAT32_PVALUE)

    df.to_csv(output_fn, sep='\t', index=False)


def main(
    observed_mutations_file,
    mutability_file,
    depths_file,
    vep_annotation_file,
    grouping_folder,
    output_fn,
    dispersion_value,
    option,
    cores,
):
    d: dict[str, str] = {
        'observed_mutations_file': observed_mutations_file,
        'mutability_file': mutability_file,
        'depths_file': depths_file,
        'vep_annotation_file': vep_annotation_file,
        'grouping_folder': grouping_folder,
    }

    if option == 'bayes':
        bayes(d, output_fn, cores= int(cores))
    if option == 'mle':
        mle(d, output_fn, dispersion_value, cores=int(cores))
