import os
import tqdm
import daiquiri
import warnings

import pandas as pd

# from enum import Enum
from multiprocessing import Pool


from omega import __logger_name__, __version__
from omega.src.estimator.assemble import Grouping, Assembler
from omega.src.estimator.omega import bayes_infer, mle_infer

logger = daiquiri.getLogger(__logger_name__ + '.estimator')

warnings.filterwarnings(module='tensorflow*', action='ignore')


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

    ground_control = Assembler(depths, vep, mut_counts, mutability, group)
    return list(ground_control.input_generator())
    

def bayes(input_json, output_fn, cores=4):
    logger.info("Running in bayes mode")

    input_grid = prepare_data(input_json)
    logger.info("Data prepared")
    res = {}
    # TODO fix that it only works when using 4 cores
    with Pool(4) as p:
        for d in tqdm.tqdm(p.imap(bayes_infer, input_grid), total=len(input_grid)):
            res = {k: res.get(k, []) + d.get(k, []) for k in d}
    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)


def mle(input_json: str, output_fn: str, cores=4):
    logger.info("Running in mle mode")

    input_grid = prepare_data(input_json)
    logger.info("Data prepared")

    res = {}
    res_learning_curve = {}
    # TODO fix that it only works when using 4 cores
    # with Pool(cores) as p:
    #     for d, d_learning_curve in tqdm.tqdm(p.imap(mle_infer, input_grid), total=len(input_grid)):
    #         res = {k: res.get(k, []) + d.get(k, []) for k in d}
    #         res_learning_curve = {k: res_learning_curve.get(k, []) + d_learning_curve.get(k, []) for k in d_learning_curve}
    for args in tqdm.tqdm(input_grid):
        d, d_learning_curve = mle_infer(args)
        res = {k: res.get(k, []) + d.get(k, []) for k in d}
        res_learning_curve = {k: res_learning_curve.get(k, []) + d_learning_curve.get(k, []) for k in d_learning_curve}

    df = pd.DataFrame(res)
    df.to_csv(output_fn, sep='\t', index=False)

    df_learning_curve_tmp = pd.DataFrame(res_learning_curve)
    learning_curve_values = pd.DataFrame(df_learning_curve_tmp["learning_curve"].tolist())
    learning_curve_values.columns = list(range(1,1001))
    df_learning_curve = pd.concat((df_learning_curve_tmp[['gene','sample','impact']], learning_curve_values), axis = 1)
    df_learning_curve.to_csv(f"{output_fn}.log", sep='\t', index=False)



# class ModelType(str, Enum):
#     bayes = "bayes"
#     mle = "mle"


# def run(input_json: str, output_fn: str, option: ModelType=ModelType.bayes, cores=4):
#     with open(input_json, 'rt') as f:
#         d = json.load(f)

#     if option == 'bayes':
#         bayes(d, output_fn, cores=cores)
#     if option == 'mle':
#         mle(d, output_fn, cores=cores)


def run_click(observed_mutations_file, mutability_file, depths_file, vep_annotation_file, grouping_folder, output_fn,
                    option,
                    cores):
    
    d = {'observed_mutations_file' : observed_mutations_file,
            'mutability_file': mutability_file,
            'depths_file': depths_file,
            'vep_annotation_file': vep_annotation_file,
            'grouping_folder' : grouping_folder
            }

    if option == 'bayes':
        bayes(d, output_fn, cores= int(cores))
    if option == 'mle':
        mle(d, output_fn, cores= int(cores))



if __name__ == "__main__":

    """
    python src/estimator/main.py --option bayes --cores 2 test/input_estimation.json test/output_estimation_bayes.tsv
    python src/estimator/main.py --option mle --cores 2 test/input_estimation.json test/output_estimation_mle.tsv
    """
    run_click()
