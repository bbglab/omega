"""Module to group and assemble input data for the omega estimator."""

import functools
import json
import operator
import os

import daiquiri
import numpy as np
import pandas as pd

import tensorflow as tf

from omega import __logger_name__
from omega.src.utils import canonical_channels


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

LOG = daiquiri.getLogger(__logger_name__ + '.assemble')
CHANNELS = canonical_channels()


class Grouping:
    """Class to collects groupings of samples, genes, and impacts."""

    def __init__(self):
        self.group = {}

    def add_group(self, key: str, fn: str) -> None:
        """Add a group to the grouping dictionary by loading JSON data from a file."""
        if key in ['samples', 'genes', 'impacts']:
            with open(fn, 'rt', encoding='utf-8') as f:
                self.group[key] = json.load(f)
        else:
            raise ValueError(f'Unknown grouping key: {key}')

    def namespace(self, key: str) -> set:
        """Return the set of unique items in the specified group."""
        return set(functools.reduce(operator.add, self.group[key].values()))


class Assembler:
    """
    Class to assemble input data for the omega estimator.

    Attributes
    ----------
    mut_counts : pd.DataFrame
        DataFrame containing mutation counts.
    mutability : pd.DataFrame
        DataFrame containing mutability information.
    group : Grouping
        Grouping object containing sample, gene, and impact groupings.
    depths : pd.DataFrame
        DataFrame containing depth information.
    genes : list
        List of genes present in both the grouping and depths data.
    lambdas : dict
        Dictionary mapping (sample, gene, impact) to lambda vectors.
    response : dict
        Dictionary mapping (sample, gene, impact) to response count vectors.

    Parameters
    ----------
    depths : pd.DataFrame
        DataFrame containing depth information.
    vep : pd.DataFrame
        DataFrame containing VEP annotation information.
    mutability : pd.DataFrame
        DataFrame containing mutability information.
    group : Grouping
        Grouping object containing sample, gene, and impact groupings.
    mutations_counts : pd.DataFrame, optional
        DataFrame containing mutation counts (default is an empty DataFrame).
    mode : str, optional
        Mode of operation, either 'mutabilities' or 'estimator' (default is 'mutabilities').
    ignore_zero_mutations : bool, optional
        Whether to ignore cases with zero observed mutations (default is True).
    """

    VEP_COLUMNS = ['CHROM', 'POS', 'CONTEXT_MUT', 'GENE', 'IMPACT']

    def __init__(
        self,
        depths: pd.DataFrame,
        vep: pd.DataFrame,
        mutability: pd.DataFrame,
        group: Grouping,
        mutations_counts: pd.DataFrame = pd.DataFrame(),
        mode: str = 'mutabilities',
        ignore_zero_mutations: bool = True,
    ):
        self.mut_counts: pd.DataFrame = mutations_counts
        self.mutability: pd.DataFrame = mutability
        self.group: Grouping = group
        self.ignore_zero_mutations: bool = ignore_zero_mutations

        self.mutabilities_per_site: pd.DataFrame = vep[self.VEP_COLUMNS].copy()

        # ** SETUP **

        # ** step 1 : merge annotated sites with depths dataframe and create depths attribute
        self.depths: pd.DataFrame = self.merge_annotated_site_with_depth(depths, self.mutabilities_per_site)

        # ** step 2: create list of genes present in both the grouping and depths data
        self.genes = list(set(self.group.namespace('genes')) & set(self.depths['GENE'].unique()))

        LOG.info('Number of samples: %d', len(self.group.namespace('samples')))
        LOG.debug('Samples:  %s', list(self.group.namespace('samples')))

        LOG.info('Number of genes in both grouping and depths data: %d', len(self.genes))
        LOG.debug('Genes: %s', self.genes)

        if mode == 'estimator':
            LOG.info('Impacts: %s', list(self.group.namespace('impacts')))

        # ** step 3: set up lookup table of lambdas
        # dict with key = sample, gene, impact
        self.lambdas: dict = self.compute_lambdas(mode)

        # ** step 4 - only estimator mode: set up the lookup table of response counts
        # dict with key = sample, gene, impact
        match mode:
            case 'mutabilities':
                self.response: dict = dict()
            case 'estimator':
                self.response: dict = self._response()
                LOG.debug('response computed')
        
        # Store a counter for warnings
        self.skip: int = 0

    @staticmethod
    def merge_annotated_site_with_depth(depths: pd.DataFrame, vep: pd.DataFrame) -> pd.DataFrame:
        """Merge annotated sites with depths dataframe."""
        depths.columns = ['CHROM', 'POS'] + list(depths.columns[2:])

        if 'CONTEXT' in depths.columns:
            depths = depths.drop('CONTEXT', axis=1)

        depths_merge_context_impact = vep.merge(depths, on=['CHROM', 'POS'], how='left')
        LOG.debug('depths annotated')

        return depths_merge_context_impact

    def _get_region_contexts(self, gene):
        df: pd.DataFrame = self.depths[self.depths['GENE'] == gene]
        return df['CONTEXT_MUT'].values

    def _depth_rescaling(self):
        """
        Depth rescaling computation.
        
        Depths fold change relative to the mean depth per sample used for mutability correction
        """
        res: dict[tuple[str, str], np.ndarray] = {}

        for g in self.genes:
            df: pd.DataFrame = self.depths[self.depths['GENE'] == g].copy()
            samples = [c for c in self.mutability.columns if c not in ['GENE', 'CONTEXT_MUT']]
            for s in samples:
                depths: np.ndarray = np.nan_to_num(df[s].values.astype(np.float32))
                res[(s, g)] = depths.astype(np.float32) / depths.mean()
        return res

    def _mutability(self):
        res = {}
        rescaling_dict = self._depth_rescaling()
        samples = [c for c in self.mutability.columns if c not in ['GENE', 'CONTEXT_MUT']]

        for g in self.genes:
            mutability_gene = self.mutability[self.mutability['GENE'] == g]
            region_contexts = self._get_region_contexts(g)
            for s in samples:
                mutability_sample_dict = dict(zip(mutability_gene['CONTEXT_MUT'].values, mutability_gene[s].values))

                # # This introduces misleading results, if no mutability for a given context, then no mutations can be randomized there.
                # # omega preprocessing already adds some pseudocount due to the potential absence of mutations in specific trinucleotide contexts
                # #     if a 0 reaches this step it is a real 0 either of depth or of expected number of mutations.
                # for c, value in mutability_sample_dict.items():
                #     mutability_sample_dict.update({c: max(value, 1e-4)})

                # TODO: revise if this .get(x, 0) is the right way of solving this
                mutability_sample = np.array(list(map(lambda x: mutability_sample_dict.get(x, 0), region_contexts)))
                fold_change = rescaling_dict[(s, g)]
                mutability_vector = mutability_sample * fold_change
                res[(s, g)] = mutability_vector.astype(np.float32)
        return res

    def _context_indicators(self):
        context_indicator_dict = {}
        for g in self.genes:
            df: pd.DataFrame = self.depths[self.depths['GENE'] == g].copy()
            for c in CHANNELS:
                context_indicator = df['CONTEXT_MUT'].apply(lambda x: x == c).values
                context_indicator_dict[(g, c)] = context_indicator.astype(np.float32)
        return context_indicator_dict

    def _impact_indicators(self):
        res = {}
        for g in self.genes:
            df = self.depths[self.depths['GENE'] == g]
            for i in self.group.namespace('impacts'):
                res[(g, i)] = df['IMPACT'].apply(lambda x: x == i).values.astype(np.float32)
        return res

    def compute_lambdas(self, mode: str = 'estimator') -> dict:
        """Compute lambda vectors for each (sample, gene, impact) combination."""
        mutability_dict = self._mutability()

        match mode:
            case 'mutabilities':
                samples = list(self.group.namespace('samples'))

                self.mutabilities_per_site[samples] = 0.0
                for g in self.genes:
                    for s in samples:
                        mutability = mutability_dict[(s, g)]
                        self.mutabilities_per_site.loc[self.mutabilities_per_site['GENE'] == g, s] = mutability
                LOG.debug('Mutabilities per site computed')
                return dict()
            case 'estimator':
                res = {}
                impact_indicator_dict = self._impact_indicators()
                context_indicator_dict = self._context_indicators()
                for g in self.genes:
                    for s in self.group.namespace('samples'):
                        mutability = mutability_dict[(s, g)]
                        for i in self.group.namespace('impacts'):
                            impact_indicator = impact_indicator_dict[(g, i)]
                            lambda_vector = np.array(
                                [
                                    np.sum(mutability * impact_indicator * context_indicator_dict[(g, c)])
                                    for c in CHANNELS
                                ]
                            )
                            res[(s, g, i)] = lambda_vector.astype(np.float32)

                LOG.debug('Lambdas computed')
                return res
            case _:
                raise ValueError("mode must be either mutabilities or estimator")

    def _response(self):
        res = {}
        self.mut_counts['CONTEXT_MUT'] = self.mut_counts['CONTEXT_MUT'].astype('category')
        self.mut_counts['CONTEXT_MUT'] = self.mut_counts['CONTEXT_MUT'].cat.set_categories(CHANNELS)
        self.mut_counts.sort_values(['GENE', 'IMPACT', 'CONTEXT_MUT'], inplace=True)

        genes_impacts = set(zip(self.mut_counts['GENE'].values, self.mut_counts['IMPACT'].values))
        for g, i in genes_impacts:
            df = self.mut_counts[(self.mut_counts['GENE'] == g) & (self.mut_counts['IMPACT'] == i)]
            for s in self.mut_counts.columns:
                if s not in ['GENE', 'IMPACT', 'CONTEXT_MUT']:
                    d = dict(zip(df['CONTEXT_MUT'].values, df[s].values))
                    res[(s, g, i)] = np.array([d.get(c, 0.0) for c in CHANNELS]).astype(np.float32)
        return res

    def input_data(self, sample_set, gene_set, impact_set):
        try:
            counts_tensors = [
                tf.convert_to_tensor(v, dtype=tf.float32)
                for (s, g, i), v in self.response.items()
                if (s in sample_set) and (g in gene_set) and (i in impact_set)
            ]
            n = functools.reduce(tf.math.add, counts_tensors)

        except TypeError as e:
            if 'reduce() of empty iterable with no initial value' in str(e):
                LOG.debug(
                    f'No mutations found for {sample_set}, {impact_set}, {gene_set}, filling the counts with 0s.'
                )
                n = [0.0] * 96
            else:
                LOG.error(f'Unknown error {e} for {sample_set}, {impact_set}, {gene_set}.')
                raise

        # Convert n to a tensor if it's not already one
        if not isinstance(n, tf.Tensor):
            n = tf.convert_to_tensor(n, dtype=tf.float32)

        try:
            lambda_tensors = [
                tf.convert_to_tensor(v, dtype=tf.float32)
                for (s, g, i), v in self.lambdas.items()
                if (s in sample_set) and (g in gene_set) and (i in impact_set)
            ]
            l = functools.reduce(tf.math.add, lambda_tensors)

        except Exception as e:
            LOG.debug(f'Lambda tensors for {sample_set}, {impact_set}, {gene_set} found error in {e}')
            l = [0.0] * 96

        # Convert l to a tensor if it's not already one
        if not isinstance(l, tf.Tensor):
            l = tf.convert_to_tensor(l, dtype=tf.float32)

        return l, n

    def input_generator(self):
        for gene_term, gene_set in self.group.group['genes'].items():
            for sample_term, sample_set in self.group.group['samples'].items():
                for impact_term, impact_set in self.group.group['impacts'].items():
                    # LOG.debug("{}{}{}".format(sample_set, gene_set, impact_set))
                    l, n = self.input_data(sample_set, gene_set, impact_set)

                    # either because of lambdas going to 0 or because of an error and then l put to 0, skip testing that group
                    if np.all(l == 0.0):
                        self.skip += 1
                        LOG.debug(f'Lambdas are 0, we are ignoring this case {sample_set}, {impact_set}, {gene_set}.')
                        continue
                    if np.all(n == 0.0):
                        if self.ignore_zero_mutations:
                            self.skip += 1
                            LOG.debug(
                                f'There is no mutation, we are ignoring this case {sample_set}, {impact_set}, {gene_set}.'
                            )
                            continue
                        LOG.debug(
                            f'There is no mutation, but we are keeping this case {sample_set}, {impact_set}, {gene_set}.'
                        )
                    yield (gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n)

        warning_reason = 'missing mutations or zero Lambdas' if self.ignore_zero_mutations else 'zero Lambdas'
        LOG.warning('Total number of warnings for %s: %i, check logs for more details', warning_reason, self.skip)
