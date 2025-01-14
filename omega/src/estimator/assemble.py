import os
import json
import functools
import operator
import daiquiri
import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

from omega.src.estimator.context_store import canonical_channels
channels = canonical_channels()


from omega import __logger_name__, __version__
logger = daiquiri.getLogger(__logger_name__ + '.estimator.assemble')



class Grouping:

    def __init__(self):
        
        self.group = {}

    def add_group(self, key, fn):

        assert(key in ['samples', 'genes', 'impacts'])
        with open(fn, 'rt') as f:
            self.group[key] = json.load(f)

    def namespace(self, key):
        
        return set(functools.reduce(operator.add, self.group[key].values()))


class Assembler:

    def __init__(self, depths, vep, mut_counts, mutability, group):

        self.mut_counts = mut_counts
        self.mutability = mutability
        
        # dict with groupings for sample, gene and impact terms, respectively
        # uses "samples", "genes" and "impacts" as keys
        self.group = group

        # ** SETUP **

        # ** step 1 : merge annotated sites with depths dataframe
        depths.columns = ["CHROM", "POS"] + list(depths.columns[2:])
        if "CONTEXT" in depths.columns: depths = depths.drop("CONTEXT", axis = 1)

        reduced_vep = vep[["CHROM", "POS", "CONTEXT_MUT", "GENE", "IMPACT"]]
        depths_merge_context_impact = reduced_vep.merge(depths, on=["CHROM", "POS"], how = "left")
        logger.debug("depths annotated")
        del depths
        del reduced_vep


        # ** step 2: create depths attribute
        self.depths = depths_merge_context_impact

        self.genes = list(set(self.group.namespace('genes')) & set(self.depths['GENE'].unique()))

        logger.debug("Samples")
        logger.debug(self.group.namespace('samples'))
        logger.debug("Genes")
        logger.debug(self.genes)
        logger.debug("Impacts")
        logger.debug(self.group.namespace('impacts'))

        # ** step 3: set up lookup table of lambdas
        # dict with key = sample, gene, impact
        self.lambdas = self._lambdas()
        logger.debug("lambdas computed")
        
        # ** step 4: set up the lookup table of response counts
        # dict with key = sample, gene, impact
        self.response = self._response()
        logger.debug("response computed")


    def _get_region_contexts(self, gene):

        df = self.depths[self.depths['GENE'] == gene]
        return df['CONTEXT_MUT'].values
        

    def _depth_rescaling(self):

        # depths fold change relative to the mean depth per sample
        # used for mutability correction

        res = {}

        for g in self.genes:
            df = self.depths[self.depths['GENE'] == g]
            samples = [c for c in self.mutability.columns if c not in ['GENE', 'CONTEXT_MUT']]
            for s in samples:
                depths = np.nan_to_num(df[s].values.astype(np.float32))
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

                # TODO
                # revise if this .get(x, 0) is the right way of solving this
                mutability_sample = np.array(list(map(lambda x: mutability_sample_dict.get(x, 0), region_contexts)))
                fold_change = rescaling_dict[(s, g)]
                mutability_vector = mutability_sample * fold_change
                res[(s, g)] = mutability_vector.astype(np.float32)
        return res


    def _context_indicators(self):

        context_indicator_dict = {}
        for g in self.genes:
            df = self.depths[self.depths['GENE'] == g]
            for c in channels:
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


    def _lambdas(self):

        res = {}
        impact_indicator_dict = self._impact_indicators()
        context_indicator_dict = self._context_indicators()
        mutability_dict = self._mutability()
        for g in self.genes:
            for s in self.group.namespace('samples'):
                mutability = mutability_dict[(s, g)]
                for i in self.group.namespace('impacts'):
                    impact_indicator = impact_indicator_dict[(g, i)]
                    lambda_vector = np.array([np.sum(mutability * impact_indicator * context_indicator_dict[(g, c)]) for c in channels])
                    res[(s, g, i)] = lambda_vector.astype(np.float32)
        return res


    def _response(self):

        res = {}
        self.mut_counts['CONTEXT_MUT'] = self.mut_counts['CONTEXT_MUT'].astype('category')
        self.mut_counts['CONTEXT_MUT'] = self.mut_counts['CONTEXT_MUT'].cat.set_categories(channels)
        self.mut_counts.sort_values(['GENE', 'IMPACT', 'CONTEXT_MUT'], inplace=True)
        
        genes_impacts = set(zip(self.mut_counts['GENE'].values, self.mut_counts['IMPACT'].values))
        for g, i in genes_impacts:
            df = self.mut_counts[(self.mut_counts['GENE'] == g) & (self.mut_counts['IMPACT'] == i)]
            for s in self.mut_counts.columns:
                if s not in ['GENE', 'IMPACT', 'CONTEXT_MUT']:
                    d = dict(zip(df['CONTEXT_MUT'].values, df[s].values))
                    res[(s, g, i)] = np.array([d.get(c, 0.) for c in channels]).astype(np.float32)
        return res


    def input_data(self, sample_set, gene_set, impact_set):

        try:
            counts_tensors = [tf.convert_to_tensor(v, dtype=tf.float32) 
                            for (s, g, i), v in self.response.items() if (s in sample_set) and (g in gene_set) and (i in impact_set)]
            n = functools.reduce(tf.math.add, counts_tensors)

        except TypeError as e:
            if "reduce() of empty iterable with no initial value" in str(e):
                logger.warning(f"No mutations found for {sample_set}, {impact_set}, {gene_set}, filling the counts with 0s.")
                n = [0.] * 96
            else:
                logger.error(f"Unknown error {e} for {sample_set}, {impact_set}, {gene_set}.")
                raise

        # Convert n to a tensor if it's not already one
        if not isinstance(n, tf.Tensor):
            n = tf.convert_to_tensor(n, dtype=tf.float32)

        try:
            lambda_tensors = [tf.convert_to_tensor(v, dtype=tf.float32) 
                            for (s, g, i), v in self.lambdas.items() if (s in sample_set) and (g in gene_set) and (i in impact_set)]
            l = functools.reduce(tf.math.add, lambda_tensors)

        except Exception as e:
            logger.warning(f"Lambda tensors for {sample_set}, {impact_set}, {gene_set} found error in {e}")
            l = [0.] * 96

        # Convert l to a tensor if it's not already one
        if not isinstance(l, tf.Tensor):
            l = tf.convert_to_tensor(l, dtype=tf.float32)

        return l, n


    def input_generator(self):
        for gene_term, gene_set in self.group.group['genes'].items():
            for sample_term, sample_set in self.group.group['samples'].items():
                for impact_term, impact_set in self.group.group['impacts'].items():
                    # logger.debug("{}{}{}".format(sample_set, gene_set, impact_set))
                    l, n = self.input_data(sample_set, gene_set, impact_set)

                    # either because of lambdas going to 0 or because of an error and then l put to 0, skip testing that group
                    if np.all(l == 0.):
                        logger.warning(f"Lambdas are 0, we are ignoring this case {sample_set}, {impact_set}, {gene_set}.")
                        continue
                    if np.all(n == 0.):
                        logger.warning(f"There is no mutation, we are ignoring this case {sample_set}, {impact_set}, {gene_set}.")
                        continue
                    yield (gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n)
