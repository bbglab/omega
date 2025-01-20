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
logger = daiquiri.getLogger(__logger_name__ + '.mutabilities.assemble')



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


    def _get_region_contexts(self, gene):

        df = self.depths[self.depths['GENE'] == gene]
        return df['CONTEXT_MUT'].values
        

    def _depths_per_gene(self):

        # depths fold change relative to the mean depth per sample
        # used for mutability correction

        res = {}

        for g in self.genes:
            df = self.depths[self.depths['GENE'] == g]
            samples = [c for c in self.mutability.columns if c not in ['GENE', 'CONTEXT_MUT']]
            for s in samples:
                depths = np.nan_to_num(df[s].values.astype(np.float32))
                res[(s, g)] = depths.astype(np.float32)
        return res


    def _mutability(self):

        res = {}
        depths_dict = self._depths_per_gene()
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
                absolute_mutability_vector = mutability_sample * depths_dict[(s, g)]
                res[(s, g)] = absolute_mutability_vector.astype(np.float32)

        return res


    def _lambdas(self):
        res = {}
        mutability_dict = self._mutability()
        return mutability_dict
