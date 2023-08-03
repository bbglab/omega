import os
import json
import functools
import operator

import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

from utils import canonical_channels


channels = canonical_channels()


class Configurator:

    def __init__(self):
        
        self.conf = {}

    def add_conf(self, key, fn):

        assert(key in ['samples', 'genes', 'impacts'])
        with open(fn, 'rt') as f:
            self.conf[key] = json.load(f)

    def namespace(self, key):
        
        return set(functools.reduce(operator.add, self.conf[key].values()))



class Assembler:

    def __init__(self, data, mut_counts, conf):

        self.data = data
        # mutation data in wide format
        # cols: chr | pos | ref | mut | context | gene | impact | samples (reads) | samples (mutabilities)
        
        self.conf = conf  
        # dict with groupings for sample, gene and impact terms, respectively
        # use "samples", "genes" and "impacts" as keys

        self.mut_counts = mut_counts
        # response counts per gene, impact and channel

        self.mutability = self._mutability()
        # correct mutability for differences in depth per site
        
        self.lambdas = self._lambdas()
        # sets up the lookup table of lambdas
        # dict with key = sample, gene, impact

        self.response = self._response()
        # sets up the lookup table of response counts
        # dict with key = sample, gene, impact


    def _depth_rescaling(self):

        # depths fold change relative to the mean depth per sample
        # used for mutability correction

        res = {}
        # for g in self.conf.namespace('genes'):
        for g in self.data['gene'].unique():
            df = self.data[self.data['gene'] == g]
            for s in self.conf.namespace('samples'):
                depths = np.nan_to_num(df[s].values.astype(np.float32))
                res[(s, g)] = depths.astype(np.float32) / depths.mean()
        return res


    def _mutability(self):

        res = {}

        ordered_genes = self.data['gene'].unique().tolist()
        d = self._depth_rescaling()

        for s in self.conf.namespace('samples'):
            rescaling = functools.reduce(operator.add, [list(d[(s, g)]) for g in ordered_genes])
            self.data[f'mutability_corrected_{s}'] = self.data[f'mutability_{s}'].values * np.array(rescaling)
        
        for s in self.conf.namespace('samples'):
            for g in self.conf.namespace('genes'):
                res[(s, g)] = self.data[self.data['gene'] == g][f'mutability_corrected_{s}'].values
        
        return res


    def _impact_indicators(self):

        res = {}
        for g in self.conf.namespace('genes'):
            df = self.data[self.data['gene'] == g]
            for i in self.conf.namespace('impacts'):
                res[(g, i)] = df['impact'].apply(lambda x: x == i).values.astype(np.float32)
        return res


    def _context_indicators(self):

        res = {}
        for g in self.conf.namespace('genes'):
            df = self.data[self.data['gene'] == g]
            for c in channels:
                res[(g, c)] = df['context_mut'].apply(lambda x: x == c).values.astype(np.float32)
        return res


    """
    def _context_vectors(self):

        res = {}
        for g in self.conf.namespace('genes'):
            df = self.data[self.data['gene'] == g]
            res[g] = df['context_mut'].values
        return res
    """
        

    def _lambdas(self):

        res = {}
        impact_indicator_dict = self._impact_indicators()
        context_indicator_dict = self._context_indicators()
        for g in self.conf.namespace('genes'):
            for s in self.conf.namespace('samples'):
                for i in self.conf.namespace('impacts'):
                    mutability = self.data[self.data['gene'] == g][f'mutability_{s}'].values
                    impact_indicator = impact_indicator_dict[(g, i)]
                    lambda_vector = [np.sum(mutability * impact_indicator * context_indicator_dict[(g, c)]) for c in channels]
                    res[(s, g, i)] = lambda_vector
        return res


    def _response(self):

        res = {}
        self.mut_counts.context_mut = self.mut_counts.context_mut.astype("category")
        self.mut_counts.context_mut = self.mut_counts.context_mut.cat.set_categories(channels)
        self.mut_counts.sort_values(["sample", "gene", "impact", "context_mut"], inplace=True)
        keys = set(zip(self.mut_counts['sample'].values, self.mut_counts['gene'].values, self.mut_counts['impact'].values))
        d = dict(zip(zip(self.mut_counts['sample'].values, 
                         self.mut_counts['gene'].values, 
                         self.mut_counts['impact'].values, 
                         self.mut_counts['context_mut'].values),
                    self.mut_counts['count'].values))
        
        for k in keys:
            res[k] = []
            for context in channels:
                res[k].append(d.get((*k, context), 0.))
        return res


    def input_data(self, sample_set, gene_set, impact_set):

        counts_tensors = [tf.convert_to_tensor(v, dtype=tf.float32) 
                         for (s, g, i), v in self.response.items() if (s in sample_set) and (g in gene_set) and (i in impact_set)]
        n = functools.reduce(tf.math.add, counts_tensors)

        lambda_tensors = [tf.convert_to_tensor(v, dtype=tf.float32) 
                         for (s, g, i), v in self.lambdas.items() if (s in sample_set) and (g in gene_set) and (i in impact_set)]
        l = functools.reduce(tf.math.add, lambda_tensors)

        return l, n


    def input_generator(self):

        for gene_term, gene_set in self.conf.conf['genes'].items():
            for sample_term, sample_set in self.conf.conf['samples'].items():
                for impact_term, impact_set in self.conf.conf['impacts'].items():
                    try:
                        l, n = self.input_data(sample_set, gene_set, impact_set)
                        yield (gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n)
                    except:
                        print(gene_term, sample_term, impact_term)
