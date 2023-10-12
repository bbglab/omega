import os
import json
import functools
import operator

import numpy as np
import pandas as pd

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors

from context_store import canonical_channels, transform_context
from impact_store import GROUPING_DICT, most_deleterious

channels = canonical_channels()


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

    def __init__(self, depths, regions, vep, mut_counts, mutability, group):

        self.mut_counts = mut_counts
        self.mutability = mutability
        
        # dict with groupings for sample, gene and impact terms, respectively
        # uses "samples", "genes" and "impacts" as keys
        self.group = group

        # ** SETUP **

        # ** step 1: annotate genes as "ELEMENTS" in depth table
        # 1.1: create a mergable table from the regions BED file
        feature_tuple = zip(regions['CHROMOSOME'], regions['START'], regions['END'], regions['ELEMENT'])
        d = {'chr': [], 'pos':[], 'ELEMENT': []}
        for chr_, start, end, elem in feature_tuple:
            span = range(start+1, end+1)
            l = len(span)
            d['chr'] += [chr_] * l
            d['pos'] += list(span)
            d['ELEMENT'] += [elem] * l
        mergable_regions_elements = pd.DataFrame(d)
        
        # 1.2: merge with depths dataframe
        depths["chr"] = depths["chr"].astype(str).str.replace("chr","")
        mergable_regions_elements["chr"] = mergable_regions_elements["chr"].astype(str).str.replace("chr","")
        depths_merge = depths.merge(mergable_regions_elements, on=['chr', 'pos'])

        # ** step 2: 
        # triplicate and annotate depths with impacts
        # 2.1: create a mergable VEP table
        vep_mergable = pd.DataFrame(columns=['chr', 'pos', 'CONTEXT_MUT', 'ELEMENT', 'IMPACT'])
        vep_mergable['chr'], vep_mergable['pos'], vep_mergable['CONTEXT_MUT'], vep_mergable['ELEMENT'], vep_mergable['IMPACT'] = \
            zip(*vep.apply(lambda r: r['#Uploaded_variation'].split('_') + [r['SYMBOL']] + [r['Consequence']], axis=1)) 
        vep_mergable['chr'] = vep_mergable['chr'].astype(str).str.replace("chr","").astype(str)
        vep_mergable['pos'] = vep_mergable['pos'].astype(int)
        vep_mergable['CONTEXT_MUT'] = vep_mergable.apply(lambda r: transform_context(r['chr'], r['pos'], r['CONTEXT_MUT']), axis=1)
        vep_mergable['IMPACT'] = vep_mergable['IMPACT'].apply(most_deleterious)
        # 2.2: merge with depths
        depths_merge_context_impact = pd.merge(depths_merge, vep_mergable, on=['chr', 'pos', 'ELEMENT'], how='left')
        depths_merge_context_impact['IMPACT'] = depths_merge_context_impact['IMPACT'].apply(lambda x: GROUPING_DICT[x])
        
        # ** step 3: create depths attribute
        self.depths = depths_merge_context_impact

        # ** step 4: set up lookup table of lambdas
        # dict with key = sample, gene, impact
        self.lambdas = self._lambdas()
        
        # ** step 5: set up the lookup table of response counts
        # dict with key = sample, gene, impact
        self.response = self._response()


    def _get_region_contexts(self, gene):

        df = self.depths[self.depths['ELEMENT'] == gene]
        return df['CONTEXT_MUT'].values
        

    def _depth_rescaling(self):

        # depths fold change relative to the mean depth per sample
        # used for mutability correction

        res = {}
        # for g in self.group.namespace('genes'):
        for g in self.depths['ELEMENT'].unique():
            df = self.depths[self.depths['ELEMENT'] == g]
            samples = [c for c in self.mutability.columns if c not in ['GENE', 'CONTEXT_MUT']]
            for s in samples:
                depths = np.nan_to_num(df[s].values.astype(np.float32))
                res[(s, g)] = depths.astype(np.float32) / depths.mean()
        return res


    def _mutability(self):

        res = {}
        rescaling_dict = self._depth_rescaling()
        samples = [c for c in self.mutability.columns if c not in ['GENE', 'CONTEXT_MUT']]
        genes_from_rescaling = set([gen for samp, gen in rescaling_dict.keys()])
        genes = [ gene for gene in self.mutability['GENE'].unique() if gene in genes_from_rescaling ]
        for g in genes:
            mutability_gene = self.mutability[self.mutability['GENE'] == g]
            region_contexts = self._get_region_contexts(g)
            for s in samples:
                mutability_sample_dict = dict(zip(mutability_gene['CONTEXT_MUT'].values, mutability_gene[s].values))

                for c, value in mutability_sample_dict.items():
                    mutability_sample_dict.update({c: max(value, 1e-4)})
    
                mutability_sample = np.array(list(map(lambda x: mutability_sample_dict[x], region_contexts)))
                fold_change = rescaling_dict[(s, g)]
                mutability_vector = mutability_sample * fold_change
                res[(s, g)] = mutability_vector.astype(np.float32)
        return res


    def _context_indicators(self):

        context_indicator_dict = {}
        genes = self.mutability['GENE'].unique()
        for g in genes:
            df = self.depths[self.depths['ELEMENT'] == g]
            for c in channels:
                context_indicator = df['CONTEXT_MUT'].apply(lambda x: x == c).values
                context_indicator_dict[(g, c)] = context_indicator.astype(np.float32)
        return context_indicator_dict


    def _impact_indicators(self):

        res = {}
        for g in self.group.namespace('genes'):
            df = self.depths[self.depths['ELEMENT'] == g]
            for i in self.group.namespace('impacts'):
                res[(g, i)] = df['IMPACT'].apply(lambda x: x == i).values.astype(np.float32)
        return res


    def _lambdas(self):

        res = {}
        impact_indicator_dict = self._impact_indicators()
        context_indicator_dict = self._context_indicators()
        mutability_dict = self._mutability()

        for g in self.group.namespace('genes'):
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

        counts_tensors = [tf.convert_to_tensor(v, dtype=tf.float32) 
                         for (s, g, i), v in self.response.items() if (s in sample_set) and (g in gene_set) and (i in impact_set)]
        n = functools.reduce(tf.math.add, counts_tensors)

        lambda_tensors = [tf.convert_to_tensor(v, dtype=tf.float32) 
                         for (s, g, i), v in self.lambdas.items() if (s in sample_set) and (g in gene_set) and (i in impact_set)]
        l = functools.reduce(tf.math.add, lambda_tensors)

        return l, n


    def input_generator(self):

        for gene_term, gene_set in self.group.group['genes'].items():
            for sample_term, sample_set in self.group.group['samples'].items():
                for impact_term, impact_set in self.group.group['impacts'].items():
                    l, n = self.input_data(sample_set, gene_set, impact_set)
                    yield (gene_term, sample_term, impact_term, gene_set, sample_set, impact_set, l, n)
