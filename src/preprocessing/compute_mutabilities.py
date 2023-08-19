import os
import re
import glob
import json, itertools
import pandas as pd
import numpy as np

from utils import *

## Input
depth_dataframe_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"

mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"

all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

## Output
table_muts_x_sample_gene_impact_context = "./test/preprocessing/mutations_per_gene_impact_context.count.tsv"
mutability_path = "./test/preprocessing/mutabilities"
# check that mutability_path exists or otherwise create it


##
# Read files
##

# Read all possible mutations annotated by VEP
all_possible_sites_annotated = pd.read_csv(all_possible_sites_annotated_file, sep = "\t", header = 0)


# Read depth matrix
depth_dataframe = pd.read_csv(depth_dataframe_file, header = 0, sep = "\t")
depth_dataframe.columns = ["CHROM", "POS"] + list(depth_dataframe.columns[2:])

# Read MAF with the mutations from all the samples
# it needs to have at least these columns:
#   'CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID'
maf = pd.read_csv(mutations_file, sep = "\t", header = 0, dtype= {"CHROM" : str, "POS": int,
                                                                            "REF" : str, "ALT" : str,
                                                                            "SAMPLE_ID" : str})
minimal_maf = maf[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID']].copy()
del maf

# select only SNVs
minimal_maf["TYPE"] = minimal_maf[['REF', 'ALT']].apply(vartype, axis = 1)
minimal_maf = minimal_maf[minimal_maf["TYPE"] == "SNV"].reset_index(drop = True)
minimal_maf = minimal_maf.drop("TYPE", axis = 1)


##
# Select all samples with available data in depth and in mutations
##
samples_muts = list(minimal_maf["SAMPLE_ID"].unique())
samples_depths = list(depth_dataframe.columns[2:])
samples = sorted(list(set(samples_muts).intersection(samples_depths)))
# samples = ['K_10_1_A_1', 'K_11_1_A_1', 'K_12_1_A_1', 'K_13_1_A_1', 'K_14_1_A_1', 'K_15_1_A_1', 'K_16_1_A_1', 'K_17_1_A_1',
#            'K_18_1_A_1', 'K_19_1_A_1', 'K_20_1_A_1', 'K_21_1_A_1', 'K_22_1_A_1', 'K_23_1_A_1', 'K_24_1_A_1', 'K_25_1_A_1',
#            'K_26_1_A_1', 'K_27_1_A_1', 'K_28_1_A_1', 'K_29_1_A_1', 'K_30_1_A_1', 'K_31_1_A_1', 'K_32_1_A_1', 'K_33_1_A_1',
#            'K_34_1_A_1', 'K_35_1_A_1', 'K_36_1_A_1', 'K_37_1_A_1', 'K_38_1_A_1', 'K_39_1_A_1', 'K_40_1_A_1', 'K_41_1_A_1',
#            'K_42_1_A_1', 'K_43_1_A_1', 'K_44_1_A_1', 'K_5_1_A_1',  'K_6_1_A_1',  'K_7_1_A_1',  'K_8_1_A_1',  'K_9_1_A_1']
print(f"{len(samples)} samples maintained starting from {len(samples_muts)} samples with mutations info and {len(samples_depths)} samples with depths info.")
# print(samples)
print(f"Removed samples with mutations info: {sorted(set(samples_muts) - set(samples))}")
print(f"Removed samples with depths info: {sorted(set(samples_depths) - set(samples))}")


##
# Annotate observed mutations
##
annotated_minimal_maf = minimal_maf.merge(all_possible_sites_annotated, on = ["CHROM", "POS", "REF", "ALT"], how = "left")

##
# Count of all observed mutations per sample, gene, impact and context
##
#   Required information:
#       Annotated mutations observed
#   Output:
#       observed mutations per:
#           sample
#           gene
#           impact
#           context
##
obs_muts_per_gene_impact_context_sample_long = annotated_minimal_maf.groupby(by = ['SAMPLE_ID', "GENE", "IMPACT", "CONTEXT_MUT"])["MUT_ID"].count()
obs_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long.reset_index()
obs_muts_per_gene_impact_context_sample_long.columns = ['SAMPLE_ID', "GENE", "IMPACT", "CONTEXT_MUT", "COUNT"]

# Store table with observed mutations 
obs_muts_per_gene_impact_context_sample_long.to_csv(table_muts_x_sample_gene_impact_context,
                                                    header = True,
                                                    index = False,
                                                    sep = "\t")

# print(obs_muts_per_gene_impact_context_sample_long.head())

# only synonymous mutations wide format
obs_muts_per_gene_impact_context_sample_wide = obs_muts_per_gene_impact_context_sample_long.pivot(
                                                            index= ['GENE', "IMPACT", "CONTEXT_MUT"], columns='SAMPLE_ID', values='COUNT').fillna(0).astype(int).reset_index()
obs_muts_per_gene_impact_context_sample_wide.columns.name = None
obs_muts_per_gene_impact_context_sample_wide.to_csv(f"{table_muts_x_sample_gene_impact_context}.wide",
                                                    header = True,
                                                    index = False,
                                                    sep = "\t")


####
## HERE STARTS THE PRE-PROCESSING
####

##
# Count of all possible sites per sample (keeping only sites with enough depth)
##
#   Required information:
#       All possible sites in the panel regions
#       Depth per site per sample
#   Output:
#       sites per:
#           sample
#           gene
#           impact
#           context
##
binary_depth_dataframe = (depth_dataframe.set_index(["CHROM", "POS"]) > 0).astype(int).reset_index()
all_possible_sites_per_sample = all_possible_sites_annotated.merge(binary_depth_dataframe, on = ["CHROM", "POS"], how = "left")

# wide format
sites_per_gene_impact_context_sample_wide = all_possible_sites_per_sample.groupby(
                                                            by = ["GENE", "IMPACT", "CONTEXT_MUT"])[samples].sum().reset_index()

# long format
sites_per_gene_impact_context_sample_long = sites_per_gene_impact_context_sample_wide.melt(id_vars = ["GENE", "IMPACT", "CONTEXT_MUT"],
                                                                                            var_name = "SAMPLE_ID",
                                                                                            value_name = "COUNT")
# sites_per_gene_impact_context_sample_long



##
# Count of all observed synonymous mutations per sample and impact
##
#   Required information:
#       Annotated mutations observed
#   Output:
#       observed mutations per:
#           sample
#           gene
#           impact
##

# only synonymous mutations
obs_syn_muts_per_gene_context_sample = obs_muts_per_gene_impact_context_sample_wide[
                                                obs_muts_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"].reset_index(
                                                    drop = True)
# obs_syn_muts_per_gene_context_sample


# already in wide format
obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_context_sample.groupby(by = ["GENE"])[samples].sum()
obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_sample.reset_index()


##
# Compute mutational profile from the input data
#       ***Remember to add some pseudocounts to the computation***
##
#   Required information:
#       Annotated all possible sites
#       Annotated mutations observed
#       Depth matrix
#   Output:
#       Mutational profile per sample, computed with pseudocounts to prevent some probabilities from being 0
##


subs = [''.join(z) for z in itertools.product('CT', 'ACGT') if z[0] != z[1]]
flanks = [''.join(z) for z in itertools.product('ACGT', repeat=2)]
contexts_unformatted = sorted([(a, b) for a, b in itertools.product(subs, flanks)], key=lambda x: (x[0], x[1]))
contexts_no_change = [b[0]+a[0]+b[1] for a, b in contexts_unformatted]
contexts_formatted = [b[0]+a[0]+b[1]+'>'+a[1] for a, b in contexts_unformatted]

# create the matrix in the desired order
empty_matrix = pd.DataFrame(index = contexts_formatted)

# count the mutations per sample and per context
counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["MUT_ID"].count().reset_index()
counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "MUT_ID")
counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
# counts_x_sample_matrix = counts_x_sample_matrix.astype(int)
counts_x_sample_matrix = counts_x_sample_matrix[samples]
counts_x_sample_matrix.index.name = "CONTEXT_MUT"

pseudocount = 0.5
counts_x_sample_matrix = counts_x_sample_matrix + pseudocount

# here we have the counts matrix for the number of mutations per sample per context
# counts_x_sample_matrix


# Compute the trinucleotide depth per sample
# merge the dataframe of all possible sites with the dataframe of the depth per site per sample
trinuc_depth_per_sample = all_possible_sites_annotated.merge(depth_dataframe[["CHROM", "POS"] + samples],
                                                                on = ["CHROM", "POS"],
                                                                how = "left").groupby(by = "CONTEXT_MUT")[samples].sum()

# print(counts_x_sample_matrix.head())
# print(trinuc_depth_per_sample.head())

# divide
mut_probability = counts_x_sample_matrix.divide( trinuc_depth_per_sample )
# print(mut_probability.head())

# normalize
mut_probability = mut_probability / mut_probability.sum()
# print(mut_probability.head())

# reindex to ensure the right order
mut_probability = mut_probability.reindex(contexts_formatted)
mut_probability.index.name = "CONTEXT_MUT"
mut_probability = mut_probability.reset_index()
# mut_probability
# print(mut_probability.head())


######
# Old way of computing or loading the mutational profile
######

# mut_probability = pd.DataFrame()
# mut_probability["CONTEXT_MUT"] = contexts_formatted

# # technically we could do this as a matrix division...
# for sample in samples:
    
#     # we select the mutation counts for the first normalization
#     Y = counts_x_sample_matrix[sample].values
    
#     trinuc_counts_96 = trinuc_depth_per_sample[sample].values
    
#     # correct by the amount of times a trinucleotide appears
#     norm_profile = [count / trinuc_r for count, trinuc_r in zip(Y,
#                                                                 trinuc_counts_96)]
#     # make the vector sum to 1
#     norm_profile = list(np.array(norm_profile) / sum(norm_profile))
    
#     # turn the vector into a dictionary
#     norm_profile_dict = dict( zip(counts_x_sample_matrix.index.values, norm_profile) )
    
#     mut_probability[sample] = [norm_profile_dict[t] for t in contexts_formatted]

# Now mut_probability has the mutational profiles per sample
# mut_probability


# for sample in samples:
#     # Open file for reading
# #    with open(f'{mut_probability_dir}/{sample}.norm_profile.kidney_bed.json', 'r') as f:
# #    with open(f'{mut_probability_dir}/{sample}.norm_profile.trinuc.json', 'r') as f:
#     with open(f'{mut_probability_dir}/{sample}.norm_profile.trinuc.depth.json', 'r') as f:
#         # Load JSON data
#         mut_probability_dict = json.load(f)
#     mut_probability[sample] = [mut_probability_dict[t] for t in contexts_formatted]
#     # break
# mut_probability.head()




## Mutability computation
## Get synonymous sites counts per gene, context sample
syn_sites_per_gene_impact_context_sample = sites_per_gene_impact_context_sample_wide[
                                                    sites_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"
                                                ].reset_index(drop = True)

syn_sites_per_gene_impact_context_sample[samples] = syn_sites_per_gene_impact_context_sample[samples].fillna(0).astype(int)
# syn_sites_per_gene_impact_context_sample
# syn_sites_per_gene_impact_context_sample["GENE"].value_counts()

# the number of synonymous sites per gene should be the same for all samples
# except if there are differences in sequencing coverage of those areas

# Merge synonymous sites with mutational probabilities
synonymous_sites_mut_probs = syn_sites_per_gene_impact_context_sample.merge(mut_probability,
                                                                            suffixes = [".sites", ".probability"],
                                                                            on = "CONTEXT_MUT")
# synonymous_sites_mut_probs.head()


synonymous_probs_gene_context = synonymous_sites_mut_probs[["GENE", "IMPACT", "CONTEXT_MUT"]].copy()
for sample in samples:
    synonymous_probs_gene_context[sample] = synonymous_sites_mut_probs[f"{sample}.sites"] * synonymous_sites_mut_probs[f"{sample}.probability"]
# synonymous_probs_gene_context



# Here we add the mutation probability of all synonymous sites in each gene and sample
#       these are the expected synonymous mutations taking into account only the mutational profile
###
# these are the values that should be compared to the observed mutations,
#    to then adjust the mutability
expected_syn_per_gene_per_sample = synonymous_probs_gene_context.drop(["IMPACT", "CONTEXT_MUT"],
                                                                        axis = 1).groupby("GENE").sum().reset_index()
# print(expected_syn_per_gene_per_sample.head())


# this comes from above and is the number of observed mutations in each sample
# print(obs_syn_muts_per_gene_sample.head())




# we compute the value of alpha per each gene-sample pair,
# by dividing the number of expected synonymous by the number of synonymous we would be generating with the original mutational profile
alpha_per_sample = obs_syn_muts_per_gene_sample.set_index("GENE").divide( expected_syn_per_gene_per_sample.set_index("GENE") )
print(alpha_per_sample)



## Compute THE mutabilities
for sample in samples:

    # iterate over all genes using the name of the gene and
    #    the alpha for which we should correct the mutability
    for gen, alpha in alpha_per_sample[sample].items():
        # print(gen, alpha)

        # take the mutation probability computed from the mutations observed
        # and the sequencing depth of each context
        mut_probability_sample_gene = mut_probability[["CONTEXT_MUT", sample]].copy()

        # adjust the probability vector by the value of alpha
        #   corresponding to that particular gene in that sample
        mut_probability_sample_gene[sample] = mut_probability_sample_gene[sample] * alpha
        mut_probability_sample_gene.to_csv(f"{mutability_path}/mutability.{sample}.{gen}.tsv",
                                            header = True,
                                            index = False,
                                            sep = "\t")

