import os
import re
import glob
import json, itertools
import pandas as pd
import numpy as np

from utils import *

depth_dataframe_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"
mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"

# check which is the difference between these two
# all_possible_sites_annotated_file = "/workspace/datasets/transfer/ferran_to_ferriol/omega_tests/KidneyPanel.all_SNVs.bed_panel.annotation_summary.tsv"
all_possible_sites_annotated_file = "/workspace/datasets/transfer/ferran_to_ferriol/omega_tests/KidneyPanel.all_SNVs.bed_panel.annotation_summary2.tsv"

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
minimal_maf = pd.read_csv(mutations_file, sep = "\t", header = 0)
minimal_maf = minimal_maf[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID']]

# select only SNVs
minimal_maf["TYPE"] = minimal_maf[['REF', 'ALT']].apply(vartype, axis = 1)
minimal_maf = minimal_maf[minimal_maf["TYPE"] == "SNV"].reset_index(drop = True)
minimal_maf = minimal_maf.drop("TYPE", axis = 1)



##
# Annotate observed mutations
##
annotated_minimal_maf = minimal_maf.merge(all_possible_sites_annotated, on = ["CHROM", "POS", "REF", "ALT"], how = "left")




####
## HERE STARTS THE PRE-PROCESSING
####
samples_muts = list(minimal_maf["SAMPLE_ID"].unique())
samples_depths = list(depth_dataframe.columns[2:])
samples = list(set(samples_muts).intersection(samples_depths))
print(f"{len(samples)} samples maintained starting from {len(samples_muts)} samples starting from the mutations file and {len(samples_depths)} samples from the depths file.")

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
# Count of all observed mutations per sample
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

# long format
obs_muts_per_gene_impact_context_sample_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "GENE", "IMPACT"])["MUT_ID"].count()
obs_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long.reset_index()
obs_muts_per_gene_impact_context_sample_long.columns = list(obs_muts_per_gene_impact_context_sample_long.columns[:-1]) + ["COUNT"]
# obs_muts_per_gene_impact_context_sample_long


# only synonymous mutations
obs_syn_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long[
                                                obs_muts_per_gene_impact_context_sample_long["IMPACT"] == "synonymous"].reset_index(
                                                                                                                        drop = True)
obs_syn_muts_per_gene_impact_context_sample_long


# only synonymous mutations wide format
obs_syn_muts_per_gene_impact_context_sample_wide = obs_syn_muts_per_gene_impact_context_sample_long.pivot(
                                                            index='GENE', columns='SAMPLE_ID', values='COUNT').fillna(0).astype(int).reset_index()
obs_syn_muts_per_gene_impact_context_sample_wide.columns.name = None


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
counts_x_sample_matrix = counts_x_sample_matrix.astype(int)

pseudocount = 0.5
counts_x_sample_matrix = counts_x_sample_matrix + pseudocount

# here we have the counts matrix for the number of mutations per sample per context
# counts_x_sample_matrix


# Compute the trinucleotide depth per sample
# merge the dataframe of all possible sites with the dataframe of the depth per site per sample
trinuc_depth_per_sample = all_possible_sites_annotated.merge(depth_dataframe,
                                                    on = ["CHROM", "POS"],
                                                    how = "left").groupby(by = "CONTEXT_MUT")[samples].sum()


mut_probability = pd.DataFrame()
mut_probability["CONTEXT_MUT"] = contexts_formatted

# technically we could do this as a matrix division...
for sample in samples:
    
    # we select the mutation counts for the first normalization
    Y = counts_x_sample_matrix[sample].values
    
    trinuc_counts_96 = trinuc_depth_per_sample[sample].values
    
    # correct by the amount of times a trinucleotide appears
    norm_profile = [count / trinuc_r for count, trinuc_r in zip(Y,
                                                                trinuc_counts_96)]
    # make the vector sum to 1
    norm_profile = list(np.array(norm_profile) / sum(norm_profile))
    
    # turn the vector into a dictionary
    norm_profile_dict = dict( zip(counts_x_sample_matrix.index.values, norm_profile) )
    
    mut_probability[sample] = [norm_profile_dict[t] for t in contexts_formatted]

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
# GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1.sites	K_6_1_A_1.sites	K_7_1_A_1.sites	K_8_1_A_1.sites	K_9_1_A_1.sites	K_10_1_A_1.sites	K_11_1_A_1.sites	...	K_35_1_A_1.probability	K_36_1_A_1.probability	K_37_1_A_1.probability	K_38_1_A_1.probability	K_39_1_A_1.probability	K_40_1_A_1.probability	K_41_1_A_1.probability	K_42_1_A_1.probability	K_43_1_A_1.probability	K_44_1_A_1.probability
# 0	ARID1A	synonymous	ACA>A	11	11	11	11	11	11	11	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
# 1	BAP1	synonymous	ACA>A	3	3	3	3	3	3	3	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
# 2	MTOR	synonymous	ACA>A	12	12	12	12	12	12	12	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
# 3	PBRM1	synonymous	ACA>A	7	7	7	7	7	7	7	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
# 4	PTEN	synonymous	ACA>A	5	5	5	5	5	5	5	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
# 5 rows × 83 columns

synonymous_probs_gene_context = synonymous_sites_mut_probs[["GENE", "IMPACT", "CONTEXT_MUT"]].copy()
for sample in samples:
    synonymous_probs_gene_context[sample] = synonymous_sites_mut_probs[f"{sample}.sites"] * synonymous_sites_mut_probs[f"{sample}.probability"]
synonymous_probs_gene_context
# GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
# 0	ARID1A	synonymous	ACA>A	0.079239	0.043751	0.121278	0.079368	0.074680	0.030908	0.238154	...	0.119234	0.111484	0.126884	0.323166	0.000000	0.000000	0.087170	0.070749	0.107296	0.105837
# 1	BAP1	synonymous	ACA>A	0.021611	0.011932	0.033076	0.021646	0.020367	0.008430	0.064951	...	0.032518	0.030405	0.034605	0.088136	0.000000	0.000000	0.023774	0.019295	0.029262	0.028865
# 2	MTOR	synonymous	ACA>A	0.086443	0.047728	0.132304	0.086583	0.081470	0.033718	0.259804	...	0.130074	0.121618	0.138419	0.352545	0.000000	0.000000	0.095095	0.077180	0.117050	0.115459
# 3	PBRM1	synonymous	ACA>A	0.050425	0.027841	0.077177	0.050507	0.047524	0.019669	0.151552	...	0.075876	0.070944	0.080745	0.205651	0.000000	0.000000	0.055472	0.045022	0.068279	0.067351
# 4	PTEN	synonymous	ACA>A	0.036018	0.019887	0.055127	0.036076	0.033946	0.014049	0.108252	...	0.054197	0.050674	0.057675	0.146894	0.000000	0.000000	0.039623	0.032158	0.048771	0.048108
# ...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
# 709	MTOR	synonymous	TTT>G	0.004681	0.013094	0.024186	0.000000	0.022425	0.000000	0.021141	...	0.023506	0.000000	0.037481	0.000000	0.017138	0.036666	0.025769	0.020620	0.000000	0.000000
# 710	PBRM1	synonymous	TTT>G	0.016383	0.045829	0.084651	0.000000	0.078489	0.000000	0.073995	...	0.082271	0.000000	0.131184	0.000000	0.059982	0.128332	0.090191	0.072169	0.000000	0.000000
# 711	SETD2	synonymous	TTT>G	0.014042	0.039282	0.072558	0.000000	0.067276	0.000000	0.063424	...	0.070518	0.000000	0.112444	0.000000	0.051413	0.109999	0.077307	0.061859	0.000000	0.000000
# 712	TP53	synonymous	TTT>G	0.001170	0.003273	0.006047	0.000000	0.005606	0.000000	0.005285	...	0.005877	0.000000	0.009370	0.000000	0.004284	0.009167	0.006442	0.005155	0.000000	0.000000
# 713	VHL	synonymous	TTT>G	0.001170	0.003273	0.006047	0.000000	0.005606	0.000000	0.005285	...	0.005877	0.000000	0.009370	0.000000	0.004284	0.009167	0.006442	0.005155	0.000000	0.000000
# 714 rows × 43 columns

# synonymous_probs_gene_context["GENE"].value_counts()
# ARID1A    92
# BAP1      92
# SETD2     92
# PBRM1     91
# MTOR      90
# TP53      87
# PTEN      82
# VHL       80
# PIK3CA     8
# Name: GENE, dtype: int64
# not all genes have all contexts with synonymous mutations and with coverage



# Here we add the mutation probability of all synonymous sites in each gene and sample
#       these are the expected synonymous mutations taking into account only the mutational profile
###
# these are the values that should be compared to the observed mutations,
#    to then adjust the mutability
expected_and_prob_computed_syn = synonymous_probs_gene_context.drop(["IMPACT", "CONTEXT_MUT"],
                                                                        axis = 1).groupby("GENE").sum().reset_index()
expected_and_prob_computed_syn
# GENE	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
# 0	ARID1A	54.793706	52.917135	50.109266	57.572528	47.963827	53.271773	55.433878	52.263707	49.920224	...	56.944989	55.439799	54.272710	48.890340	56.173718	53.130103	55.251997	49.777512	54.313869	56.310034
# 1	BAP1	16.773881	16.170691	16.017419	16.723305	14.190383	16.052720	17.128361	16.985088	15.132463	...	17.302510	16.406146	16.860433	14.180926	15.076421	16.630995	17.027691	15.289855	15.412378	17.836026
# 2	MTOR	33.127952	32.106874	30.822802	32.969478	28.638224	32.399202	33.218407	33.727436	30.148829	...	34.333079	32.347531	32.490993	29.618702	29.072793	31.547714	32.565024	29.987956	30.232525	33.530165
# 3	PBRM1	32.418355	27.718065	28.462600	28.333935	30.721493	30.575938	28.855236	30.528460	27.148238	...	29.193736	30.960471	26.773498	26.547287	24.721556	24.065745	30.628855	25.441169	27.523253	28.391990
# 4	PIK3CA	0.035731	0.018042	0.021542	0.050147	0.019004	0.023569	0.008799	0.021202	0.019769	...	0.017563	0.029722	0.031269	0.023145	0.026200	0.004054	0.027480	0.019302	0.038808	0.028054
# 5	PTEN	8.791261	6.817721	7.247144	7.066167	8.408613	7.909935	6.894515	7.404436	7.225386	...	6.798051	7.480017	6.574486	7.283108	6.329155	5.577903	7.538939	6.348524	6.929577	7.359451
# 6	SETD2	57.529999	45.629159	48.641280	46.831496	55.388107	52.188578	49.268098	52.233740	47.307296	...	50.005824	53.339215	44.308965	45.267539	41.134096	38.108487	52.018100	41.882543	49.721392	49.656200
# 7	TP53	8.712561	8.204975	8.014765	8.630891	7.538078	8.666331	8.774577	8.603508	7.994692	...	9.186853	8.575601	8.213504	7.575752	7.842056	8.093661	8.514121	7.612064	8.045224	9.090149
# 8	VHL	5.797906	6.137848	5.607701	6.551065	4.943952	5.277301	6.226692	5.307039	5.549250	...	5.799991	6.186401	6.415235	5.377828	6.978716	6.311090	6.067287	5.915788	6.137951	5.550146
# 9 rows × 41 columns

# this comes from above and is the number of observed mutations in each sample
# obs_syn_muts_per_gene_impact_context_sample_wide
# # GENE	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	K_14_1_A_1	K_15_1_A_1	K_16_1_A_1	K_17_1_A_1	K_18_1_A_1	...	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1
# # 0	ARID1A	17	15	33	11	27	25	19	12	25	...	11	11	8	6	3	32	18	15	4	31
# # 1	BAP1	11	7	7	5	15	3	7	6	4	...	3	2	2	2	1	18	6	7	2	7
# # 2	MTOR	7	11	17	9	25	12	16	13	11	...	7	7	9	5	5	61	18	15	6	23
# # 3	PBRM1	12	10	14	15	11	12	11	11	8	...	10	14	8	6	5	22	17	15	6	16
# # 4	PIK3CA	0	0	0	0	0	0	0	0	0	...	0	0	0	0	0	0	0	0	0	0
# # 5	PTEN	1	4	4	1	5	2	3	4	3	...	1	1	0	0	1	6	1	3	0	8
# # 6	SETD2	18	18	33	11	25	16	20	24	23	...	15	22	7	6	7	48	31	33	5	54
# # 7	TP53	4	2	4	6	5	4	5	2	2	...	0	0	1	2	0	6	3	4	2	5
# # 8	VHL	0	1	3	2	0	2	4	0	2	...	1	2	0	1	0	5	1	2	1	3
# # 9 rows × 41 columns



# we compute the value of alpha per each gene-sample pair,
# by dividing the number of expected synonymous by the number of synonymous we would be generating with the original mutational profile

alpha_per_sample = obs_syn_muts_per_gene_impact_context_sample_wide.set_index("GENE").divide( expected_and_prob_computed_syn.set_index("GENE") )
alpha_per_sample
# K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	K_14_1_A_1	K_15_1_A_1	K_16_1_A_1	K_17_1_A_1	K_18_1_A_1	K_19_1_A_1	...	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1
# GENE																					
# ARID1A	0.319118	0.270593	0.631413	0.220352	0.491175	0.423037	0.345675	0.233281	0.475029	0.217975	...	0.207039	0.199088	0.160715	0.110469	0.053276	0.584009	0.340154	0.299346	0.069478	0.646320
# BAP1	0.685242	0.408679	0.412126	0.330415	0.839990	0.180090	0.418829	0.370798	0.242614	0.251714	...	0.180386	0.117456	0.130806	0.129766	0.056066	1.073097	0.371042	0.437024	0.119594	0.493292
# MTOR	0.216055	0.331142	0.504041	0.298519	0.742295	0.370108	0.498483	0.395545	0.349431	0.367259	...	0.221886	0.214955	0.300120	0.165385	0.149119	1.841345	0.560628	0.486653	0.181987	0.803122
# PBRM1	0.392465	0.346558	0.458588	0.552522	0.415940	0.427571	0.371387	0.372318	0.310242	0.367423	...	0.415528	0.457085	0.314451	0.217997	0.176106	0.678628	0.613318	0.527007	0.211760	0.520808
# PIK3CA	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	...	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
# PTEN	0.126423	0.580171	0.540217	0.138401	0.817150	0.279950	0.413115	0.556311	0.492572	0.000000	...	0.179279	0.132645	0.000000	0.000000	0.135880	0.682496	0.146677	0.413956	0.000000	0.951405
# SETD2	0.344903	0.365348	0.631776	0.232522	0.586490	0.342656	0.400639	0.480603	0.563413	0.493897	...	0.393613	0.422930	0.167134	0.120672	0.140969	0.834347	0.679390	0.678436	0.106766	0.974939
# TP53	0.461556	0.227931	0.464927	0.750498	0.580811	0.451192	0.578537	0.236718	0.254192	0.126046	...	0.000000	0.000000	0.131370	0.248595	0.000000	0.688661	0.365632	0.499079	0.231726	0.663299
# VHL	0.000000	0.160599	0.565287	0.360409	0.000000	0.277992	0.634276	0.000000	0.317066	0.302499	...	0.158451	0.329637	0.000000	0.162921	0.000000	0.862380	0.162924	0.356652	0.152647	0.606802
# 9 rows × 40 columns




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
        mut_probability_sample_gene.to_csv(f"/home/fcalvet/projects/omega/omega/tests_ferriol/mutabilities2/mutability.{sample}.{gen}.tsv",
                                            header = True,
                                            index = False,
                                            sep = "\t")

