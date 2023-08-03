
import os
import re
import glob
import json, itertools
import pandas as pd
import numpy as np
def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 100
            ):
        
    if ">" in (x["REF"] + x["ALT"]) or "<" in (x["REF"] + x["ALT"]):
        return "SV"
​
    elif len(x["REF"]) > (len_SV_lim+1) or len(x["ALT"]) > (len_SV_lim+1) :
        return "SV"
    
    elif x["REF"] in letters and x["ALT"] in letters:
        return "SNV"
    
    elif len(x["REF"]) == len(x["ALT"]):
        return "MNV"
    
    elif x["REF"] == "-" or ( len(x["REF"]) == 1 and x["ALT"].startswith(x["REF"]) ):
        return "INSERTION"
    
    elif x["ALT"] == "-" or ( len(x["ALT"]) == 1 and x["REF"].startswith(x["ALT"]) ):
        return "DELETION"
    
    return "COMPLEX"
​
​
Generate all mutation possiblities
I take the BED file defining the regions that we target/sequence.

I generate all possible SNVs in those coordinates and then I use either dnds or EnsemblVEP to annotate the file.

This is loaded in a notebook again and then I add the information corresponding to the mutation context of each site.

This is being done in the 2023-08-02_GeneratePossibleMutations.ipynb notebook
all_possible_muts = pd.read_csv(f"/workspace/datasets/transfer/ferran_to_ferriol/omega_tests/KidneyPanel.all_SNVs.bed_panel.annotation_summary.tsv",
                                sep = "\t", header = 0)
all_possible_muts.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT']
all_possible_muts.head()
CHROM	POS	REF	ALT	MUT_ID	GENE	IMPACT	CONTEXT_MUT
0	chr1	11108180	C	A	chr1_11108180_C/A	MTOR	essential_splice	ACC>A
1	chr1	11108180	C	G	chr1_11108180_C/G	MTOR	essential_splice	ACC>G
2	chr1	11108180	C	T	chr1_11108180_C/T	MTOR	essential_splice	ACC>T
3	chr1	11108181	C	A	chr1_11108181_C/A	MTOR	missense	CCA>A
4	chr1	11108181	C	G	chr1_11108181_C/G	MTOR	missense	CCA>G
Read depth matrix
depth_dataframe = pd.read_csv("/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz",
                              header = 0, sep = "\t")
depth_dataframe.columns = ["CHROM", "POS"] + list(depth_dataframe.columns[2:])
depth_dataframe.head()
CHROM	POS	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	...	K_58_1_A_1	K_59_1_A_1	K_60_1_A_1	K_61_1_A_1	K_62_1_A_1	K_63_1_A_1	K_64_1_A_1	K_65_1_A_1	K_66_1_A_1	K_67_1_A_1
0	chr1	11107485	33407	25356	29543	1431.0	31416	17149	22168	12006	...	16599	15914	14000	17330	17743	13844	14319	20856	20046	14330
1	chr1	11107486	33382	25292	29520	1417.0	31348	17115	22139	12014	...	16570	15868	13972	17300	17707	13827	14300	20831	20035	14309
2	chr1	11107487	33417	25311	29559	1420.0	31381	17126	22151	12038	...	16592	15879	14006	17308	17710	13849	14303	20852	20037	14317
3	chr1	11107488	33645	25392	29665	1431.0	31497	17112	22334	12129	...	16713	16022	14087	17452	17869	13894	14337	21049	20217	14357
4	chr1	11107489	33718	25390	29682	1430.0	31466	17069	22370	12181	...	16732	16045	14098	17504	17918	13902	14307	21127	20299	14332
5 rows × 65 columns

Read mutations file (MAF)
maf_df = pd.read_csv(f"/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz",
                     sep = "\t",
                     header = 0)
/tmp/jobs/fcalvet/8928335/ipykernel_2508/504074552.py:1: DtypeWarning: Columns (37) have mixed types. Specify dtype option on import or set low_memory=False.
  maf_df = pd.read_csv(f"/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz",
maf_df = maf_df[maf_df["SAMPLE_ID"].isin(pd.unique(sample_information_df["IRB_subsample_id"]))].reset_index(drop = True)
maf_df.head()
CHROM	POS	REF	ALT	FILTER	INFO	FORMAT	SAMPLE	DEPTH	ALT_DEPTH	...	CANONICAL	ENSP	CHROM_ensembl	POS_ensembl	REF_ensembl	ALT_ensembl	protein_affecting	Consequence_single	Consequence_broader	CONTEXT
0	chr1	11108243	T	A	v2;pSTD	SAMPLE=A1;TYPE=SNV;DP=32564;VD=1;AF=0;BIAS=2:0...	GT:DP:VD:AD:AF:RD:ALD	0/1:32564:1:32563,1:0:14785,17778:0,1	32564	1	...	YES	ENSP00000354558	1	11108243	T	A	True	missense_variant	missense	CTT>A
1	chr1	11108922	A	<INV>	PASS	SAMPLE=A1;TYPE=INV;DP=34652;VD=1;AF=0;BIAS=0:0...	GT:DP:VD:AD:AF:RD:ALD	0/1:34652:1:0,1:0:0,0:1,0	34652	1	...	YES	ENSP00000354558	1	11108922	A	<INV>	True	splice_polypyrimidine_tract_variant	splice_region	-
2	chr1	11109301	T	A	v2;pSTD	SAMPLE=A1;TYPE=SNV;DP=32387;VD=1;AF=0;BIAS=2:0...	GT:DP:VD:AD:AF:RD:ALD	0/1:32387:1:32386,1:0:18886,13500:1,0	32387	1	...	YES	ENSP00000354558	1	11109301	T	A	True	missense_variant	missense	ATC>A
3	chr1	11109311	T	A	v2;pSTD	SAMPLE=A1;TYPE=SNV;DP=32192;VD=1;AF=0;BIAS=2:0...	GT:DP:VD:AD:AF:RD:ALD	0/1:32192:1:32190,1:0:19035,13155:1,0	32192	1	...	YES	ENSP00000354558	1	11109311	T	A	True	missense_variant	missense	CTG>A
4	chr1	11109311	T	C	v2;pSTD	SAMPLE=A1;TYPE=SNV;DP=32192;VD=1;AF=0;BIAS=2:0...	GT:DP:VD:AD:AF:RD:ALD	0/1:32192:1:32190,1:0:19035,13155:0,1	32192	1	...	YES	ENSP00000354558	1	11109311	T	C	True	missense_variant	missense	CTG>C
5 rows × 45 columns

minimal_maf = maf_df[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID']].copy()
minimal_maf["TYPE"] = minimal_maf[['REF', 'ALT']].apply(vartype, axis = 1)
print(minimal_maf.shape)
minimal_maf = minimal_maf[minimal_maf["TYPE"] == "SNV"].reset_index(drop = True)
minimal_maf = minimal_maf.drop("TYPE", axis = 1)
print(minimal_maf.shape)
minimal_maf
(13604, 6)
(11489, 5)
CHROM	POS	REF	ALT	SAMPLE_ID
0	chr1	11108243	T	A	K_5_1_A_1
1	chr1	11109301	T	A	K_5_1_A_1
2	chr1	11109311	T	A	K_5_1_A_1
3	chr1	11109311	T	C	K_5_1_A_1
4	chr1	11109695	G	C	K_5_1_A_1
...	...	...	...	...	...
11484	chr10	87931056	A	G	K_44_1_A_1
11485	chr10	87957876	C	T	K_44_1_A_1
11486	chr10	87957895	C	G	K_44_1_A_1
11487	chr10	87960967	A	G	K_44_1_A_1
11488	chr17	7676154	G	C	K_44_1_A_1
11489 rows × 5 columns

HERE STARTS THE PROCESSING
samples = list(minimal_maf["SAMPLE_ID"].unique())
minimal_maf.head()
CHROM	POS	REF	ALT	SAMPLE_ID
0	chr1	11108243	T	A	K_5_1_A_1
1	chr1	11109301	T	A	K_5_1_A_1
2	chr1	11109311	T	A	K_5_1_A_1
3	chr1	11109311	T	C	K_5_1_A_1
4	chr1	11109695	G	C	K_5_1_A_1
depth_dataframe.head()
CHROM	POS	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	...	K_58_1_A_1	K_59_1_A_1	K_60_1_A_1	K_61_1_A_1	K_62_1_A_1	K_63_1_A_1	K_64_1_A_1	K_65_1_A_1	K_66_1_A_1	K_67_1_A_1
0	chr1	11107485	33407	25356	29543	1431.0	31416	17149	22168	12006	...	16599	15914	14000	17330	17743	13844	14319	20856	20046	14330
1	chr1	11107486	33382	25292	29520	1417.0	31348	17115	22139	12014	...	16570	15868	13972	17300	17707	13827	14300	20831	20035	14309
2	chr1	11107487	33417	25311	29559	1420.0	31381	17126	22151	12038	...	16592	15879	14006	17308	17710	13849	14303	20852	20037	14317
3	chr1	11107488	33645	25392	29665	1431.0	31497	17112	22334	12129	...	16713	16022	14087	17452	17869	13894	14337	21049	20217	14357
4	chr1	11107489	33718	25390	29682	1430.0	31466	17069	22370	12181	...	16732	16045	14098	17504	17918	13902	14307	21127	20299	14332
5 rows × 65 columns

all_possible_muts.head()
CHROM	POS	REF	ALT	MUT_ID	GENE	IMPACT	CONTEXT_MUT
0	chr1	11108180	C	A	chr1_11108180_C/A	MTOR	essential_splice	ACC>A
1	chr1	11108180	C	G	chr1_11108180_C/G	MTOR	essential_splice	ACC>G
2	chr1	11108180	C	T	chr1_11108180_C/T	MTOR	essential_splice	ACC>T
3	chr1	11108181	C	A	chr1_11108181_C/A	MTOR	missense	CCA>A
4	chr1	11108181	C	G	chr1_11108181_C/G	MTOR	missense	CCA>G
binary_depth_dataframe = (depth_dataframe.set_index(["CHROM", "POS"]) > 0).astype(int).reset_index()
binary_depth_dataframe
CHROM	POS	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	...	K_58_1_A_1	K_59_1_A_1	K_60_1_A_1	K_61_1_A_1	K_62_1_A_1	K_63_1_A_1	K_64_1_A_1	K_65_1_A_1	K_66_1_A_1	K_67_1_A_1
0	chr1	11107485	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
1	chr1	11107486	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
2	chr1	11107487	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
3	chr1	11107488	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
4	chr1	11107489	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
29194	chr17	7676590	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
29195	chr17	7676591	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
29196	chr17	7676592	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
29197	chr17	7676593	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
29198	chr17	7676594	1	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
29199 rows × 65 columns

all_possible_muts_per_sample = all_possible_muts.merge(binary_depth_dataframe, on = ["CHROM", "POS"], how = "left")
all_possible_muts_per_sample
CHROM	POS	REF	ALT	MUT_ID	GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1	K_6_1_A_1	...	K_58_1_A_1	K_59_1_A_1	K_60_1_A_1	K_61_1_A_1	K_62_1_A_1	K_63_1_A_1	K_64_1_A_1	K_65_1_A_1	K_66_1_A_1	K_67_1_A_1
0	chr1	11108180	C	A	chr1_11108180_C/A	MTOR	essential_splice	ACC>A	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
1	chr1	11108180	C	G	chr1_11108180_C/G	MTOR	essential_splice	ACC>G	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
2	chr1	11108180	C	T	chr1_11108180_C/T	MTOR	essential_splice	ACC>T	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN
3	chr1	11108181	C	A	chr1_11108181_C/A	MTOR	missense	CCA>A	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
4	chr1	11108181	C	G	chr1_11108181_C/G	MTOR	missense	CCA>G	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
87994	chr3	179234297	A	G	chr3_179234297_A/G	PIK3CA	missense	ATG>C	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
87995	chr3	179234297	A	T	chr3_179234297_A/T	PIK3CA	missense	ATG>A	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
87996	chr3	179234298	T	A	chr3_179234298_T/A	PIK3CA	missense	ATC>A	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
87997	chr3	179234298	T	C	chr3_179234298_T/C	PIK3CA	synonymous	ATC>C	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
87998	chr3	179234298	T	G	chr3_179234298_T/G	PIK3CA	missense	ATC>G	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
87999 rows × 71 columns

# wide format
muts_per_gene_impact_context_sample_wide = all_possible_muts_per_sample.groupby(by = ["GENE", "IMPACT", "CONTEXT_MUT"])[samples].sum().reset_index()
muts_per_gene_impact_context_sample_wide
GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
0	-	non_genic_variant	ACA>A	0.0	0.0	0.0	0.0	0.0	0.0	0.0	...	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
1	-	non_genic_variant	ACA>G	0.0	0.0	0.0	0.0	0.0	0.0	0.0	...	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
2	-	non_genic_variant	ACA>T	0.0	0.0	0.0	0.0	0.0	0.0	0.0	...	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
3	-	non_genic_variant	ACT>A	0.0	0.0	0.0	0.0	0.0	0.0	0.0	...	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
4	-	non_genic_variant	ACT>G	0.0	0.0	0.0	0.0	0.0	0.0	0.0	...	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
2305	VHL	synonymous	TTG>C	7.0	7.0	7.0	7.0	7.0	7.0	7.0	...	7.0	7.0	7.0	7.0	7.0	7.0	7.0	7.0	7.0	7.0
2306	VHL	synonymous	TTG>G	4.0	4.0	4.0	4.0	4.0	4.0	4.0	...	4.0	4.0	4.0	4.0	4.0	4.0	4.0	4.0	4.0	4.0
2307	VHL	synonymous	TTT>A	2.0	2.0	2.0	2.0	2.0	2.0	2.0	...	2.0	2.0	2.0	2.0	2.0	2.0	2.0	2.0	2.0	2.0
2308	VHL	synonymous	TTT>C	3.0	3.0	3.0	3.0	3.0	3.0	3.0	...	3.0	3.0	3.0	3.0	3.0	3.0	3.0	3.0	3.0	3.0
2309	VHL	synonymous	TTT>G	1.0	1.0	1.0	1.0	1.0	1.0	1.0	...	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0
2310 rows × 43 columns

muts_per_gene_impact_context_sample_long = muts_per_gene_impact_context_sample_wide.melt(id_vars = ["GENE", "IMPACT", "CONTEXT_MUT"],
                                                                                         var_name = "SAMPLE_ID",
                                                                                         value_name= "COUNT")
muts_per_gene_impact_context_sample_long
GENE	IMPACT	CONTEXT_MUT	SAMPLE_ID	COUNT
0	-	non_genic_variant	ACA>A	K_5_1_A_1	0.0
1	-	non_genic_variant	ACA>G	K_5_1_A_1	0.0
2	-	non_genic_variant	ACA>T	K_5_1_A_1	0.0
3	-	non_genic_variant	ACT>A	K_5_1_A_1	0.0
4	-	non_genic_variant	ACT>G	K_5_1_A_1	0.0
...	...	...	...	...	...
92395	VHL	synonymous	TTG>C	K_44_1_A_1	7.0
92396	VHL	synonymous	TTG>G	K_44_1_A_1	4.0
92397	VHL	synonymous	TTT>A	K_44_1_A_1	2.0
92398	VHL	synonymous	TTT>C	K_44_1_A_1	3.0
92399	VHL	synonymous	TTT>G	K_44_1_A_1	1.0
92400 rows × 5 columns

muts_per_gene_impact_context_sample_long["GENE"].value_counts()
PBRM1        12640
MTOR         12280
ARID1A       11800
SETD2        11720
BAP1         11000
TP53         10800
PTEN         10320
VHL           9440
PIK3CA         960
GNL3           600
-              480
MTOR-AS1       280
DNAH1           40
RNU6-856P       40
Name: GENE, dtype: int64
Annotate mutations
annotated_minimal_maf = minimal_maf.merge(all_possible_muts, on = ["CHROM", "POS", "REF", "ALT"], how = "left")
annotated_minimal_maf
CHROM	POS	REF	ALT	SAMPLE_ID	MUT_ID	GENE	IMPACT	CONTEXT_MUT
0	chr1	11108243	T	A	K_5_1_A_1	chr1_11108243_T/A	MTOR	missense	CTT>A
1	chr1	11109301	T	A	K_5_1_A_1	chr1_11109301_T/A	MTOR	missense	ATC>A
2	chr1	11109311	T	A	K_5_1_A_1	chr1_11109311_T/A	MTOR	missense	CTG>A
3	chr1	11109311	T	C	K_5_1_A_1	chr1_11109311_T/C	MTOR	missense	CTG>C
4	chr1	11109695	G	C	K_5_1_A_1	chr1_11109695_G/C	MTOR	synonymous	CCC>G
...	...	...	...	...	...	...	...	...	...
11484	chr10	87931056	A	G	K_44_1_A_1	chr10_87931056_A/G	PTEN	missense	CTT>C
11485	chr10	87957876	C	T	K_44_1_A_1	chr10_87957876_C/T	PTEN	synonymous	GCT>T
11486	chr10	87957895	C	G	K_44_1_A_1	chr10_87957895_C/G	PTEN	missense	TCC>G
11487	chr10	87960967	A	G	K_44_1_A_1	chr10_87960967_A/G	PTEN	missense	ATT>C
11488	chr17	7676154	G	C	K_44_1_A_1	chr17_7676154_G/C	TP53	missense	CCC>G
11489 rows × 9 columns

# wide format
obs_muts_per_gene_impact_context_sample_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "GENE", "IMPACT"])["MUT_ID"].count()
obs_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long.reset_index()
obs_muts_per_gene_impact_context_sample_long.columns = list(obs_muts_per_gene_impact_context_sample_long.columns[:-1]) + ["COUNT"]
obs_muts_per_gene_impact_context_sample_long
SAMPLE_ID	GENE	IMPACT	COUNT
0	K_10_1_A_1	ARID1A	missense	63
1	K_10_1_A_1	ARID1A	nonsense	2
2	K_10_1_A_1	ARID1A	synonymous	17
3	K_10_1_A_1	BAP1	missense	16
4	K_10_1_A_1	BAP1	synonymous	11
...	...	...	...	...
936	K_9_1_A_1	TP53	splice_region	1
937	K_9_1_A_1	TP53	synonymous	5
938	K_9_1_A_1	VHL	missense	11
939	K_9_1_A_1	VHL	nonsense	1
940	K_9_1_A_1	VHL	synonymous	3
941 rows × 4 columns

obs_syn_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long[
                                                obs_muts_per_gene_impact_context_sample_long["IMPACT"] == "synonymous"].reset_index(
                                                                                                                        drop = True)
obs_syn_muts_per_gene_impact_context_sample_long
SAMPLE_ID	GENE	IMPACT	COUNT
0	K_10_1_A_1	ARID1A	synonymous	17
1	K_10_1_A_1	BAP1	synonymous	11
2	K_10_1_A_1	MTOR	synonymous	7
3	K_10_1_A_1	PBRM1	synonymous	12
4	K_10_1_A_1	PTEN	synonymous	1
...	...	...	...	...
292	K_9_1_A_1	PBRM1	synonymous	16
293	K_9_1_A_1	PTEN	synonymous	8
294	K_9_1_A_1	SETD2	synonymous	54
295	K_9_1_A_1	TP53	synonymous	5
296	K_9_1_A_1	VHL	synonymous	3
297 rows × 4 columns

obs_syn_muts_per_gene_impact_context_sample_wide = obs_syn_muts_per_gene_impact_context_sample_long.pivot(
                                                            index='GENE', columns='SAMPLE_ID', values='COUNT').fillna(0).astype(int).reset_index()
obs_syn_muts_per_gene_impact_context_sample_wide.columns.name = None
obs_syn_muts_per_gene_impact_context_sample_wide
GENE	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	K_14_1_A_1	K_15_1_A_1	K_16_1_A_1	K_17_1_A_1	K_18_1_A_1	...	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1
0	ARID1A	17	15	33	11	27	25	19	12	25	...	11	11	8	6	3	32	18	15	4	31
1	BAP1	11	7	7	5	15	3	7	6	4	...	3	2	2	2	1	18	6	7	2	7
2	MTOR	7	11	17	9	25	12	16	13	11	...	7	7	9	5	5	61	18	15	6	23
3	PBRM1	12	10	14	15	11	12	11	11	8	...	10	14	8	6	5	22	17	15	6	16
4	PIK3CA	0	0	0	0	0	0	0	0	0	...	0	0	0	0	0	0	0	0	0	0
5	PTEN	1	4	4	1	5	2	3	4	3	...	1	1	0	0	1	6	1	3	0	8
6	SETD2	18	18	33	11	25	16	20	24	23	...	15	22	7	6	7	48	31	33	5	54
7	TP53	4	2	4	6	5	4	5	2	2	...	0	0	1	2	0	6	3	4	2	5
8	VHL	0	1	3	2	0	2	4	0	2	...	1	2	0	1	0	5	1	2	1	3
9 rows × 41 columns

​
Load mutation probabilities per sample
IDEALLY I SHOULD COMPUTE THEM HERE
Remember to add some pseudocounts to the computation
CONTEXT_MUT
mut_probability = pd.DataFrame()
mut_probability["CONTEXT_MUT"] = contexts_formatted
​
for sample in samples:
    
    # Open file for reading
#    with open(f'{mut_probability_dir}/{sample}.norm_profile.kidney_bed.json', 'r') as f:
#    with open(f'{mut_probability_dir}/{sample}.norm_profile.trinuc.json', 'r') as f:
    with open(f'{mut_probability_dir}/{sample}.norm_profile.trinuc.depth.json', 'r') as f:
        # Load JSON data
        mut_probability_dict = json.load(f)
        
    mut_probability[sample] = [mut_probability_dict[t] for t in contexts_formatted]
    
    # break
mut_probability.head()
CONTEXT_MUT	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
0	ACA>A	0.007204	0.003977	0.011025	0.007215	0.006789	0.002810	0.021650	0.009884	0.015882	...	0.010839	0.010135	0.011535	0.029379	0.000000	0.000000	0.007925	0.006432	0.009754	0.009622
1	ACC>A	0.005383	0.017462	0.009137	0.008406	0.008465	0.003464	0.015257	0.007971	0.023420	...	0.013701	0.021478	0.010860	0.009140	0.006619	0.014541	0.014927	0.032274	0.011780	0.000000
2	ACG>A	0.011189	0.000000	0.014459	0.025754	0.005364	0.022310	0.024060	0.012364	0.000000	...	0.029206	0.013656	0.000000	0.028733	0.042187	0.015060	0.015896	0.050635	0.000000	0.039434
3	ACT>A	0.000000	0.018742	0.006462	0.000000	0.006377	0.009941	0.000000	0.005842	0.007474	...	0.008499	0.003953	0.000000	0.017133	0.000000	0.009042	0.004640	0.007500	0.000000	0.011335
4	CCA>A	0.004949	0.006892	0.003782	0.009004	0.003735	0.009543	0.006307	0.014221	0.012882	...	0.010121	0.007155	0.006014	0.005054	0.010999	0.002678	0.019286	0.008932	0.000000	0.013645
5 rows × 41 columns

Mutability computation
Get synonymous counts
syn_muts_per_gene_impact_context_sample = muts_per_gene_impact_context_sample_wide[
                                                muts_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"].reset_index(
                                                                                                                        drop = True)
syn_muts_per_gene_impact_context_sample[samples] = syn_muts_per_gene_impact_context_sample[samples].fillna(0).astype(int)
syn_muts_per_gene_impact_context_sample
GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
0	ARID1A	synonymous	ACA>A	11	11	11	11	11	11	11	...	11	11	11	11	11	11	11	11	11	11
1	ARID1A	synonymous	ACA>G	11	11	11	11	11	11	11	...	11	11	11	11	11	11	11	11	11	11
2	ARID1A	synonymous	ACA>T	66	66	66	66	66	66	66	...	66	66	66	66	66	66	66	66	66	66
3	ARID1A	synonymous	ACC>A	9	9	9	9	9	9	9	...	9	9	9	9	9	9	9	9	9	9
4	ARID1A	synonymous	ACC>G	9	9	9	9	9	9	9	...	9	9	9	9	9	9	9	9	9	9
...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
709	VHL	synonymous	TTG>C	7	7	7	7	7	7	7	...	7	7	7	7	7	7	7	7	7	7
710	VHL	synonymous	TTG>G	4	4	4	4	4	4	4	...	4	4	4	4	4	4	4	4	4	4
711	VHL	synonymous	TTT>A	2	2	2	2	2	2	2	...	2	2	2	2	2	2	2	2	2	2
712	VHL	synonymous	TTT>C	3	3	3	3	3	3	3	...	3	3	3	3	3	3	3	3	3	3
713	VHL	synonymous	TTT>G	1	1	1	1	1	1	1	...	1	1	1	1	1	1	1	1	1	1
714 rows × 43 columns

syn_muts_per_gene_impact_context_sample["GENE"].value_counts()
ARID1A    92
BAP1      92
SETD2     92
PBRM1     91
MTOR      90
TP53      87
PTEN      82
VHL       80
PIK3CA     8
Name: GENE, dtype: int64
Merge synonymous sites with mutational probabilities
synonymous_sites_mut_probs = syn_muts_per_gene_impact_context_sample.merge(mut_probability,
                                                                       suffixes = [".sites", ".probability"],
                                                                       on = "CONTEXT_MUT")
synonymous_sites_mut_probs.head()
GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1.sites	K_6_1_A_1.sites	K_7_1_A_1.sites	K_8_1_A_1.sites	K_9_1_A_1.sites	K_10_1_A_1.sites	K_11_1_A_1.sites	...	K_35_1_A_1.probability	K_36_1_A_1.probability	K_37_1_A_1.probability	K_38_1_A_1.probability	K_39_1_A_1.probability	K_40_1_A_1.probability	K_41_1_A_1.probability	K_42_1_A_1.probability	K_43_1_A_1.probability	K_44_1_A_1.probability
0	ARID1A	synonymous	ACA>A	11	11	11	11	11	11	11	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
1	BAP1	synonymous	ACA>A	3	3	3	3	3	3	3	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
2	MTOR	synonymous	ACA>A	12	12	12	12	12	12	12	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
3	PBRM1	synonymous	ACA>A	7	7	7	7	7	7	7	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
4	PTEN	synonymous	ACA>A	5	5	5	5	5	5	5	...	0.010839	0.010135	0.011535	0.029379	0.0	0.0	0.007925	0.006432	0.009754	0.009622
5 rows × 83 columns

synonymous_probs_gene_context = synonymous_sites_mut_probs[["GENE", "IMPACT", "CONTEXT_MUT"]].copy()
for sample in samples:
    synonymous_probs_gene_context[sample] = synonymous_sites_mut_probs[f"{sample}.sites"] * synonymous_sites_mut_probs[f"{sample}.probability"]
synonymous_probs_gene_context
GENE	IMPACT	CONTEXT_MUT	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
0	ARID1A	synonymous	ACA>A	0.079239	0.043751	0.121278	0.079368	0.074680	0.030908	0.238154	...	0.119234	0.111484	0.126884	0.323166	0.000000	0.000000	0.087170	0.070749	0.107296	0.105837
1	BAP1	synonymous	ACA>A	0.021611	0.011932	0.033076	0.021646	0.020367	0.008430	0.064951	...	0.032518	0.030405	0.034605	0.088136	0.000000	0.000000	0.023774	0.019295	0.029262	0.028865
2	MTOR	synonymous	ACA>A	0.086443	0.047728	0.132304	0.086583	0.081470	0.033718	0.259804	...	0.130074	0.121618	0.138419	0.352545	0.000000	0.000000	0.095095	0.077180	0.117050	0.115459
3	PBRM1	synonymous	ACA>A	0.050425	0.027841	0.077177	0.050507	0.047524	0.019669	0.151552	...	0.075876	0.070944	0.080745	0.205651	0.000000	0.000000	0.055472	0.045022	0.068279	0.067351
4	PTEN	synonymous	ACA>A	0.036018	0.019887	0.055127	0.036076	0.033946	0.014049	0.108252	...	0.054197	0.050674	0.057675	0.146894	0.000000	0.000000	0.039623	0.032158	0.048771	0.048108
...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
709	MTOR	synonymous	TTT>G	0.004681	0.013094	0.024186	0.000000	0.022425	0.000000	0.021141	...	0.023506	0.000000	0.037481	0.000000	0.017138	0.036666	0.025769	0.020620	0.000000	0.000000
710	PBRM1	synonymous	TTT>G	0.016383	0.045829	0.084651	0.000000	0.078489	0.000000	0.073995	...	0.082271	0.000000	0.131184	0.000000	0.059982	0.128332	0.090191	0.072169	0.000000	0.000000
711	SETD2	synonymous	TTT>G	0.014042	0.039282	0.072558	0.000000	0.067276	0.000000	0.063424	...	0.070518	0.000000	0.112444	0.000000	0.051413	0.109999	0.077307	0.061859	0.000000	0.000000
712	TP53	synonymous	TTT>G	0.001170	0.003273	0.006047	0.000000	0.005606	0.000000	0.005285	...	0.005877	0.000000	0.009370	0.000000	0.004284	0.009167	0.006442	0.005155	0.000000	0.000000
713	VHL	synonymous	TTT>G	0.001170	0.003273	0.006047	0.000000	0.005606	0.000000	0.005285	...	0.005877	0.000000	0.009370	0.000000	0.004284	0.009167	0.006442	0.005155	0.000000	0.000000
714 rows × 43 columns

synonymous_probs_gene_context["GENE"].value_counts()
ARID1A    92
BAP1      92
SETD2     92
PBRM1     91
MTOR      90
TP53      87
PTEN      82
VHL       80
PIK3CA     8
Name: GENE, dtype: int64
# these are the values that should be compared to the observed mutations,
#    to then adjust the mutability profile
expected_and_prob_computed_syn = synonymous_probs_gene_context.drop(["IMPACT", "CONTEXT_MUT"],
                                                                      axis = 1).groupby("GENE").sum().reset_index()
expected_and_prob_computed_syn
GENE	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	...	K_35_1_A_1	K_36_1_A_1	K_37_1_A_1	K_38_1_A_1	K_39_1_A_1	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1
0	ARID1A	54.793706	52.917135	50.109266	57.572528	47.963827	53.271773	55.433878	52.263707	49.920224	...	56.944989	55.439799	54.272710	48.890340	56.173718	53.130103	55.251997	49.777512	54.313869	56.310034
1	BAP1	16.773881	16.170691	16.017419	16.723305	14.190383	16.052720	17.128361	16.985088	15.132463	...	17.302510	16.406146	16.860433	14.180926	15.076421	16.630995	17.027691	15.289855	15.412378	17.836026
2	MTOR	33.127952	32.106874	30.822802	32.969478	28.638224	32.399202	33.218407	33.727436	30.148829	...	34.333079	32.347531	32.490993	29.618702	29.072793	31.547714	32.565024	29.987956	30.232525	33.530165
3	PBRM1	32.418355	27.718065	28.462600	28.333935	30.721493	30.575938	28.855236	30.528460	27.148238	...	29.193736	30.960471	26.773498	26.547287	24.721556	24.065745	30.628855	25.441169	27.523253	28.391990
4	PIK3CA	0.035731	0.018042	0.021542	0.050147	0.019004	0.023569	0.008799	0.021202	0.019769	...	0.017563	0.029722	0.031269	0.023145	0.026200	0.004054	0.027480	0.019302	0.038808	0.028054
5	PTEN	8.791261	6.817721	7.247144	7.066167	8.408613	7.909935	6.894515	7.404436	7.225386	...	6.798051	7.480017	6.574486	7.283108	6.329155	5.577903	7.538939	6.348524	6.929577	7.359451
6	SETD2	57.529999	45.629159	48.641280	46.831496	55.388107	52.188578	49.268098	52.233740	47.307296	...	50.005824	53.339215	44.308965	45.267539	41.134096	38.108487	52.018100	41.882543	49.721392	49.656200
7	TP53	8.712561	8.204975	8.014765	8.630891	7.538078	8.666331	8.774577	8.603508	7.994692	...	9.186853	8.575601	8.213504	7.575752	7.842056	8.093661	8.514121	7.612064	8.045224	9.090149
8	VHL	5.797906	6.137848	5.607701	6.551065	4.943952	5.277301	6.226692	5.307039	5.549250	...	5.799991	6.186401	6.415235	5.377828	6.978716	6.311090	6.067287	5.915788	6.137951	5.550146
9 rows × 41 columns

this comes from above and is the number of observed mutations in each sample
obs_syn_muts_per_gene_impact_context_sample_wide
GENE	K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	K_14_1_A_1	K_15_1_A_1	K_16_1_A_1	K_17_1_A_1	K_18_1_A_1	...	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1
0	ARID1A	17	15	33	11	27	25	19	12	25	...	11	11	8	6	3	32	18	15	4	31
1	BAP1	11	7	7	5	15	3	7	6	4	...	3	2	2	2	1	18	6	7	2	7
2	MTOR	7	11	17	9	25	12	16	13	11	...	7	7	9	5	5	61	18	15	6	23
3	PBRM1	12	10	14	15	11	12	11	11	8	...	10	14	8	6	5	22	17	15	6	16
4	PIK3CA	0	0	0	0	0	0	0	0	0	...	0	0	0	0	0	0	0	0	0	0
5	PTEN	1	4	4	1	5	2	3	4	3	...	1	1	0	0	1	6	1	3	0	8
6	SETD2	18	18	33	11	25	16	20	24	23	...	15	22	7	6	7	48	31	33	5	54
7	TP53	4	2	4	6	5	4	5	2	2	...	0	0	1	2	0	6	3	4	2	5
8	VHL	0	1	3	2	0	2	4	0	2	...	1	2	0	1	0	5	1	2	1	3
9 rows × 41 columns

we compute the value of alpha per each gene-sample pair, by dividing the number of expected synonymous by the number of synonymous we would be generating with the original mutational profile
alpha_per_sample = obs_syn_muts_per_gene_impact_context_sample_wide.set_index("GENE").divide( expected_and_prob_computed_syn.set_index("GENE") )
alpha_per_sample
K_10_1_A_1	K_11_1_A_1	K_12_1_A_1	K_13_1_A_1	K_14_1_A_1	K_15_1_A_1	K_16_1_A_1	K_17_1_A_1	K_18_1_A_1	K_19_1_A_1	...	K_40_1_A_1	K_41_1_A_1	K_42_1_A_1	K_43_1_A_1	K_44_1_A_1	K_5_1_A_1	K_6_1_A_1	K_7_1_A_1	K_8_1_A_1	K_9_1_A_1
GENE																					
ARID1A	0.319118	0.270593	0.631413	0.220352	0.491175	0.423037	0.345675	0.233281	0.475029	0.217975	...	0.207039	0.199088	0.160715	0.110469	0.053276	0.584009	0.340154	0.299346	0.069478	0.646320
BAP1	0.685242	0.408679	0.412126	0.330415	0.839990	0.180090	0.418829	0.370798	0.242614	0.251714	...	0.180386	0.117456	0.130806	0.129766	0.056066	1.073097	0.371042	0.437024	0.119594	0.493292
MTOR	0.216055	0.331142	0.504041	0.298519	0.742295	0.370108	0.498483	0.395545	0.349431	0.367259	...	0.221886	0.214955	0.300120	0.165385	0.149119	1.841345	0.560628	0.486653	0.181987	0.803122
PBRM1	0.392465	0.346558	0.458588	0.552522	0.415940	0.427571	0.371387	0.372318	0.310242	0.367423	...	0.415528	0.457085	0.314451	0.217997	0.176106	0.678628	0.613318	0.527007	0.211760	0.520808
PIK3CA	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	...	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000
PTEN	0.126423	0.580171	0.540217	0.138401	0.817150	0.279950	0.413115	0.556311	0.492572	0.000000	...	0.179279	0.132645	0.000000	0.000000	0.135880	0.682496	0.146677	0.413956	0.000000	0.951405
SETD2	0.344903	0.365348	0.631776	0.232522	0.586490	0.342656	0.400639	0.480603	0.563413	0.493897	...	0.393613	0.422930	0.167134	0.120672	0.140969	0.834347	0.679390	0.678436	0.106766	0.974939
TP53	0.461556	0.227931	0.464927	0.750498	0.580811	0.451192	0.578537	0.236718	0.254192	0.126046	...	0.000000	0.000000	0.131370	0.248595	0.000000	0.688661	0.365632	0.499079	0.231726	0.663299
VHL	0.000000	0.160599	0.565287	0.360409	0.000000	0.277992	0.634276	0.000000	0.317066	0.302499	...	0.158451	0.329637	0.000000	0.162921	0.000000	0.862380	0.162924	0.356652	0.152647	0.606802
9 rows × 40 columns

Compute THE mutabilities
for sample in samples:
​
    mutations_in_sample = pd.DataFrame()
​
    # iterate over all genes using the name of the gene and
    #    the alpha for which we should correct the mutability
    for gen, alpha in alpha_per_sample[sample].items():
        # print(gen, alpha)
​
        # take the mutation probability computed from the mutations observed
        # and the sequencing depth of each context
        mut_probability_sample_gene = mut_probability[["CONTEXT_MUT", sample]].copy()
​
        # adjust the probability vector by the value of alpha
        #   corresponding to that particular gene in that sample
        mut_probability_sample_gene[sample] = mut_probability_sample_gene[sample] * alpha
        mut_probability_sample_gene.to_csv(f"/workspace/datasets/transfer/ferran_to_ferriol/omega_tests/mutabilities/mutability.{sample}.{gen}.tsv",
                                           header = True,
                                           index = False,
                                           sep = "\t")
​
