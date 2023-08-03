
￼
Read Ensembl VEP annotation
Load functions
[62]:
￼
CONSEQUENCES_LIST = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
]
​
consequence_rank_dict = { consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST) }
rank_consequence_dict = { rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST) }
[63]:
￼
# consequence_rank_dict
def get_single_annotation(annotations):
    all_consequences = annotations.split(",")
    all_consequences_ranks = map(lambda x: consequence_rank_dict[x], all_consequences)
    return rank_consequence_dict[min(all_consequences_ranks)]
[64]:
￼
def is_protein_affecting(dat):
    """
    dat must have Consequence, IMPACT columns    
    """
    if dat["IMPACT"] in ["HIGH", "MODERATE"]:
        return True
    elif dat["IMPACT"] == "LOW" and (not "synonymous_variant" in dat["Consequence"]) :
        return True
    elif "coding_sequence_variant" in dat["Consequence"]:
        return True
    return False
[65]:
￼
GROUPING_CONSEQUENCE_DICT = {
​
    'transcript_ablation': 'nonsense',
    
    'splice_acceptor_variant': 'nonsense',
    'splice_donor_variant': 'nonsense',
    'stop_gained': 'nonsense',
    'frameshift_variant': 'nonsense',
    'stop_lost': 'nonsense',
    'start_lost': 'nonsense',
​
    'missense_variant': 'missense',
    'inframe_insertion': 'missense',
    'inframe_deletion': 'missense',
    
    'protein_altering_variant' : 'protein_altering_variant', ##
    'transcript_amplification' : 'transcript_amplification', ##
    'coding_sequence_variant': 'coding_sequence_variant', ##
    
    'splice_donor_variant': 'essential_splice',
    'splice_acceptor_variant': 'essential_splice',
    'splice_region_variant': 'essential_splice',
​
​
    'splice_region_variant': 'splice_region',
    'splice_donor_5th_base_variant': 'splice_region',
    'splice_donor_region_variant': 'splice_region',
    'splice_polypyrimidine_tract_variant': 'splice_region',
​
    'synonymous_variant': 'synonymous',
    'incomplete_terminal_codon_variant': 'synonymous',
    'start_retained_variant': 'synonymous',
    'stop_retained_variant': 'synonymous',
​
​
    'mature_miRNA_variant': 'non_coding_exon_region',
    '5_prime_UTR_variant': 'non_coding_exon_region',
    '3_prime_UTR_variant': 'non_coding_exon_region',
    'non_coding_transcript_exon_variant': 'non_coding_exon_region',
​
    'NMD_transcript_variant': 'non_coding_exon_region',
​
    'intron_variant': 'intron_variant',
​
    'non_coding_transcript_variant' : 'non_coding_transcript_variant',
    'upstream_gene_variant': 'non_genic_variant',
    'downstream_gene_variant': 'non_genic_variant',
    'TFBS_ablation': 'non_genic_variant',
    'TFBS_amplification': 'non_genic_variant',
    'TF_binding_site_variant': 'non_genic_variant',
    'regulatory_region_ablation': 'non_genic_variant',
    'regulatory_region_amplification': 'non_genic_variant',
    'feature_elongation': 'non_genic_variant',
    'regulatory_region_variant': 'non_genic_variant',
    'feature_truncation': 'non_genic_variant',
    'intergenic_variant': 'non_genic_variant'
​
}
[66]:
￼
from FunctionsContextCounts import *
from FunctionsAnnotateMutations import *
def build_context_mut_simple(x, context_size = 3, nucl_dict = { "A":"T", "C":"G", "G":"C", "T":"A" }): 
    if x["TYPE"] != "SNV":
        return "-"
​
    # since the context is already corrected to only C and T in the middle we do not need to translate here
    # we are only translating the mutation in case the reference was an A or a G
    if x['REF'] == "A" or x['REF'] == "G":
        return f"{ x['CONTEXT'] }>{ nucl_dict[x['ALT']] }"
    
    return f"{x['CONTEXT']}>{x['ALT']}"
[67]:
￼
all_possible_muts = pd.read_csv(f"/home/fcalvet/projects/omega/KidneyPanel.sites.VEP_annotated.tsv", sep = "\t", header = 0)
[68]:
￼
all_possible_muts[["CHROM", "POS", "MUT" ]] = all_possible_muts.iloc[:,0].str.split("_", expand = True)
all_possible_muts[["REF", "ALT"]] = all_possible_muts["MUT"].str.split("/", expand = True)
all_possible_muts["POS"] = all_possible_muts["POS"].astype(int)
[69]:
￼
all_possible_muts = all_possible_muts[['#Uploaded_variation', 'Location', 'Allele', 'Consequence',
                                       'IMPACT', 'SYMBOL', 'CHROM', 'POS', 'MUT', 'REF', 'ALT']]
[70]:
￼
all_possible_muts["TYPE"] = all_possible_muts[["REF", "ALT"]].apply(vartype, axis = 1)
all_possible_muts = all_possible_muts[all_possible_muts["TYPE"] == "SNV"].reset_index(drop = True)
The processing of the output from multiple genes per variant to a single gene per variant can be optimized
Prioritize a set of genes
[71]:
￼
annotated_variants = VEP_annotation_to_single_row(all_possible_muts, canonical_only = False)
Initial number of rows:	(143805, 12)
Initial number without duplicates:	(140850, 12)
Selecting specific columns and removing duplicates:	(140850, 12)
Selecting row with highest impact per variant:	(87999, 5)
Selecting row with highest impact per variant:	(87999, 12)
[72]:
￼
annotated_variants["Consequence"].value_counts()
[72]:
missense_variant                                                                    62139
synonymous_variant                                                                  19058
stop_gained                                                                          3681
missense_variant,splice_region_variant                                               1825
splice_region_variant,synonymous_variant                                              562
splice_donor_variant                                                                  333
stop_gained,splice_region_variant                                                     124
splice_acceptor_variant                                                                87
start_lost                                                                             63
stop_lost                                                                              55
upstream_gene_variant                                                                  19
5_prime_UTR_variant                                                                    17
regulatory_region_variant                                                              12
stop_retained_variant                                                                   8
3_prime_UTR_variant                                                                     7
splice_acceptor_variant,non_coding_transcript_variant                                   6
downstream_gene_variant                                                                 2
splice_polypyrimidine_tract_variant,intron_variant,non_coding_transcript_variant        1
Name: Consequence, dtype: int64
[73]:
￼
annotated_variants["protein_affecting"] = annotated_variants[
                                                            ["IMPACT","Consequence"]].apply(is_protein_affecting, axis = 1)
​
# add a new column containing a single consequence per variant
annotated_variants["Consequence_single"] = annotated_variants["Consequence"].apply(get_single_annotation)
​
# add a new column containing a broader  consequence per variant
annotated_variants["Consequence_broader"] = annotated_variants["Consequence_single"].map(GROUPING_CONSEQUENCE_DICT)
​
​
# add context type to all SNVs
# remove context from the other substitution types
annotated_variants_context = getContext_from_df(annotated_variants)
annotated_variants["CONTEXT"] = annotated_variants_context.apply(build_context_mut_simple, axis = 1)
​
[74]:
￼
annotated_variants.head()
[74]:
MUT_ID	Location	Allele	Consequence	IMPACT	SYMBOL	CHROM	POS	MUT	REF	ALT	TYPE	protein_affecting	Consequence_single	Consequence_broader	CONTEXT
56532	chr10_87864469_C/A	10:87864469-87864469	A	regulatory_region_variant	MODIFIER	-	chr10	87864469	C/A	C	A	SNV	False	regulatory_region_variant	non_genic_variant	ACA>A
56536	chr10_87864469_C/G	10:87864469-87864469	G	regulatory_region_variant	MODIFIER	-	chr10	87864469	C/G	C	G	SNV	False	regulatory_region_variant	non_genic_variant	ACA>G
56540	chr10_87864469_C/T	10:87864469-87864469	T	regulatory_region_variant	MODIFIER	-	chr10	87864469	C/T	C	T	SNV	False	regulatory_region_variant	non_genic_variant	ACA>T
56546	chr10_87864470_A/C	10:87864470-87864470	C	start_lost	HIGH	PTEN	chr10	87864470	A/C	A	C	SNV	True	start_lost	nonsense	ATG>G
56551	chr10_87864470_A/G	10:87864470-87864470	G	start_lost	HIGH	PTEN	chr10	87864470	A/G	A	G	SNV	True	start_lost	nonsense	ATG>C
[75]:
￼
annotated_variants["Consequence_broader"].value_counts()
[75]:
missense                  63964
synonymous                19066
nonsense                   3923
splice_region               563
essential_splice            426
non_genic_variant            33
non_coding_exon_region       24
Name: Consequence_broader, dtype: int64
[ ]:
￼
​
[76]:
￼
annotated_variants.columns
[76]:
Index(['MUT_ID', 'Location', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL',
       'CHROM', 'POS', 'MUT', 'REF', 'ALT', 'TYPE', 'protein_affecting',
       'Consequence_single', 'Consequence_broader', 'CONTEXT'],
      dtype='object')
[82]:
￼
annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT',
                                                 'MUT_ID', 'SYMBOL',
                                                 'Consequence_broader', 'CONTEXT']]
annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] ).reset_index(drop = True)
annotated_variants_reduced.head()
# chr   pos ref mut gene    impact  context_mut mut_id
[82]:
CHROM	POS	REF	ALT	MUT_ID	SYMBOL	Consequence_broader	CONTEXT
0	chr1	11108180	C	A	chr1_11108180_C/A	MTOR	essential_splice	ACC>A
1	chr1	11108180	C	G	chr1_11108180_C/G	MTOR	essential_splice	ACC>G
2	chr1	11108180	C	T	chr1_11108180_C/T	MTOR	essential_splice	ACC>T
3	chr1	11108181	C	A	chr1_11108181_C/A	MTOR	missense	CCA>A
4	chr1	11108181	C	G	chr1_11108181_C/G	MTOR	missense	CCA>G
[83]:
￼
annotated_variants_reduced.to_csv(f"/home/fcalvet/projects/omega/KidneyPanel.all_SNVs.bed_panel.annotation_summary.tsv",
                                  header = True,
                                  index = False,
                                  sep = "\t")