import os
import re
import glob
import json, itertools
import pandas as pd
import numpy as np


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
consequence_rank_dict = { consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST) }
rank_consequence_dict = { rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST) }


# consequence_rank_dict
def get_single_annotation(annotations):
    all_consequences = annotations.split(",")
    all_consequences_ranks = map(lambda x: consequence_rank_dict[x], all_consequences)
    return rank_consequence_dict[min(all_consequences_ranks)]


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


GROUPING_CONSEQUENCE_DICT = {
    'transcript_ablation': 'nonsense',
    
    'splice_acceptor_variant': 'nonsense',
    'splice_donor_variant': 'nonsense',
    'stop_gained': 'nonsense',
    'frameshift_variant': 'nonsense',
    'stop_lost': 'nonsense',
    'start_lost': 'nonsense',
    'missense_variant': 'missense',
    'inframe_insertion': 'missense',
    'inframe_deletion': 'missense',
    
    'protein_altering_variant' : 'protein_altering_variant', ##
    'transcript_amplification' : 'transcript_amplification', ##
    'coding_sequence_variant': 'coding_sequence_variant', ##
    
    'splice_donor_variant': 'essential_splice',
    'splice_acceptor_variant': 'essential_splice',
    'splice_region_variant': 'essential_splice',
    'splice_region_variant': 'splice_region',
    'splice_donor_5th_base_variant': 'splice_region',
    'splice_donor_region_variant': 'splice_region',
    'splice_polypyrimidine_tract_variant': 'splice_region',
    'synonymous_variant': 'synonymous',
    'incomplete_terminal_codon_variant': 'synonymous',
    'start_retained_variant': 'synonymous',
    'stop_retained_variant': 'synonymous',
    'mature_miRNA_variant': 'non_coding_exon_region',
    '5_prime_UTR_variant': 'non_coding_exon_region',
    '3_prime_UTR_variant': 'non_coding_exon_region',
    'non_coding_transcript_exon_variant': 'non_coding_exon_region',
    'NMD_transcript_variant': 'non_coding_exon_region',
    'intron_variant': 'intron_variant',
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
}


from FunctionsContextCounts import *
from FunctionsAnnotateMutations import *
def build_context_mut_simple(x, context_size = 3, nucl_dict = { "A":"T", "C":"G", "G":"C", "T":"A" }): 
    if x["TYPE"] != "SNV":
        return "-"
    # since the context is already corrected to only C and T in the middle we do not need to translate here
    # we are only translating the mutation in case the reference was an A or a G
    if x['REF'] == "A" or x['REF'] == "G":
        return f"{ x['CONTEXT'] }>{ nucl_dict[x['ALT']] }"
    
    return f"{x['CONTEXT']}>{x['ALT']}"


all_possible_sites = pd.read_csv(f"/home/fcalvet/projects/omega/omega/tests_ferriol/KidneyPanel.sites.VEP_annotated.tsv",
                                sep = "\t", header = 0)

all_possible_sites[["CHROM", "POS", "MUT" ]] = all_possible_sites.iloc[:,0].str.split("_", expand = True)
all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split("/", expand = True)
all_possible_sites["POS"] = all_possible_sites["POS"].astype(int)


all_possible_sites = all_possible_sites[['#Uploaded_variation', 'Location', 'Allele', 'Consequence',
                                        'IMPACT', 'SYMBOL', 'CHROM', 'POS', 'MUT', 'REF', 'ALT']]

all_possible_sites["TYPE"] = all_possible_sites[["REF", "ALT"]].apply(vartype, axis = 1)
all_possible_sites = all_possible_sites[all_possible_sites["TYPE"] == "SNV"].reset_index(drop = True)


# work to improve this function
annotated_variants = VEP_annotation_to_single_row(all_possible_sites, canonical_only = False)


# add a new column containing a single consequence per variant
annotated_variants["Consequence"] = annotated_variants["Consequence"].apply(get_single_annotation)

annotated_variants["protein_affecting"] = annotated_variants[["IMPACT","Consequence"]].apply(is_protein_affecting, axis = 1)

# add a new column containing a broader  consequence per variant
annotated_variants["Consequence_broader"] = annotated_variants["Consequence"].map(GROUPING_CONSEQUENCE_DICT)

# add context type to all SNVs
# remove context from the other substitution types
annotated_variants_context = getContext_from_df(annotated_variants)
annotated_variants["CONTEXT"] = annotated_variants_context.apply(build_context_mut_simple, axis = 1)


annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT',
                                                    'MUT_ID', 'SYMBOL',
                                                    'Consequence_broader', 'CONTEXT']]
all_possible_sites.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT']
annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] ).reset_index(drop = True)
annotated_variants_reduced.head()


annotated_variants_reduced.to_csv(f"/home/fcalvet/projects/omega/KidneyPanel.all_SNVs.bed_panel.annotation_summary.tsv",
                                    header = True,
                                    index = False,
                                    sep = "\t")