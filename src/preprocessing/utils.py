import pandas as pd
from bgreference import hg38


def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 100
            ):
    """
    Define the TYPE of a variant
    """
    if ">" in (x["REF"] + x["ALT"]) or "<" in (x["REF"] + x["ALT"]):
        return "SV"
    
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


##########
# These are functions used for computing the contexts of a dataframe of mutations or sites
##########
def getContext_from_row(x, context_size = 3):
    """
    Use bgreference to get the context around a specific position
    """
    return hg38(x["CHROM"], x["POS"] - context_size//2, size=context_size)


def translate_context(s, mapping):
    """
    Transcribe a string using the mapping dictionary
    Used for translating the DNA trinucleotide context when
    changing G or A centered contexts to their complement
    """
    translation_table = str.maketrans(mapping)
    return s.translate(translation_table)


def curate_context(x, context_size = 3, nucl_dict = { "A":"T", "C":"G", "G":"C", "T":"A" }):
    """
    Curate trinucleotide context to ensure that the middle position is a C or a T
    """

    mid_pos = context_size // 2
    
    # if the mutated position is an A or a G invert the trinucleotide context
    if x[mid_pos] == "A" or x[mid_pos] == "G":
        upd_context = translate_context(x, nucl_dict)[::-1]
        return f"{ upd_context }"
    return f"{x}"


def getContext_from_df(df_mutations, context_size = 3):
    """
    Provide a dataframe with mutations
    
    return a pandas series with the number of mutations in each context
    """
    
    cols_of_interest = df_mutations[["CHROM","POS"]]
    df_mutations["CONTEXT"] = cols_of_interest.apply(getContext_from_row, axis = 1,
                                                    context_size = context_size
                                                    )
    
    df_mutations["CONTEXT"] = df_mutations["CONTEXT"].apply(curate_context, context_size = context_size)

    return df_mutations

def build_context_mut_simple(x, nucl_dict = { "A":"T", "C":"G", "G":"C", "T":"A" }): 
    """
    Build the context>MUT string ensuring that the mutation is in the correct format
    the context format is enforced in previous functions
    """

    if x["TYPE"] != "SNV":
        return "-"
    # since the context is already corrected to only C and T in the middle we do not need to translate here
    # we are only translating the mutation in case the reference was an A or a G
    if x['REF'] == "A" or x['REF'] == "G":
        return f"{ x['CONTEXT'] }>{ nucl_dict[x['ALT']] }"
    
    return f"{x['CONTEXT']}>{x['ALT']}"







# These functions are for annotating the mutations
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

consequence_rank_dict = { consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST) }
rank_consequence_dict = { rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST) }


# consequence_rank_dict
def get_single_annotation(annotations):
    all_consequences = annotations.split(",")
    all_consequences_ranks = map(lambda x: consequence_rank_dict[x], all_consequences)
    return rank_consequence_dict[min(all_consequences_ranks)]




def VEP_annotation_to_single_row(df_annotation,
                                    canonical_only = True):
    """
    [['MUT_ID',
     'Location', 'Allele',
     'Gene', 'Feature', 'Feature_type',
     'Consequence',
     'IMPACT', 
     'cDNA_position', 'CDS_position', 'Protein_position',
     'Amino_acids', 'Codons',
     'DISTANCE', 'STRAND',
     'Existing_variation',
     'SYMBOL', 'CANONICAL',
     'ENSP'
    ]]
    """
    
    # print(f"Initial number of rows:\t{df_annotation.shape}")
    df_annotation = df_annotation.drop_duplicates().reset_index(drop = True)
    # print(f"Initial number without duplicates:\t{df_annotation.shape}")
    
    
    # update the first column name to ID
    df_annotation.columns = ["MUT_ID"] + list(df_annotation.columns[1:])
    
    if canonical_only:
        df_annotation = df_annotation[df_annotation["CANONICAL"] == "YES"]
        df_annotation = df_annotation.drop("CANONICAL", axis = 1)

        # select a subset of the columns
        df_annotation_small = df_annotation[['MUT_ID',
                                            # 'Location', 'Allele',
                                            # 'Gene', 'Feature', 'Feature_type',
                                            'Consequence',
                                            'IMPACT', 
                                            # 'cDNA_position', 'CDS_position', 'Protein_position',
                                            # 'Amino_acids', 'Codons',
                                            # 'DISTANCE', 'STRAND',
                                            # 'Existing_variation',
                                            'SYMBOL', 'CANONICAL',
                                            # 'ENSP'
                                            ]]


        df_annotation_small = df_annotation_small.drop_duplicates()
        # print(f"Selecting specific columns and removing duplicates:\t{df_annotation.shape}")
    
    else:
        # select a subset of the columns
        df_annotation_small = df_annotation[['MUT_ID',
                                            # 'Location', 'Allele',
                                            # 'Gene', 'Feature', 'Feature_type',
                                            'Consequence',
                                            'IMPACT', 
                                            # 'cDNA_position', 'CDS_position', 'Protein_position',
                                            # 'Amino_acids', 'Codons',
                                            # 'DISTANCE', 'STRAND',
                                            # 'Existing_variation',
                                            'SYMBOL',
                                            # 'CANONICAL',
                                            # 'ENSP'
                                            ]]


        df_annotation_small = df_annotation_small.drop_duplicates()
        # print(f"Selecting specific columns and removing duplicates:\t{df_annotation.shape}")
    
    
    
    # add a new column containing a single consequence per variant
    df_annotation_small["Consequence"] = df_annotation_small["Consequence"].apply(get_single_annotation)

    # assign a numerical value to each consequence according to its rank in damaging consequence 
    df_annotation_small["NUM_Consequence"] = df_annotation_small["Consequence"].map(consequence_rank_dict)

    # # get the IMPACT field in a numerical scale
    # df_annotation_small["NUM_IMPACT"] = df_annotation_small["IMPACT"].map(impact2numbers)

    # sort the data by the ID and the IMPACT in numerical format.
    # The smaller the value of NUM_IMPACT the bigger the impact.
    df_annotation_small_sorted = df_annotation_small.sort_values(by = ['MUT_ID', "NUM_Consequence"],
                                                                    ascending = (True, True)
                                                                )

    df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID'
                                                                                            # , 'SYMBOL'
                                                                                            ],
                                                                                    keep='first')
    
    # print(f"Selecting row with highest impact per variant:\t{df_annotation_small_highest_impact.shape}")
    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].copy()
    returned_df = returned_df.reset_index(drop = True)
    # print(f"Selecting row with highest impact per variant:\t{returned_df.shape}")

    # update the Consequence column to containing a single consequence per variant
    returned_df["Consequence"] = returned_df["Consequence"].apply(get_single_annotation)
    
    # we return the dataframe with all the original columns of the VEP file
    return returned_df




