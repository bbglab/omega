import pandas as pd
from bgreference import hg38


def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 100
            ):
        
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



# These are functions used for computing the contexts or the counts of context
def getContext_from_row(x, context_size = 3):
    return hg38(x["CHROM"], x["POS"] - context_size//2, size=context_size)


def translate_context(s, mapping):
    translation_table = str.maketrans(mapping)
    return s.translate(translation_table)


def curate_context(x, context_size = 3, nucl_dict = { "A":"T", "C":"G", "G":"C", "T":"A" }):
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










# These functions are for annotating the mutations
numbers2impact = {1: 'HIGH', 2: 'MODERATE', 3: 'LOW', 4: 'MODIFIER'}
impact2numbers = {"HIGH" : 1, "MODERATE" : 2, "LOW" : 3, "MODIFIER" : 4}

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
    
    print(f"Initial number of rows:\t{df_annotation.shape}")
    df_annotation = df_annotation.drop_duplicates().reset_index(drop = True)
    print(f"Initial number without duplicates:\t{df_annotation.shape}")
    
    
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
        print(f"Selecting specific columns and removing duplicates:\t{df_annotation.shape}")
    
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
        print(f"Selecting specific columns and removing duplicates:\t{df_annotation.shape}")
    
    
    
        
    # get the IMPACT field in a numerical scale
    df_annotation_small["NUM_IMPACT"] = df_annotation_small["IMPACT"].map(impact2numbers)

    # sort the data by the ID and the IMPACT in numerical format.
    # The smaller the value of NUM_IMPACT the bigger the impact.
    df_annotation_small_sorted = df_annotation_small.sort_values(by = ['MUT_ID', "NUM_IMPACT"],
                                                                 # ascending = True,
                                                                 ascending = (True, True)
                                                                )

    df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID'
                                                                                           # , 'SYMBOL'
                                                                                           ]
                                                                                    ,
                                                                                    keep='first')
    
    print(f"Selecting row with highest impact per variant:\t{df_annotation_small_highest_impact.shape}")
    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:]
    print(f"Selecting row with highest impact per variant:\t{returned_df.shape}")
    
    # we return the dataframe with all the original columns of the VEP file
    return returned_df
    
    







