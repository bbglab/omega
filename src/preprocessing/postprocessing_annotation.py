import pandas as pd
import numpy as np

from utils import *

# TODO
# define the context and impact stores in a single place,
# not two different files
from context_store import transform_context
from impact_store import GROUPING_DICT, consequence_rank_dict, most_deleterious



def VEP_annotation_to_single_row(df_annotation,
                                    canonical_only = True):
    """
    Process Ensembl VEP output to get a single consequence per gene
    Select always the most deleterious
    """
    # update the first column name to ID
    df_annotation.columns = ["MUT_ID"] + list(df_annotation.columns[1:])
    
    if canonical_only:
        df_annotation = df_annotation[df_annotation["CANONICAL"] == "YES"]
        df_annotation = df_annotation.drop("CANONICAL", axis = 1)

    # select a subset of the columns
    df_annotation_small = df_annotation[['MUT_ID',
                                        'Consequence',
                                        'SYMBOL']]
    df_annotation_small = df_annotation_small.drop_duplicates()


    # add a new column containing a single consequence per variant
    df_annotation_small["Consequence"] = df_annotation_small["Consequence"].apply(most_deleterious)

    # assign a numerical value to each consequence according to its rank in damaging consequence 
    df_annotation_small["NUM_Consequence"] = df_annotation_small["Consequence"].map(consequence_rank_dict)

    # sort the data by the ID and the IMPACT in numerical format.
    # The smaller the value of NUM_IMPACT the bigger the impact.
    df_annotation_small_sorted = df_annotation_small.sort_values(by = ['MUT_ID', "NUM_Consequence"],
                                                                    ascending = (True, True)
                                                                )

    df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID', 'SYMBOL'],
                                                                                    keep='first')
    returned_df = df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].copy()
    returned_df = returned_df.reset_index(drop = True)
    
    # we return the dataframe with all the original columns of the VEP file
    return returned_df







def vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file):

    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t", header = 0)

    all_possible_sites[["CHROM", "POS", "MUT" ]] = all_possible_sites.iloc[:,0].str.split("_", expand = True)
    all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split("/", expand = True)
    all_possible_sites["POS"] = all_possible_sites["POS"].astype(int)

    # TODO: Is it robust enough to use columns names here?
    all_possible_sites = all_possible_sites[['#Uploaded_variation', 'Consequence', 'SYMBOL',
                                                'CHROM', 'POS', 'REF', 'ALT', 'MUT']]

    all_possible_sites["TYPE"] = all_possible_sites[["REF", "ALT"]].apply(vartype, axis = 1)
    all_possible_sites = all_possible_sites[all_possible_sites["TYPE"] == "SNV"].reset_index(drop = True)

    # if the annotation already contains a single consequence per gene this function does not do much
    # but if it contains multiple variants per gene it keeps only the most deleterious
    annotated_variants = VEP_annotation_to_single_row(all_possible_sites, canonical_only = False)

    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader  consequence per variant
    annotated_variants['IMPACT'] = annotated_variants['Consequence'].apply(most_deleterious)
    annotated_variants['IMPACT'] = annotated_variants['IMPACT'].apply(lambda x: GROUPING_DICT[x])

    # add context type to all SNVs
    # remove context from the other substitution types
    annotated_variants["CONTEXT"] = annotated_variants.apply(lambda x: transform_context(x["CHROM"], x["POS"], x["MUT"]), axis = 1)

    annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'SYMBOL', 'IMPACT', 'CONTEXT']]
    annotated_variants_reduced.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT']
    annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] ).reset_index(drop = True)

    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")


if __name__ == '__main__':
    # Input
    VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"

    # Output
    all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

    vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file)