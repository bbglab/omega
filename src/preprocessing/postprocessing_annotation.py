import pandas as pd
import numpy as np

from utils import *


def vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file):

    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t", header = 0)

    all_possible_sites[["CHROM", "POS", "MUT" ]] = all_possible_sites.iloc[:,0].str.split("_", expand = True)
    all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split("/", expand = True)
    all_possible_sites["POS"] = all_possible_sites["POS"].astype(int)

    # TODO: Is it robust enough to use columns names here?
    all_possible_sites = all_possible_sites[['#Uploaded_variation', 'Location', 'Allele', 'Consequence', 'IMPACT',
                                            'SYMBOL', 'CHROM', 'POS', 'MUT', 'REF', 'ALT']]

    all_possible_sites["TYPE"] = all_possible_sites[["REF", "ALT"]].apply(vartype, axis = 1)
    all_possible_sites = all_possible_sites[all_possible_sites["TYPE"] == "SNV"].reset_index(drop = True)


    # work to improve this function
    annotated_variants = VEP_annotation_to_single_row(all_possible_sites, canonical_only = False)

    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader  consequence per variant
    annotated_variants["Consequence_broader"] = annotated_variants["Consequence"].map(GROUPING_CONSEQUENCE_DICT)

    # add context type to all SNVs
    # remove context from the other substitution types
    annotated_variants_context = getContext_from_df(annotated_variants)
    annotated_variants["CONTEXT"] = annotated_variants_context.apply(build_context_mut_simple, axis = 1)


    annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'SYMBOL', 'Consequence_broader', 'CONTEXT']]
    annotated_variants_reduced.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT']
    annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] ).reset_index(drop = True)
    # annotated_variants_reduced.head()


    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")
    
    # return annotated_variants_reduced


if __name__ == '__main__':
    # Input
    VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"

    # Output
    all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

    vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file)