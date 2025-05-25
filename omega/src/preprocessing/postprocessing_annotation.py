
import sys
import pandas as pd

# TODO
# define the context and impact stores in a single place,
# not two different files
from omega.src.preprocessing.context_store import transform_context
from omega.src.preprocessing.impact_store import GROUPING_DICT, consequence_rank_dict, most_deleterious



def VEP_annotation_to_single_row(df_annotation):
    """
    Process Ensembl VEP output to get a single consequence per gene
    Select always the most deleterious
    """

    # add a new column containing a single consequence per variant
    df_annotation["Consequence"] = df_annotation["Consequence"].apply(most_deleterious)

    # # assign a numerical value to each consequence according to its rank in damaging consequence 
    df_annotation["NUM_Consequence"] = df_annotation["Consequence"].map(consequence_rank_dict)


    # sort the data by the ID and the IMPACT in numerical format.
    # The smaller the value of NUM_IMPACT the bigger the impact.
    df_annotation_small_sorted = df_annotation.sort_values(by = ['MUT_ID', "NUM_Consequence"],
                                                                    ascending = (True, True)
                                                                )

    df_annotation_small_highest_impact = df_annotation_small_sorted.drop_duplicates(subset=['MUT_ID', 'SYMBOL'],
                                                                                    keep='first')

    # we return the dataframe with all the original columns of the VEP file
    return df_annotation.iloc[df_annotation_small_highest_impact.index.values,:].reset_index(drop = True)




def vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file):

    # TODO: Document this or add some subprocess commands to do it in bash
    # this requires a bit more of preprocessing in bash:
    # See here:
    # cd /workspace/projects/prominent/analysis/omega/data/duplexome/
    # zcat Exome.no_header.tab.gz | cut -f 1,7,18 | awk '$3!="-"' | gzip > Exome.no_header.min_columns.tab.gz
    # zcat Exome.no_header.min_columns.tab.gz | tail -n +2 | \
    #       awk -F'\t' 'BEGIN {OFS = "\t"} {split($1, a, "[_/]"); print a[1], a[2], a[3], $2, $3}' | \
    #       gzip > Exome.no_header.min_columns.processed.no_labels.tab.gz

    # You need to put the columns in this order without header
    # all_possible_sites = all_possible_sites[['CHROM', 'POS', 'REF', 'ALT', '#Uploaded_variation',
    #                                             'Consequence', 'SYMBOL', 'MUT']]

    all_possible_sites = pd.read_csv(VEP_output_file, sep = "\t",
                                        header = None,
                                        # dtype={0: str},
                                        # usecols = ['#Uploaded_variation', 'Consequence', 'SYMBOL']
                                        )
    print("all possible sites loaded")

    # all_possible_sites[["CHROM", "POS", "MUT" ]] = all_possible_sites.iloc[:,0].str.split("_", expand = True)
    # print("mut_id splitted")

    # all_possible_sites[["REF", "ALT"]] = all_possible_sites["MUT"].str.split("/", expand = True)
    # print("mutation splitted")

    # all_possible_sites["POS"] = all_possible_sites["POS"].astype(int)
    # print("POS as int")

    # all_possible_sites = all_possible_sites[['CHROM', 'POS', 'REF', 'ALT', '#Uploaded_variation',
    #                                             'Consequence', 'SYMBOL', 'MUT']]
    # print("all possible sites columns selected")
    all_possible_sites.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'Consequence', 'SYMBOL', 'MUT']

    # if the annotation already contains a single consequence per gene this function does not do much
    # but if it contains multiple variants per gene it keeps only the most deleterious
    annotated_variants = VEP_annotation_to_single_row(all_possible_sites)
    del all_possible_sites
    print("VEP to single row working")

    # TODO: agree on a consensus for these broader consequence types
    # add a new column containing a broader  consequence per variant
    annotated_variants['IMPACT'] = annotated_variants['Consequence'].map(GROUPING_DICT)
    print("Consequence to IMPACT working")


    # add context type to all SNVs
    # remove context from the other substitution types
    annotated_variants["CONTEXT"] = annotated_variants.apply(lambda x: transform_context(x["CHROM"], x["POS"], x["MUT"]), axis = 1)
    print("Context added")

    annotated_variants_reduced = annotated_variants[['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'SYMBOL', 'IMPACT', 'CONTEXT']]
    annotated_variants_reduced.columns = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'GENE', 'IMPACT', 'CONTEXT_MUT']
    annotated_variants_reduced = annotated_variants_reduced.sort_values(by = ['CHROM', 'POS', 'REF', 'ALT'] )
    print("Annotation sorted")

    annotated_variants_reduced.to_csv(all_possible_sites_annotated_file,
                                        header = True,
                                        index = False,
                                        sep = "\t")


if __name__ == '__main__':
    # Input
    # VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
    VEP_output_file = sys.argv[1]

    # Output
    # all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
    all_possible_sites_annotated_file = sys.argv[2]

    vep2summarizedannotation(VEP_output_file, all_possible_sites_annotated_file)

