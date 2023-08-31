import sys
import pandas as pd

from bgreference import hg38


def get_pos_in_row(x, context_size = 3):
    return hg38(x["CHROM"], x["POS"], size = 1)

def get_non_ref(l, letters = {"A", "C", "G", "T"}):
    return letters - set(l)

def to_int_if_possible(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

def generate_all_sites_4VEP(input_bedfile, output_file_with_sites): 


    positions_df = pd.read_csv(input_bedfile, sep = "\t", header = None)

    first_coord = positions_df.iloc[0,1]
    if to_int_if_possible(first_coord):    
        positions_df = positions_df.iloc[:,:3]

    # it means there is a header, and we don't want it
    else:
        positions_df = positions_df.iloc[1:,:3]

    positions_df.columns = ["CHROM", "START", "END"]
    positions_df["CHROM"] = positions_df["CHROM"].astype(str).str.replace("chr", "")
    positions_df[["START", "END"]] = positions_df[["START", "END"]].astype(int)


    positions_df["POS"] = [ list(range(x, y+1)) for x, y in positions_df[["START", "END"]].values ]
    positions_df = positions_df.explode("POS").reset_index(drop = True)
    positions_df = positions_df[["CHROM", "POS"]]

    positions_df["REF"] = positions_df.apply(get_pos_in_row, axis = 1)
    positions_df["ALT"] = positions_df["REF"].apply(get_non_ref)

    positions_df = positions_df.explode("ALT").reset_index(drop = True)

    # 1   881906    881906    T/C   +
    positions_df["MUTATION"] = positions_df["REF"].astype(str) + "/" + positions_df["ALT"].astype(str)
    positions_df["STRAND"] = "+"
    positions_df_end = positions_df[['CHROM', 'POS', 'POS', 'MUTATION', 'STRAND']]


    positions_df_end.to_csv(output_file_with_sites,
                            header = False,
                            index = False,
                            sep = "\t")


if __name__ == '__main__':
    # Input
    # input_bedfile = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
    input_bedfile = sys.argv[1]

    # Output
    # output_file_with_sites = ./test/preprocessing/KidneyPanel.sites4VEP.tsv"
    output_file_with_sites = sys.argv[2]

    generate_all_sites_4VEP(input_bedfile, output_file_with_sites)