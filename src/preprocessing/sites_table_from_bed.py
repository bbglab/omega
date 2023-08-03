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
        return int(string)
    except ValueError:
        return None


# input_bedfile = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
input_bedfile = sys.argv[1]
# output_file_with_sites = f"{positive_sel_dir}/KidneyPanel.all_SNVs.bed.panel.4VEP.tsv"
output_file_with_sites = sys.argv[2]


positions = input_bedfile
positions_df = pd.read_csv(f"{positions}", sep = "\t", header = None)

first_coord = positions_df.iloc[0,1]
if to_int_if_possible(first_coord):    
	positions_df = positions_df.iloc[:,:3]

# it means there is a header, and we don't want it
else:
    positions_df = positions_df.iloc[1:,:3]

positions_df.columns = ["CHROM", "START", "END"]


positions_df["POS"] = [ list(range(x, y+1)) for x, y in positions_df[["START", "END"]].values ]
positions_df = positions_df.explode("POS").reset_index(drop = True)
positions_df = positions_df[["CHROM", "POS"]]
positions_df.columns = ["CHROM", "POS"]


positions_df["REF"] = positions_df.apply(get_pos_in_row, axis = 1)
positions_df["ALT"] = positions_df["REF"].apply(get_non_ref)

positions_df = positions_df.explode("ALT").reset_index(drop = True)
positions_df.index.name = "SAMPLE"


positions_df = positions_df.reset_index(drop = True)


# 1   881907    881906    -/C   +
positions_df["MUTATION"] = positions_df["REF"].astype(str) + "/" + positions_df["ALT"].astype(str)
positions_df["STRAND"] = "+"
positions_df_end = positions_df[['CHROM', 'POS', 'POS', 'MUTATION', 'STRAND']]


positions_df_end.to_csv(output_file_with_sites,
						header = False,
						index = False,
						sep = "\t")
