depths_file = ""
mutations_file = ""

run_vep = False
postprocess_vep = False
compute_mutabilities = False

bed_regions_file = ""
vep_output_file = ""
vep_postprocessed_file = ""

mutabilities_dir = ""

python sites_table_from_bed.py
# input_bedfile = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
input_bedfile = sys.argv[1]
# output_file_with_sites = f"{positive_sel_dir}/KidneyPanel.all_SNVs.bed.panel.4VEP.tsv"
output_file_with_sites = sys.argv[2]


## run Ensembl VEP from their web server


python postprocessing_annotation.py 
# Input
# VEP_output_file = "/home/fcalvet/projects/omega/omega/tests_ferriol/KidneyGenes.canonical_transcripts_CDS.VEPannotated.tsv"
VEP_output_file = f"/home/fcalvet/projects/omega/omega/tests_ferriol/KidneyPanel.sites.VEP_annotated.tsv"

# Output
all_possible_sites_annotated_file = "/workspace/datasets/transfer/ferran_to_ferriol/omega_tests/KidneyPanel.all_SNVs.bed_panel.annotation_summary2.tsv"




python compute_mutabilities.py 

depth_dataframe_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"

mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"

all_possible_sites_annotated_file = "/workspace/datasets/transfer/ferran_to_ferriol/omega_tests/KidneyPanel.all_SNVs.bed_panel.annotation_summary.tsv"



## Output
# 
table_muts_x_sample_gene_impact_context = "/workspace/datasets/transfer/ferran_to_ferriol/v2023-08-02_data/mutations_per_gene_impact_context.count.tsv"
mutability_path = "/home/fcalvet/projects/omega/omega/tests_ferriol/mutabilities2"
# check that mutability_path exists or otherwise create it
