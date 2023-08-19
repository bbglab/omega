depths_file = ""
mutations_file = ""

run_vep = False
postprocess_vep = False
compute_mutabilities = False

bed_regions_file = ""
vep_output_file = ""
vep_postprocessed_file = ""

mutabilities_dir = ""

python src/preprocessing/sites_table_from_bed.py
# # input_bedfile = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
# input_bedfile = sys.argv[1]
# # output_file_with_sites = ./test/preprocessing/KidneyPanel.sites4VEP.tsv"
# output_file_with_sites = sys.argv[2]


## run Ensembl VEP from their web server


python src/preprocessing/postprocessing_annotation.py 
# # Input
# VEP_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
# # Output
# all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"




python src/preprocessing/compute_mutabilities.py 
# ## Input
# depth_dataframe_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"

# mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"

# all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

# ## Output
# table_muts_x_sample_gene_impact_context = "./test/preprocessing/mutations_per_gene_impact_context.count.tsv"
# mutability_path = "./test/preprocessing/mutabilities"
# # check that mutability_path exists or otherwise create it