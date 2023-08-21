depths_file = ""                # depths_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"
mutations_file = ""             # mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"

run_vep = False
postprocess_vep = False
compute_mutabilities = False

bed_regions_file = ""           # input_bedfile = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
vep_input_file = ""             # output_file_with_sites = ./test/preprocessing/KidneyPanel.sites4VEP.tsv"
vep_output_file = ""            # vep_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
vep_postprocessed_file = ""     # vep_postprocessed_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

# Outputs
table_observed_muts = ""        # table_observed_muts = "./test/preprocessing/mutations_per_gene_impact_context.count.tsv"
mutabilities_dir = ""           # mutabilities_dir = "./test/preprocessing/mutabilities"
# TODO: check that mutabilities_dir exists or otherwise create it,
##         we could also check if the mutabilities have already been computed previously



from src.preprocessing.sites_table_from_bed import generate_all_sites_4VEP
from src.preprocessing.postprocessing_annotation import vep2summarizedannotation
from src.preprocessing.compute_mutabilities import compute_mutabilities_wrapper



if run_vep:
    generate_all_sites_4VEP(bed_regions_file, vep_input_file)
# else check if vep_input_file, or vep_output_file or vep_postprocessed_file exists

## run Ensembl VEP from their web server
# TODO add action point here where the user needs to do things
### this would imply using "vep_input_file" to generate "vep_output_file"
### the user should introduce the path to vep_output_file

if postprocess_vep:
    vep2summarizedannotation(vep_output_file, vep_postprocessed_file)
# else check if vep_postprocessed_file exists

if compute_mutabilities:
    compute_mutabilities_wrapper(vep_postprocessed_file, depths_file, mutations_file, table_observed_muts, mutabilities_dir)