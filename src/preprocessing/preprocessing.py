# depths_file = ""                # depths_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"
# mutations_file = ""             # mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"

# run_vep = True
# postprocess_vep = True
# compute_mutabilities = True

# bed_regions_file = ""           # input_bedfile = "/workspace/datasets/prominent/metadata/regions/data/oncodrivefml/kidneypanel4oncodrivefml.bed5.bed"
# vep_input_file = ""             # output_file_with_sites = ./test/preprocessing/KidneyPanel.sites4VEP.tsv"
# vep_output_file = ""            # vep_output_file = f"./test/preprocessing/KidneyPanel.sites.VEP_annotated.tsv"
# vep_postprocessed_file = ""     # vep_postprocessed_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

# # Outputs
# table_observed_muts = ""        # table_observed_muts = "./test/preprocessing/mutations_per_sample_gene_impact_context.count.tsv"
# mutabilities_table = ""         # mutabilities_table = "./test/preprocessing/mutability_per_sample_gene_context.tsv"
# # TODO: we could check if the mutabilities have already been computed previously



depths_file = "./test/input/full_case/2023-06-30.kidney_panel.chr.633.tsv.gz"
mutations_file = "./test/input/full_case/2023-06-30.kidney.633.maf.annot.tsv.gz"
bed_regions_file = "./test/input/full_case/kidneypanel4oncodrivefml.bed5.bed"

vep_input_file = "./test/input/full_case/KidneyPanel.sites4VEP.tsv"
vep_postprocessed_file = "./test/input/full_case/KidneyPanel.sites.bed_panel.annotation_summary.tsv"
vep_output_file = ""

run_vep = True
postprocess_vep = True
compute_mutabilities = True


# Outputs
table_observed_muts = "./test/input/full_case/mutations_per_sample_gene_impact_context.count.tsv"
mutabilities_table = "./test/input/full_case/mutability_per_sample_gene_context.tsv"


from src.preprocessing.sites_table_from_bed import generate_all_sites_4VEP
from src.preprocessing.postprocessing_annotation import vep2summarizedannotation
from src.preprocessing.compute_mutabilities import compute_mutabilities_wrapper


# TODO we could check for the presence of the chr prefix in the files we read just to make sure this is not causing problems


if run_vep or vep_output_file == "":
    generate_all_sites_4VEP(bed_regions_file, vep_input_file)
# else check if vep_input_file, or vep_output_file or vep_postprocessed_file exists

## run Ensembl VEP from their web server
# TODO add action point here where the user needs to do things
### this would imply using "vep_input_file" to generate "vep_output_file"
### the user should introduce the path to vep_output_file

    print(f"Use this file to run Ensembl VEP:\n{vep_input_file}")
    print( "You can run it here: https://www.ensembl.org/Homo_sapiens/Tools/VEP, and then download the results." )
    vep_output_file = input("Introduce the path to the EnsemblVEP annotated file:\n")
    # TODO check if the file exists
    vep_output_file = vep_output_file.strip()

if postprocess_vep:
    vep2summarizedannotation(vep_output_file, vep_postprocessed_file)
# else check if vep_postprocessed_file exists

if compute_mutabilities:
    compute_mutabilities_wrapper(vep_postprocessed_file, depths_file, mutations_file, table_observed_muts, mutabilities_table)