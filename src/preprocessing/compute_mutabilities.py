import itertools
import pandas as pd
import numpy as np

from utils import *

def read_inputs(all_possible_sites_annotated_file, depth_dataframe_file, mutations_file):
    """
    This function reads the three files needed for running the preprocessing
    """

    # Read all possible mutations annotated by VEP
    all_possible_sites_annotated = pd.read_csv(all_possible_sites_annotated_file, sep = "\t", header = 0)

    # Read depth matrix
    depth_dataframe = pd.read_csv(depth_dataframe_file, header = 0, sep = "\t")
    depth_dataframe.columns = ["CHROM", "POS"] + list(depth_dataframe.columns[2:])
    depth_dataframe["CHROM"] = depth_dataframe["CHROM"].astype(str)
    depth_dataframe["POS"] = depth_dataframe["POS"].astype(int)

    # Read MAF with the mutations from all the samples
    # it needs to have at least these columns:
    #   'CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID'
    maf = pd.read_csv(mutations_file, sep = "\t", header = 0, dtype= {"CHROM" : str, "POS": int,
                                                                        "REF" : str, "ALT" : str,
                                                                        "SAMPLE_ID" : str}
                                                                        )

    return all_possible_sites_annotated, depth_dataframe, maf


def compute_mutations_per_sample_gene_impact_context_table(maf, all_possible_sites_annotated):
    """
    This function receives:
        - a mutations dataframe
        - a dataframe with all the possible SNV sites to be mutated
    and returns:
        - The observed mutations annotated
        - a table with the number of mutations per sample, gene, impact and context
    """

    minimal_maf = maf[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID']].copy()

    # select only SNVs
    minimal_maf["TYPE"] = minimal_maf[['REF', 'ALT']].apply(vartype, axis = 1)
    minimal_maf = minimal_maf[minimal_maf["TYPE"] == "SNV"].reset_index(drop = True)
    minimal_maf = minimal_maf.drop("TYPE", axis = 1)

    ##
    # Annotate observed mutations
    ##
    annotated_minimal_maf = minimal_maf.merge(all_possible_sites_annotated, on = ["CHROM", "POS", "REF", "ALT"], how = "left")

    ##
    # Count of all observed mutations per sample, gene, impact and context
    ##
    #   Required information:
    #       Annotated mutations observed
    #   Output:
    #       observed mutations per:
    #           sample
    #           gene
    #           impact
    #           context
    ##
    obs_muts_per_gene_impact_context_sample_long = annotated_minimal_maf.groupby(by = ['SAMPLE_ID', "GENE", "IMPACT", "CONTEXT_MUT"])["MUT_ID"].count()
    obs_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long.reset_index()
    obs_muts_per_gene_impact_context_sample_long.columns = ['SAMPLE_ID', "GENE", "IMPACT", "CONTEXT_MUT", "COUNT"]

    # wide format
    obs_muts_per_gene_impact_context_sample_wide = obs_muts_per_gene_impact_context_sample_long.pivot(
                                                                index= ['GENE', "IMPACT", "CONTEXT_MUT"], columns='SAMPLE_ID', values='COUNT').fillna(0).astype(int).reset_index()
    obs_muts_per_gene_impact_context_sample_wide.columns.name = None
    
    return annotated_minimal_maf, obs_muts_per_gene_impact_context_sample_wide






def define_samples(depth_dataframe, annotated_minimal_maf):
    ##
    # Select all samples with available data in depth and in mutations
    ##
    samples_muts = list(annotated_minimal_maf["SAMPLE_ID"].unique())
    samples_depths = list(depth_dataframe.columns[2:])
    samples = sorted(list(set(samples_muts).intersection(samples_depths)))
    
    print(f"{len(samples)} samples maintained starting from {len(samples_muts)} samples with mutations info and {len(samples_depths)} samples with depths info.")
    
    if len(set(samples_muts) - set(samples)) > 0:
        print(f"Removed samples with mutations info: {sorted(set(samples_muts) - set(samples))}")
    if len(set(samples_depths) - set(samples)) > 0:
        print(f"Removed samples with depths info: {sorted(set(samples_depths) - set(samples))}")
    
    return samples





def compute_mutational_profile(annotated_minimal_maf, all_possible_sites_annotated, depth_dataframe, samples, pseudocount = 0):
    """
    Compute mutational profile from the input data
          ***Remember to add some pseudocounts to the computation***
    
        Required information:
            Annotated all possible sites
            Annotated mutations observed
            Depth matrix
        Output:
            Mutational profile per sample, possibility to add pseudocounts to prevent some probabilities from being 0
    """
    subs = [''.join(z) for z in itertools.product('CT', 'ACGT') if z[0] != z[1]]
    flanks = [''.join(z) for z in itertools.product('ACGT', repeat=2)]
    contexts_unformatted = sorted([(a, b) for a, b in itertools.product(subs, flanks)], key=lambda x: (x[0], x[1]))
    contexts_no_change = [b[0]+a[0]+b[1] for a, b in contexts_unformatted]
    contexts_formatted = [b[0]+a[0]+b[1]+'>'+a[1] for a, b in contexts_unformatted]

    # create the matrix in the desired order
    empty_matrix = pd.DataFrame(index = contexts_formatted)

    # count the mutations per sample and per context
    counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["MUT_ID"].count().reset_index()
    counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "MUT_ID")
    counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
    counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
    counts_x_sample_matrix.index.name = "CONTEXT_MUT"

    # add a pseudocount if desired
    counts_x_sample_matrix = counts_x_sample_matrix + pseudocount
    # here we have the counts matrix for the number of mutations per sample per context


    # Compute the trinucleotide depth per sample
    # merge the dataframe of all possible sites with the dataframe of the depth per site per sample
    trinuc_depth_per_sample = all_possible_sites_annotated.merge(depth_dataframe,
                                                                    on = ["CHROM", "POS"],
                                                                    how = "left").groupby(by = "CONTEXT_MUT")[samples].sum()

    # divide
    mut_probability = counts_x_sample_matrix.divide( trinuc_depth_per_sample )
    # print(mut_probability.head())

    # normalize
    mut_probability = mut_probability / mut_probability.sum()
    # print(mut_probability.head())

    # reindex to ensure the right order
    mut_probability = mut_probability.reindex(contexts_formatted)
    mut_probability.index.name = "CONTEXT_MUT"
    mut_probability = mut_probability.reset_index()
    # mut_probability
    # print(mut_probability.head())

    return mut_probability




def compute_expected_synonymous_mutations(all_possible_sites_annotated, depth_dataframe, mut_probability, samples):
    """
    Required information:
            All possible sites in the panel regions
            Depth per site per sample (only used as a binary dataframe, covered or not covered)
            mutational probabilities per sample
    Output:
        sites per:
            sample
            gene
            impact
            context
    1. Count of all possible sites per sample (keeping only sites with enough depth)
    2. Use the mutational probability information to estimate the number of observed synonymous mutations using relative values
    """
    binary_depth_dataframe = (depth_dataframe.set_index(["CHROM", "POS"]) > 0).astype(int).reset_index()
    all_possible_sites_per_sample = all_possible_sites_annotated.merge(binary_depth_dataframe, on = ["CHROM", "POS"], how = "left")

    # wide format
    sites_per_gene_impact_context_sample_wide = all_possible_sites_per_sample.groupby(
                                                                by = ["GENE", "IMPACT", "CONTEXT_MUT"])[samples].sum().reset_index()

    # # long format
    # sites_per_gene_impact_context_sample_long = sites_per_gene_impact_context_sample_wide.melt(id_vars = ["GENE", "IMPACT", "CONTEXT_MUT"],
    #                                                                                             var_name = "SAMPLE_ID",
    #                                                                                             value_name = "COUNT")
    # # sites_per_gene_impact_context_sample_long


    ## Get synonymous sites counts per gene, impact, context and sample
    syn_sites_per_gene_impact_context_sample = sites_per_gene_impact_context_sample_wide[
                                                        sites_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"
                                                    ].reset_index(drop = True)

    syn_sites_per_gene_impact_context_sample[samples] = syn_sites_per_gene_impact_context_sample[samples].fillna(0).astype(int)

    # the number of synonymous sites per gene should be the same for all samples
    # except if there are differences in sequencing coverage of those areas

    # Merge synonymous sites with mutational probabilities
    synonymous_sites_mut_probs = syn_sites_per_gene_impact_context_sample.merge(mut_probability,
                                                                                suffixes = [".sites", ".probability"],
                                                                                on = "CONTEXT_MUT")
    # synonymous_sites_mut_probs.head()


    synonymous_probs_gene_context = synonymous_sites_mut_probs[["GENE", "IMPACT", "CONTEXT_MUT"]].copy()
    for sample in samples:
        synonymous_probs_gene_context[sample] = synonymous_sites_mut_probs[f"{sample}.sites"] * synonymous_sites_mut_probs[f"{sample}.probability"]
    # synonymous_probs_gene_context



    # Here we add the mutation probability of all synonymous sites in each gene and sample
    #       these are the expected synonymous mutations taking into account only the mutational profile
    ###
    # these are the values that should be compared to the observed mutations,
    #    to then adjust the mutability
    expected_syn_per_gene_per_sample = synonymous_probs_gene_context.drop(["IMPACT", "CONTEXT_MUT"],
                                                                            axis = 1).groupby("GENE").sum().reset_index()
    # print(expected_syn_per_gene_per_sample.head())
    
    return expected_syn_per_gene_per_sample


def compute_mutabilities(alpha_per_sample, mut_probability, samples):

    mutability_all_samples = pd.DataFrame()

    #####
    # Compute the mutabilities
    #####
    for sample in samples:
        mutability_sample = pd.DataFrame()

        # iterate over all genes using the name of the gene and
        #    the alpha for which we should correct the mutability
        for gen, alpha in alpha_per_sample[sample].items():
            # print(gen, alpha)

            # take the mutation probability computed from the mutations observed
            # and the sequencing depth of each context
            mut_probability_sample_gene = mut_probability[["CONTEXT_MUT", sample]].copy()

            # adjust the probability vector by the value of alpha
            #   corresponding to that particular gene in that sample
            mut_probability_sample_gene[sample] = mut_probability_sample_gene[sample] * alpha
            mut_probability_sample_gene["GENE"] = gen

            mutability_sample = pd.concat( (mutability_sample, mut_probability_sample_gene), axis = 0)
            # mut_probability_sample_gene.to_csv(f"{mutability_path}/mutability.{sample}.{gen}.tsv",
            #                                     header = True,
            #                                     index = False,
            #                                     sep = "\t")
        
        mutability_sample = mutability_sample.set_index(["GENE", "CONTEXT_MUT"])
        mutability_all_samples = pd.concat( (mutability_all_samples, mutability_sample), axis = 1)

    return mutability_all_samples.reset_index()



def compute_mutabilities_wrapper(all_possible_sites_annotated_file,
                                    depth_dataframe_file,
                                    mutations_file, 
                                    table_muts_x_sample_gene_impact_context,
                                    mutability_path
                                    ):
    """
    Wrapper for all the steps required to compute the mutabilities per sample, gene and context
    """
    # Read files
    all_possible_sites_annotated, depth_dataframe, maf = read_inputs(all_possible_sites_annotated_file,
                                                                        depth_dataframe_file,
                                                                        mutations_file
                                                                        )

    # Annotate mutations and compute table of observed mutations
    annotated_minimal_maf, obs_muts_per_gene_impact_context_sample_wide = compute_mutations_per_sample_gene_impact_context_table(
                                                                                maf, all_possible_sites_annotated
                                                                                )

    # Define for which samples we have enough data
    samples = define_samples(depth_dataframe, annotated_minimal_maf)
    
    # focus the mutations and depth dataframes into the selected samples
    annotated_minimal_maf = annotated_minimal_maf[annotated_minimal_maf["SAMPLE_ID"].isin(samples)].copy().reset_index(drop = True)
    depth_dataframe = depth_dataframe[["CHROM", "POS"] + samples].copy()


    # select only the samples that will be part of the analysis and store the table
    # obs_muts_per_gene_impact_context_sample_long = long
    obs_muts_per_gene_impact_context_sample_wide = obs_muts_per_gene_impact_context_sample_wide[
                                                                    ["GENE", "IMPACT", "CONTEXT_MUT"] + samples
                                                                ].copy()


    # *** OUTPUT 1 ***
    ## Store table with observed mutations 
    obs_muts_per_gene_impact_context_sample_wide.to_csv(f"{table_muts_x_sample_gene_impact_context}",
                                                        header = True,
                                                        index = False,
                                                        sep = "\t")


    # Compute mutational profile from the input data
    mut_probability = compute_mutational_profile(annotated_minimal_maf,
                                                    all_possible_sites_annotated,
                                                    depth_dataframe,
                                                    samples,
                                                    pseudocount = 0.5)

    expected_syn_per_gene_per_sample = compute_expected_synonymous_mutations(all_possible_sites_annotated,
                                                                                depth_dataframe,
                                                                                mut_probability,
                                                                                samples)

    # Count of observed synonymous mutations per sample and gene
    obs_syn_muts_per_gene_context_sample = obs_muts_per_gene_impact_context_sample_wide[
                                                    obs_muts_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"].reset_index(
                                                        drop = True)
    obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_context_sample.groupby(by = ["GENE"])[samples].sum()
    obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_sample.reset_index()

    #####
    # Alpha computation
    #####
    # What is alpha?
    # Alpha is a factor that allows us to update mutability per context
    # from relative to absolute in each specific gene-sample pair
    #####

    # we compute the value of alpha per each gene-sample pair,
    # by dividing the number of expected synonymous
    # by the number of synonymous we would be generating with the original mutational profile
    alpha_per_sample = obs_syn_muts_per_gene_sample.set_index("GENE").divide( expected_syn_per_gene_per_sample.set_index("GENE") )
    # print(alpha_per_sample)

    # compute the mutabilities by adjusting the mutational profile (mut_probability)
    # by the value of alpha, to obtain an absolute mutability per context
    mutability_all_samples = compute_mutabilities(alpha_per_sample, mut_probability, samples)


    # *** OUTPUT 2 *** 
    # create a single table with all the mutabilities per sample, gene, context
    mutability_all_samples.to_csv(f"{mutability_path}",
                                        header = True,
                                        index = False,
                                        sep = "\t")






if __name__ == '__main__':

    ## Input
    depth_dataframe_file = "/workspace/datasets/prominent/data/kidney/depth/2023-06-30.kidney_panel.chr.633.tsv.gz"
    mutations_file = "/workspace/datasets/prominent/data/kidney/mutations/2023-06-30.kidney.633.maf.annot.tsv.gz"
    all_possible_sites_annotated_file = "./test/preprocessing/KidneyPanel.sites.bed_panel.annotation_summary.tsv"

    ## Output
    table_muts_x_sample_gene_impact_context = "./test/preprocessing/mutations_per_sample_gene_impact_context.count.tsv"
    mutabilities_table = "./test/preprocessing/mutability_per_sample_gene_context.tsv"

    compute_mutabilities_wrapper(all_possible_sites_annotated_file,
                                    depth_dataframe_file,
                                    mutations_file, 
                                    table_muts_x_sample_gene_impact_context,
                                    mutabilities_table
                                    )


