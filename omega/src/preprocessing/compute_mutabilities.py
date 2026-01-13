import os

import daiquiri
import pandas as pd

from omega import __logger_name__
from omega.src.utils import canonical_channels
from omega.src.preprocessing.utils import vartype

logger = daiquiri.getLogger(__logger_name__ + '.preprocessing.comp_mutabs')

CHANNELS = canonical_channels()

def read_inputs(all_possible_sites_annotated_file, depth_dataframe_file, mutations_file):
    """
    This function reads the three files needed for running the preprocessing
    """
    # logger.debug("{} {} {}".format(all_possible_sites_annotated_file, depth_dataframe_file, mutations_file))
    # Read all possible mutations annotated by VEP
    all_possible_sites_annotated = pd.read_csv(all_possible_sites_annotated_file,
                                                sep = "\t", header = 0,
                                                dtype = {"CHROM" : str, "POS": int,
                                                            "REF" : str, "ALT" : str,
                                                            "MUT_ID" : str, "GENE" : str,
                                                            "IMPACT" : str, "CONTEXT_MUT" : str} )
    logger.debug("All sites loaded")

    # Read depth matrix
    depth_dataframe = pd.read_csv(depth_dataframe_file, header = 0, sep = "\t")
    depth_dataframe.columns = ["CHROM", "POS"] + list(depth_dataframe.columns[2:])
    if "CONTEXT" in depth_dataframe.columns: depth_dataframe = depth_dataframe.drop("CONTEXT", axis = 1)
    depth_dataframe["CHROM"] = depth_dataframe["CHROM"].astype(str)
    depth_dataframe["POS"] = depth_dataframe["POS"].astype(int)
    logger.debug("Depths loaded")

    # Read MAF with the mutations from all the samples
    # it needs to have at least these columns:
    #   'CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID'
    maf = pd.read_csv(mutations_file, sep = "\t", header = 0, dtype= {"CHROM" : str, "POS": int,
                                                                        "REF" : str, "ALT" : str,
                                                                        "SAMPLE_ID" : str}
                                                                        )
    logger.debug("MAF loaded")

    if "EFFECTIVE_MUTS" in maf:
        maf["EFFECTIVE_MUTS"] = maf["EFFECTIVE_MUTS"].astype(float)
        logger.debug("Using EFFECTIVE_MUTS provided.")
    else:
        maf["EFFECTIVE_MUTS"] = 1.
        logger.debug("Counting each mutation only once")

    # FIXME revise this part of counting mutations only once,
    # maybe since we keep the sample and count once per sample this is a smaller problem


    # make sure that all the files have the chr prefix in the files
    if not all_possible_sites_annotated["CHROM"].iloc[0].startswith("chr"):
        all_possible_sites_annotated["CHROM"] = "chr" + all_possible_sites_annotated["CHROM"]
    
    if not depth_dataframe["CHROM"].iloc[0].startswith("chr"):
        depth_dataframe["CHROM"] = "chr" + depth_dataframe["CHROM"]

    if not maf["CHROM"].iloc[0].startswith("chr"):
        maf["CHROM"] = "chr" + maf["CHROM"]

    return all_possible_sites_annotated, depth_dataframe, maf


def annotate_mutations_using_vep(maf, all_possible_sites_annotated):
    """
    This function receives:
        - a mutations dataframe
        - a dataframe with all the possible SNV sites to be mutated
    and returns:
        - The observed mutations annotated
    """

    minimal_maf = maf[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID', 'EFFECTIVE_MUTS']].copy()

    # select only SNVs
    minimal_maf["TYPE"] = minimal_maf[['REF', 'ALT']].apply(vartype, axis = 1)
    minimal_maf = minimal_maf[minimal_maf["TYPE"] == "SNV"].reset_index(drop = True)
    minimal_maf = minimal_maf.drop("TYPE", axis = 1)

    ###
    ## Annotate observed mutations
    # I selected the inner strategy so that variants that fall outside of the regions
    # of possible sites are not included in any of the analysis nor the counts
    ###
    # TODO revise whether we want to keep it as inner or we want left
    # I prefer inner since if the user provides mutations outside
    # the areas of interest those are excluded from the analysis
    annotated_minimal_maf = minimal_maf.merge(all_possible_sites_annotated, on = ["CHROM", "POS", "REF", "ALT"], how = "inner")

    return annotated_minimal_maf


def compute_mutations_per_sample_gene_impact_context_table(annotated_minimal_maf,
                                                            impacts_to_exclude = ["non_genic_variant", "intron_variant"]
                                                            ):
    # TODO: revise the default list of excluded impacts
    """
    This function receives:
        - The observed mutations annotated
    and returns:
        - a table with the number of mutations per sample, gene, impact and context
    
    here the total number of mutations per sample, gene, impact and context can be higher than the number of sites when counting them 1x
    if the effective_muts are being computed with the ALT_DEPTH
    """
    # TODO: revise whether we are interested in doing it this way or not filter annotations/positions that we are not interested in
    annotated_minimal_maf = annotated_minimal_maf[annotated_minimal_maf["GENE"] != '-']
    annotated_minimal_maf = annotated_minimal_maf[~annotated_minimal_maf["IMPACT"].isin(impacts_to_exclude)].reset_index(drop = True)

    obs_muts_per_gene_impact_context_sample_long = annotated_minimal_maf.groupby(by = ['SAMPLE_ID', "GENE", "IMPACT", "CONTEXT_MUT"])["EFFECTIVE_MUTS"].sum()
    obs_muts_per_gene_impact_context_sample_long = obs_muts_per_gene_impact_context_sample_long.reset_index()
    obs_muts_per_gene_impact_context_sample_long.columns = ['SAMPLE_ID', "GENE", "IMPACT", "CONTEXT_MUT", "COUNT"]
    
    # wide format
    obs_muts_per_gene_impact_context_sample_wide = obs_muts_per_gene_impact_context_sample_long.pivot(
                                                                index= ['GENE', "IMPACT", "CONTEXT_MUT"], columns='SAMPLE_ID', values='COUNT').fillna(0).astype(float).reset_index()
    obs_muts_per_gene_impact_context_sample_wide.columns.name = None

    return obs_muts_per_gene_impact_context_sample_wide






def define_samples(depth_dataframe, annotated_minimal_maf):
    ##
    # Select all samples with available data in depth and in mutations
    ##
    samples_muts = list(annotated_minimal_maf["SAMPLE_ID"].unique())
    samples_depths = list(depth_dataframe.columns[2:])
    samples = sorted(list(set(samples_muts).intersection(samples_depths)))

    logger.info(f"{len(samples)} samples maintained starting from {len(samples_muts)} samples with mutations info and {len(samples_depths)} samples with depths info.")
    
    if len(set(samples_muts) - set(samples)) > 0:
        logger.info(f"Removed samples with mutations info: {sorted(set(samples_muts) - set(samples))}")
    if len(set(samples_depths) - set(samples)) > 0:
        logger.info(f"Removed samples with depths info: {sorted(set(samples_depths) - set(samples))}")
    
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

    # create the matrix in the desired order
    empty_matrix = pd.DataFrame(index = CHANNELS)

    # make sure to count each mutation only once (avoid annotation issues)
    annotated_minimal_maf = annotated_minimal_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "EFFECTIVE_MUTS"]].drop_duplicates().reset_index(drop = True)
    # TODO: maybe we should make sure that no mutation is counted
    # if it falls in a position outside the ones for which we have values of depth?
    # it should not happen but who knows...
    ###
    # it could work by defining the intersection of the three dataframes
    # in terms of positions the same way as we do with the samples



    # count the mutations per sample and per context
    counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["EFFECTIVE_MUTS"].sum().reset_index()
    counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "EFFECTIVE_MUTS")
    counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
    counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
    counts_x_sample_matrix.index.name = "CONTEXT_MUT"

    # add a pseudocount if desired
    counts_x_sample_matrix = counts_x_sample_matrix + pseudocount
    # here we have the counts matrix for the number of mutations per sample per context


    # Compute the trinucleotide depth per sample
    # merge the dataframe of all possible sites with the dataframe of the depth per site per sample
    
    ## TODO make a decision here (use depth or counts?)
    # I think depth makes more sense
    # in case there are regions with very small coverage we can make sure
    # that they don't have the same contribution to the background counts...
    ###
    # We could also normalize by the trinucleotide counts, not trinucleotide depth
    # # how = "left").groupby(by = "CONTEXT_MUT")[samples].count()
    # # revise how it handles the 0s in coverage, are they counted? if so, we should use the binary_depth_dataframe
    ###

    # make sure to count each site only once (avoid annotation issues)
    all_possible_sites_annotated = all_possible_sites_annotated[["CHROM", "POS", "CONTEXT_MUT"]].drop_duplicates().reset_index(drop = True)
    trinuc_depth_per_sample = all_possible_sites_annotated.merge(depth_dataframe,
                                                                    on = ["CHROM", "POS"],
                                                                    how = "left").groupby(by = "CONTEXT_MUT")[samples].sum()
                                                                    # how = "left").groupby(by = "CONTEXT_MUT")[samples].count()

    # divide
    mut_probability = counts_x_sample_matrix.divide( trinuc_depth_per_sample )
    # logger.debug(mut_probability.head())

    # normalize
    mut_probability = mut_probability / mut_probability.sum()
    # logger.debug(mut_probability.head())

    # reindex to ensure the right order
    mut_probability = mut_probability.reindex(CHANNELS)
    mut_probability.index.name = "CONTEXT_MUT"
    mut_probability = mut_probability.reset_index()
    # mut_probability
    # logger.debug(mut_probability.head())

    return mut_probability




def compute_expected_synonymous_mutations(all_possible_sites_annotated, depth_dataframe, mut_probability, samples, single_sample = False):
    """
    Goal:
        Based on the number of sites, depth, and mutational profile, infer how many synonymous mutation we would expect to see.
        This comes from absolute vales, such as the depth, but also relative values such as the mutation probablity,
        however we are not using the absolute value of depth at any time, only as a measure of a position being covered or not.

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
    if single_sample: # this is wrong, here if a group of samples has 25 samples each site should be counted 25 times
        logger.debug("Running in single sample mode.")
        logger.debug(all_possible_sites_per_sample.columns)
        samples_names = [x for x in all_possible_sites_per_sample.columns if x not in ["CHROM", "POS", "REF", "ALT", "GENE", "IMPACT", "CONTEXT_MUT"] ]
        logger.debug(samples_names)

        # this number should be at most, the total number of samples in the group that is being run
        # whenever you are running a group this should never be 1...
        sites_per_gene_impact_context_sample_wide = all_possible_sites_per_sample.groupby(
                                                                    by = ["GENE", "IMPACT", "CONTEXT_MUT"])[samples_names].sum().sum(axis=1).reset_index()
        sites_per_gene_impact_context_sample_wide.columns = ["GENE", "IMPACT", "CONTEXT_MUT", single_sample]

    else:
        sites_per_gene_impact_context_sample_wide = all_possible_sites_per_sample.groupby(
                                                                    by = ["GENE", "IMPACT", "CONTEXT_MUT"])[samples].sum().reset_index()


    ## Get synonymous sites counts per gene, impact, context and sample
    syn_sites_per_gene_impact_context_sample = sites_per_gene_impact_context_sample_wide[
                                                        sites_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"
                                                    ].reset_index(drop = True)

    syn_sites_per_gene_impact_context_sample[samples] = syn_sites_per_gene_impact_context_sample[samples].fillna(0).astype(float)

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
    # logger.debug(expected_syn_per_gene_per_sample.head())
    
    return expected_syn_per_gene_per_sample

def compute_sample_gene_specific_differences(all_possible_sites_annotated, depth_dataframe,
                                                mut_probability_total,
                                                mut_probability,
                                                samples, single_sample = False):
    """
    Here we compute this: 
        (depth * (sample_mut_profile/all_samples_mut_profile) )

    Required information:
            All possible sites in the panel regions
            Depth per site per sample
            Mutational profile of the cohort in which the global mutation rates were computed
            Mutational profile of the sample or samples used here
    Output:
        context corrected depth per trinucleotide per gene
            (corrected by differences in the mutational profile of the sample and the cohort)
    """
    all_possible_synonymous_sites = all_possible_sites_annotated[all_possible_sites_annotated["IMPACT"] == "synonymous"][["CHROM", "POS", "GENE", "CONTEXT_MUT"]].reset_index(drop = True)
    
    # annotate the synonymous sites with depth information
    all_synonymous_sites_depth = all_possible_synonymous_sites.merge(depth_dataframe, on = ["CHROM", "POS"], how = "left")
    # contains information about: gene, chromosome, positions, trinucleotide and depth
    # logger.debug("Depth in synonymous sites of the gene")
    # print(all_synonymous_sites_depth.head())
    # print(all_synonymous_sites_depth.sum())



    # # wide format
    # if single_sample:
    #     logger.debug("Running in single sample mode.")
    #     logger.debug(all_synonymous_sites_depth.columns)
    #     samples_names = [x for x in all_synonymous_sites_depth.columns if x not in ["CHROM", "POS", "REF", "ALT", "GENE", "IMPACT", "CONTEXT_MUT"] ]
    #     logger.debug(samples_names)
    #     depth_per_context_gene_sample_wide = all_synonymous_sites_depth.groupby(
    #                                                                 by = ["GENE", "CONTEXT_MUT"])[samples_names].sum().sum(axis=1).reset_index()
    #     depth_per_context_gene_sample_wide.columns = ["GENE", single_sample]
    
    samples_names = [x for x in all_synonymous_sites_depth.columns if x not in ["CHROM", "POS", "REF", "ALT", "GENE", "IMPACT", "CONTEXT_MUT"] ]
    logger.debug(samples_names)

    # group all positions by GENE n CONTEXT_MUT,
    # we do not expect any differences within trinucleotides of the same gene, except for the depth,
    # but we are not looking at specific position differences here, but more at the gene level
    depth_per_context_gene_sample_wide = all_synonymous_sites_depth.groupby(
                                                                by = ["GENE", "CONTEXT_MUT"])[samples].sum().reset_index()
    depth_per_context_gene_sample_wide[samples] = depth_per_context_gene_sample_wide[samples].fillna(0).astype(float)
    
    print(depth_per_context_gene_sample_wide.head())

    mut_probability_ind = mut_probability.set_index("CONTEXT_MUT")
    mut_probability_total_ind = mut_probability_total.set_index("CONTEXT_MUT")

    print(mut_probability_ind.sum())
    print(mut_probability_total_ind.sum())

    # TODO
    # revise if this division makes sense or we would need to center the resulting vector or something

    ## URGENT

    # normalize sample specific mutational profile compared to the all_samples one
    mut_probability_norm = (mut_probability_ind / mut_probability_total_ind).reset_index()
    # mut_probability_norm = ( mut_probability_total_ind / mut_probability_ind).reset_index()
    # print(mut_probability_total_ind / mut_probability_ind)
    # print((mut_probability_total_ind / mut_probability_ind).sum())
    # print(mut_probability_ind / mut_probability_total_ind)
    # print((mut_probability_ind / mut_probability_total_ind).sum())
    # print(mut_probability_norm)

    # merge the depth per context gene with the correction of the mutation probability
    probability_depth_per_context_gene_sample_wide = depth_per_context_gene_sample_wide.merge(mut_probability_norm,
                                                                                  on = 'CONTEXT_MUT',
                                                                                  suffixes = [".sites", ".probability_correction"],
                                                                                  )
    print(probability_depth_per_context_gene_sample_wide.head())

    # apply the correction of difference in mutational profile
    corrected_depth_per_context_gene_sample_wide = probability_depth_per_context_gene_sample_wide[["GENE", "CONTEXT_MUT"]].copy()
    for sample in samples:
        corrected_depth_per_context_gene_sample_wide[sample] = probability_depth_per_context_gene_sample_wide[f"{sample}.sites"] * probability_depth_per_context_gene_sample_wide[f"{sample}.probability_correction"]

    # add up all the corrections to gene level depth
    context_corrected_depth_per_gene_sample_wide = corrected_depth_per_context_gene_sample_wide.groupby(
                                                                    by = ["GENE"])[samples].sum()

    # print(context_corrected_depth_per_gene_sample_wide)
    # then we should multiply this by the gene mutation rate per bp
    
    # gene_mutrate_per_bp_allsamples
    # should be a GENE indexed pandas.Series that contains what it's name suggests

    # result = context_corrected_depth_per_gene_sample_wide * gene_mutrate_per_bp_allsamples
    # result # is potentially already a candidate number of synonymous mutations per gene in the sample

    # result2 = num_synonymous * ( result /result.sum() )


    return context_corrected_depth_per_gene_sample_wide




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
        
        mutability_sample = mutability_sample.set_index(["GENE", "CONTEXT_MUT"])
        mutability_all_samples = pd.concat( (mutability_all_samples, mutability_sample), axis = 1)

    return mutability_all_samples.reset_index()



def adapt_mutational_profile(mut_profile_file, samples = None):
    """
    Process a custom mutational profile provided by the user to ensure:
        - Proper format.
        - Any zeros are replaced with a pseudocount equal to the minimum non-zero value, 
            followed by renormalization.
    
    Args:
        mut_profile_file (str): Path to the mutational profile file (tab-separated).
        samples (list): List of sample names to include in the process.

    Returns:
        pd.DataFrame: Processed mutational profile with appropriate pseudocounts and normalization.
    
    Note: this function can receive a dataframe containing a column per sample and should be able to handle it
    """
    # Load the mutational profile
    mut_probability = pd.read_csv(mut_profile_file, sep="\t", header=0, index_col=0)
    mut_probability.columns = [x.split(".")[0] for x in mut_probability.columns]

    if samples is not None:
        # Filter for the specified samples
        mut_probability = mut_probability[samples].copy()

    # Create an empty matrix with all contexts, filling missing entries with zeros
    empty_matrix = pd.DataFrame(index=CHANNELS)
    mut_probability = pd.concat((empty_matrix, mut_probability), axis=1)
    mut_probability = mut_probability.fillna(0)

    # TODO
    # revise that this way of adding a pseudocount makes sense
    # If there are zeros, add a pseudocount
    if (mut_probability == 0).any().any():
        # Add pseudocount: minimum non-zero value
        min_value = mut_probability[mut_probability > 0].min().min()
        pseudocount = min_value if min_value > 0 else 1e-6  # Fallback to a small value if necessary
        logger.info("Adding a pseudocount of {}".format(pseudocount))
        mut_probability += pseudocount

    # Renormalize each column (sample) independently
    mut_probability = mut_probability.div(mut_probability.sum(axis=0), axis=1)

    # Add the context mutation name to the index
    mut_probability.index.name = "CONTEXT_MUT"
    return mut_probability.reset_index()





def compute_mutabilities_wrapper(all_possible_sites_annotated_file,
                                    depth_dataframe_file,
                                    mutations_file, 
                                    table_muts_x_sample_gene_impact_context,
                                    mutability_table,
                                    syn_muts_table,
                                    mut_profile = None,
                                    mut_profile_global = None,
                                    single_sample = None,
                                    absent_synonymous = 'ignore',
                                    gene_mutation_rates_file = None
                                    ):
    """
    Wrapper for all the steps required to compute the mutabilities per sample, gene and context
    """
    # Read files
    all_possible_sites_annotated, depth_dataframe, maf = read_inputs(all_possible_sites_annotated_file,
                                                                        depth_dataframe_file,
                                                                        mutations_file
                                                                        )

    if absent_synonymous == 'infer_global_custom':
        # check if exists
        if not os.path.isfile(gene_mutation_rates_file):
            logger.debug(f"Gene synonymous mutations file : {gene_mutation_rates_file} does not exist")
        if not os.path.isfile(mut_profile_global):
            logger.debug(f"Global mutational profile file : {mut_profile_global} does not exist")
            

    logger.debug("Inputs loaded")
    # Annotate mutations
    annotated_minimal_maf = annotate_mutations_using_vep(maf, all_possible_sites_annotated)
    logger.debug("Mutations annotated")


    # checked until here

    # Define for which samples we have enough data
    if single_sample:
        samples = [single_sample]
        logger.info("Running in single sample mode with {}".format(single_sample))
        logger.warning("Be careful that there is no subselection of mutations when running in single sample mode.")

    else:
        samples = define_samples(depth_dataframe, annotated_minimal_maf)
        logger.info("Samples selected")

        # focus the mutations and depth dataframes into the selected samples
        annotated_minimal_maf = annotated_minimal_maf[annotated_minimal_maf["SAMPLE_ID"].isin(samples)].copy().reset_index(drop = True)
        logger.debug("Mutations subsetted")

    depth_dataframe = depth_dataframe[["CHROM", "POS"] + samples].copy()
    depth_dataframe_small = depth_dataframe[["CHROM", "POS"] + samples].copy()
    depth_dataframe_small[samples] = depth_dataframe_small[samples]/3
    logger.debug("Depths subsetted")


    # compute table of observed mutations
    obs_muts_per_gene_impact_context_sample_wide = compute_mutations_per_sample_gene_impact_context_table(annotated_minimal_maf)
    logger.debug("Mutations table produced")

    if single_sample:
        obs_muts_per_gene_impact_context_sample_wide_indexed = obs_muts_per_gene_impact_context_sample_wide.set_index(["GENE", "IMPACT", "CONTEXT_MUT"])
        obs_muts_per_gene_impact_context_sample_wide = obs_muts_per_gene_impact_context_sample_wide_indexed.sum(axis = 1).reset_index()
        obs_muts_per_gene_impact_context_sample_wide.columns = ["GENE", "IMPACT", "CONTEXT_MUT"] + samples
        logger.debug("Mutations table compressed for single sample")

    ## TODO: We could try to do something similar with the depths
    # if single_sample:
    #     depth_dataframe
    #     obs_muts_per_gene_impact_context_sample_wide = obs_muts_per_gene_impact_context_sample_wide_indexed.sum(axis = 1).reset_index()
    #     obs_muts_per_gene_impact_context_sample_wide.columns = ["GENE", "IMPACT", "CONTEXT_MUT"] + samples
    #     logger.debug("Mutations table compressed for single sample")



    # *** OUTPUT 1 ***
    ## Store table with observed mutations 
    obs_muts_per_gene_impact_context_sample_wide.to_csv(f"{table_muts_x_sample_gene_impact_context}",
                                                        header = True,
                                                        index = False,
                                                        sep = "\t")



    # Compute mutational profile from the input data
    if mut_profile:
        mut_probability = adapt_mutational_profile(mut_profile, samples)
        logger.info("Mutational profile loaded")       

    else:
        mut_probability = compute_mutational_profile(annotated_minimal_maf,
                                                        all_possible_sites_annotated,
                                                        depth_dataframe,
                                                        samples,
                                                        pseudocount = 0.5 # FIXME revise if this value makes sense
                                                        )
        logger.info("Mutational profile computed")


    # mut_probability: contains the weights per context (96 values)



    # I want to revise this one in more detail.
    # Compute expected synonymous mutations
    weights_syn_per_gene_per_sample = compute_expected_synonymous_mutations(all_possible_sites_annotated,
                                                                                depth_dataframe,
                                                                                mut_probability,
                                                                                samples)
    logger.debug("Theoretical expected synonymous computed")



    if absent_synonymous == 'infer_global_custom':
        ## FIXME: this will keep the original samples' name, even multiple columns if single sample is not activated
        mut_probability_global = adapt_mutational_profile(mut_profile_global)

        # this might not work if there are more than one column            
        mut_probability_global.columns = ['CONTEXT_MUT'] + samples
        logger.info("Global mutational profile loaded")

        #######
        ## the following lines follow this correction strategy applied for each gene
        #######
        # syn_mutrate * (depth * (sample_mut_profile/all_samples_mut_profile) )
        #######

        # this term is loaded from the cohort information:
        # syn_mutrate

        # Read mutation rates per MB
        gene_mutation_rates = pd.read_table(gene_mutation_rates_file)
        gene_mutation_rates = gene_mutation_rates.set_index('GENE')["MUTDENSITY"]
        logger.info("Gene mutation rates loaded")

        if len(gene_mutation_rates[gene_mutation_rates == 0]) > 0:
            logger.info("The following genes have a synonymous mutation rate of 0, that needs to be filled.")
            logger.info(gene_mutation_rates[gene_mutation_rates == 0].index)
            # print(gene_mutation_rates.mean())
            mean_syn_mutrate = gene_mutation_rates[gene_mutation_rates != 0].mean()
            logger.info("They will be filled with " + str(mean_syn_mutrate))
            gene_mutation_rates[gene_mutation_rates == 0] = mean_syn_mutrate
            # print(gene_mutation_rates)

        gene_mutation_rates = pd.DataFrame(gene_mutation_rates).reset_index()
        gene_mutation_rates.columns = ["GENE", "MUTDENSITY"]
        gene_mutation_rates = gene_mutation_rates.set_index('GENE')

        # remove the per MB correction
        gene_mutation_rates = gene_mutation_rates / 1e6
        logger.debug(gene_mutation_rates)

        ## FIXME: revise if this should be changed and mutrates should be provided without scaling to Mb

        # this term is computed below starting now.
        # (depth * (sample_mut_profile/all_samples_mut_profile) )
        # we get a value per gene
        sample_specific_biases = compute_sample_gene_specific_differences(all_possible_sites_annotated,
                                                                            depth_dataframe_small,
                                                                            mut_probability_global,
                                                                            mut_probability,                                                
                                                                            samples)
        sample_columnsss = [x for x in sample_specific_biases.columns]
        sample_specific_biases = sample_specific_biases.reset_index()

        logger.info("Sample specific biases computed based on mutational profile and sequencing depth")
        logger.info(sample_specific_biases)

        # Extract the gene name before the -- (used as separator) for subgene-level data
        sample_specific_biases['GENE_BASE'] = sample_specific_biases['GENE'].str.split('--').str[0]

        # Add info from the synoymous mutation rates defined for each gene 
        weighted_depth_n_mutrate = sample_specific_biases.merge(gene_mutation_rates,
                                                                left_on='GENE_BASE',
                                                                right_on='GENE',
                                                                how = 'left'
                                                                ).fillna(0)
        weighted_depth_n_mutrate_ind = weighted_depth_n_mutrate[["GENE", 'MUTDENSITY'] + sample_columnsss].set_index("GENE")
        
        # compute mutation_rate * corrected_depth product
        mutation_numbers = (weighted_depth_n_mutrate_ind.iloc[:,0] * weighted_depth_n_mutrate_ind.iloc[:,1]).reset_index()
        mutation_numbers.columns = ["GENE", "mutations"]
        mutation_numbers = mutation_numbers.merge(weighted_depth_n_mutrate[["GENE", "GENE_BASE"]],
                                                    on = "GENE",
                                                    how = 'left')
        logger.debug("computed synonymous mutation numbers from global mutrate.")
        logger.debug("This is a potential output of the number of synonymous mutations per gene-sample")
        logger.debug(mutation_numbers)
        # mutation_numbers 
        # is a potential output number of mutations per gene per sample
        # using the global synonymous mutation rates


        logger.debug("[not used] real observed synonymous mutations per gene")
        logger.debug(obs_muts_per_gene_impact_context_sample_wide[
                                                                obs_muts_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"
                                                            ].reset_index(drop = True).groupby(by = 'GENE')[samples].sum()
        )


        #####
        ## Alternative option
        ##   keep the total number of observed synonymous mutations constant
        ##   and redistribute them between genes according to the recomputed probabilities
        #####
        mutation_numbers_full_genes = mutation_numbers[mutation_numbers["GENE"] == mutation_numbers["GENE_BASE"]]
        logger.debug("computed synonymous mutations per full gene from global mutrates")
        logger.debug(mutation_numbers_full_genes)


        # Handle the subgenic elements
        # total number of mutations comes only from full genes
        # these are just portions of the bigger ones and need to be handled differently
        mutation_numbers_gene_regions = mutation_numbers[mutation_numbers["GENE"] != mutation_numbers["GENE_BASE"]]
        mutation_numbers_gene_regions_with_gene_info = mutation_numbers_gene_regions.merge(mutation_numbers_full_genes[["GENE_BASE", "mutations"]],
                                                                                           on = 'GENE_BASE',
                                                                                           suffixes = ("_regions", "_full")
                                                                                           )
        mutation_numbers_gene_regions_with_gene_info = mutation_numbers_gene_regions_with_gene_info[
                                                            ["GENE", "GENE_BASE",
                                                             "mutations_regions", "mutations_full"]
                                                             ]

        # compute proportion of mutations seen in subgenic element
        # to be able to normalize the new total with this value
        mutation_numbers_gene_regions_with_gene_info["proportion"] = (mutation_numbers_gene_regions_with_gene_info["mutations_regions"] \
                                                                        / mutation_numbers_gene_regions_with_gene_info["mutations_full"]).fillna(0)
        mutation_numbers_gene_regions_with_gene_info = mutation_numbers_gene_regions_with_gene_info[["GENE", "GENE_BASE","proportion"]]
        # logger.debug("number of mutations")
        # logger.debug(mutation_numbers_gene_regions_with_gene_info)

        mutation_numbers_full_genes = mutation_numbers_full_genes.drop("GENE_BASE",
                                                                       axis='columns').set_index("GENE")


        # convert the total number of mutations "predicted" for each gene-sample to a relative value
        relative_syn_muts_per_gene_allsamples = mutation_numbers_full_genes / mutation_numbers_full_genes.sum()


        # Count how many synonymous mutations are there in each sample irrespective of the gene
        # restrict to full genes only to avoid counting mutations more than once
        full_genes_names = sorted(mutation_numbers_full_genes.index)
        obs_muts_per_gene_impact_context_sample_wide_only_full_genes = obs_muts_per_gene_impact_context_sample_wide[obs_muts_per_gene_impact_context_sample_wide["GENE"].isin(full_genes_names)]
        syn_muts_per_sample = obs_muts_per_gene_impact_context_sample_wide_only_full_genes[
                                                            obs_muts_per_gene_impact_context_sample_wide_only_full_genes["IMPACT"] == "synonymous"
                                                        ].reset_index(drop = True)[samples].sum()
        syn_muts_per_sample_df = pd.DataFrame(syn_muts_per_sample).T

        # FIXME : this multiplication should have some more safety measures
        # Multiply the two "vectors" to get a value of synonymous mutations per sample per gene
        result_array = relative_syn_muts_per_gene_allsamples.values * syn_muts_per_sample_df.values


        # put the right names to the rows and columns
        obs_syn_muts_per_gene_sample = pd.DataFrame(result_array,
                                                            columns= syn_muts_per_sample_df.columns,
                                                            index=relative_syn_muts_per_gene_allsamples.index)
        obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_sample.reset_index()
        logger.debug("computed synonymous mutation numbers per gene after normalizing")
        logger.debug(obs_syn_muts_per_gene_sample)


        # get the numbers including the subgenic regions
        mutation_numbers_gene_regions_with_gene_info_final_numbers = mutation_numbers_gene_regions_with_gene_info.merge(obs_syn_muts_per_gene_sample,
                                                                                                                        left_on = 'GENE_BASE',
                                                                                                                        right_on = 'GENE',                                                                                                                        
                                                                                                                        how = 'left',
                                                                                                                        suffixes = ("", "_gene")
                                                                                                                        )

        # get the numbers of synoymous mutations for the subgenic areas
        mutation_numbers_gene_regions_with_gene_info_final_numbers[samples[0]] = mutation_numbers_gene_regions_with_gene_info_final_numbers[samples[0]] \
                                                                                * mutation_numbers_gene_regions_with_gene_info_final_numbers["proportion"]
        logger.debug("computed synonymous mutation numbers per subgenic region after normalizing")
        logger.debug(mutation_numbers_gene_regions_with_gene_info_final_numbers)

        # concat gene and sub-genic level metrics
        obs_syn_muts_per_gene_sample = pd.concat((obs_syn_muts_per_gene_sample,
                                                  mutation_numbers_gene_regions_with_gene_info_final_numbers[["GENE"] + samples])
                                                  ).reset_index(drop = True)
        logger.debug("Synonynmous computed from the total number of observed synonymous and distributed according to the relative counts in the custom file provided.")
        # print('CV', obs_syn_muts_per_gene_sample)

    elif absent_synonymous == 'infer_covariates':
        pass

    # else:
    elif absent_synonymous == 'ignore':
    # Count of observed synonymous mutations per sample and gene
        obs_syn_muts_per_gene_context_sample = obs_muts_per_gene_impact_context_sample_wide[
                                                        obs_muts_per_gene_impact_context_sample_wide["IMPACT"] == "synonymous"].reset_index(
                                                            drop = True)
        obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_context_sample.groupby(by = ["GENE"])[samples].sum()
        obs_syn_muts_per_gene_sample = obs_syn_muts_per_gene_sample.reset_index()
        # print('LOC', obs_syn_muts_per_gene_sample1)
        logger.debug("Observed synonymous computed")

    if syn_muts_table:
        syn_muts2store = obs_syn_muts_per_gene_sample.copy()
        syn_muts2store = syn_muts2store.groupby(by = ["GENE"]).sum()[samples].sum(axis = 1).reset_index()
#        syn_muts2store.columns = ["GENE", "SYNONYMOUS_MUTS"]
        syn_muts2store.to_csv(f"{syn_muts_table}",
                                header = True,
                                index = False,
                                sep = "\t")
        logger.debug(f"Synonymous mutations used stored into: {syn_muts_table}")


    obs_muts_per_gene_context_sample = obs_muts_per_gene_impact_context_sample_wide.reset_index(drop = True)
    obs_muts_per_gene_sample = obs_muts_per_gene_context_sample.groupby(by = ["GENE"])[samples].sum()
    obs_muts_per_gene_sample = obs_muts_per_gene_sample.reset_index()



    #####
    # Alpha computation
    #####
    # What is alpha?
    # Alpha is a factor that allows us to update mutability per context
    # from relative to absolute in each specific gene-sample pair
    #####

    # we compute the value of alpha per each gene-sample pair,
    # by dividing the number of observed synonymous
    # by the number of synonymous we would be generating with the original mutational profile
    alpha_per_sample = obs_syn_muts_per_gene_sample.set_index("GENE").divide( weights_syn_per_gene_per_sample.set_index("GENE") )
    # print(alpha_per_sample)


    # print("alpha per sample all mutations")
    # print(obs_syn_muts_per_gene_sample.set_index("GENE"))
    # print(obs_muts_per_gene_sample.set_index("GENE"))
    # print(obs_muts_per_gene_sample.set_index("GENE").divide( expected_syn_per_gene_per_sample.set_index("GENE") ))


    # compute the mutabilities by adjusting the mutational profile (mut_probability)
    # by the value of alpha, to obtain an absolute mutability per context
    mutability_all_samples = compute_mutabilities(alpha_per_sample, mut_probability, samples)
    logger.debug("Mutabilities table computed.")

    # *** OUTPUT 2 *** 
    # create a single table with all the mutabilities per sample, gene, context
    # TODO revise what happens here with this fillna(0)
    # I think that it goes to NA when there is no synonymous mutation
    # but we should check properly
    mutability_all_samples.fillna(0).to_csv(f"{mutability_table}",
                                        header = True,
                                        index = False,
                                        sep = "\t")

