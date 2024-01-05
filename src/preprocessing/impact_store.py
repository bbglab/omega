GROUPING_DICT = {

    'transcript_ablation': 'nonsense',
    'splice_acceptor_variant': 'nonsense',
    'splice_donor_variant': 'nonsense',
    'stop_gained': 'nonsense',
    'frameshift_variant': 'nonsense',
    'stop_lost': 'nonsense',
    'start_lost': 'nonsense',

    'missense_variant': 'missense',
    'inframe_insertion': 'missense',
    'inframe_deletion': 'missense',

    'splice_donor_variant': 'essential_splice',
    'splice_acceptor_variant': 'essential_splice',


    'splice_region_variant': 'splice_region',
    'splice_donor_5th_base_variant': 'splice_region',
    'splice_donor_region_variant': 'splice_region',
    'splice_polypyrimidine_tract_variant': 'splice_region',

    'synonymous_variant': 'synonymous',
    'incomplete_terminal_codon_variant': 'synonymous',
    'start_retained_variant': 'synonymous',
    'stop_retained_variant': 'synonymous',

    'mature_miRNA_variant': 'non_coding_exon_region',
    '5_prime_UTR_variant': 'non_coding_exon_region',
    '3_prime_UTR_variant': 'non_coding_exon_region',
    'non_coding_transcript_exon_variant': 'non_coding_exon_region',

    'NMD_transcript_variant': 'non_coding_exon_region',

    'intron_variant': 'intron_variant',

    'non_coding_transcript_variant' : 'non_coding_transcript_variant',
    'upstream_gene_variant': 'non_genic_variant',
    'downstream_gene_variant': 'non_genic_variant',
    'TFBS_ablation': 'non_genic_variant',
    'TFBS_amplification': 'non_genic_variant',
    'TF_binding_site_variant': 'non_genic_variant',
    'regulatory_region_ablation': 'non_genic_variant',
    'regulatory_region_amplification': 'non_genic_variant',
    'feature_elongation': 'non_genic_variant',
    'regulatory_region_variant': 'non_genic_variant',
    'feature_truncation': 'non_genic_variant',
    'intergenic_variant': 'non_genic_variant',

    'coding_sequence_variant' : 'coding_sequence_variant'

}


CONSEQUENCES_LIST = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'splice_donor_5th_base_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
]

consequence_rank_dict = {consequence : rank for rank, consequence in enumerate(CONSEQUENCES_LIST)}
rank_consequence_dict = {rank : consequence for rank, consequence in enumerate(CONSEQUENCES_LIST)}


def most_deleterious(impact_vep_string):
    all_consequences = impact_vep_string.split(",")
    all_consequences_ranks = map(lambda x: consequence_rank_dict[x], all_consequences)
    return rank_consequence_dict[min(all_consequences_ranks)]