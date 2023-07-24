# bayes run cli
# --input: table with columns chr | pos | ref | mut | context | gene | impact | samples (reads) | samples (mutabilities)
# --impact: grouping of impacts json dict
# --genes: grouping of genes json dict
# --samples: grouping of samples json dict

# mle run cli
# --input: table with columns chr | pos | ref | mut | context | gene | impact | samples (reads) | samples (mutabilities)
# --impact: grouping of impacts json dict
# --genes: grouping of genes json dict
# --samples: grouping of samples json dict

# notes: 
# use typer for cli

# input mutations should be given in some canonical order, e.g.
# chr, pos, alt

# mutability col in the input
# is calculated using the no. syn mutations and depth per position
