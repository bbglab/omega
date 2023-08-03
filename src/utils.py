from itertools import product


def canonical_channels(): 

    subs = [''.join(z) for z in product('CT', 'ACGT') if z[0] != z[1]]
    flanks = [''.join(z) for z in product('ACGT', repeat=2)]
    contexts_tuples = [(a, b) for a, b in product(subs, flanks)]
    sorted_contexts_tuples = sorted(contexts_tuples, key=lambda x: (x[0], x[1]))
    sorted_contexts = [b[0] + a[0] + b[1] + '>' + a[1] for a, b in sorted_contexts_tuples]
    
    return sorted_contexts


def dict_append(d1, d2):
    
    return {k: d1.get(k, []) + d2.get(k, []) for k in d2}
    