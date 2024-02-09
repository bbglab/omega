from itertools import product

from bgreference import hg38, hg19, mm10


cb = dict(zip('ACGT', 'TGCA'))


def canonical_channels():

    subs = [''.join(z) for z in product('CT', 'ACGT') if z[0] != z[1]]
    flanks = [''.join(z) for z in product('ACGT', repeat=2)]
    contexts_tuples = [(a, b) for a, b in product(subs, flanks)]
    sorted_contexts_tuples = sorted(contexts_tuples, key=lambda x: (x[0], x[1]))
    sorted_contexts = [b[0] + a[0] + b[1] + '>' + a[1] for a, b in sorted_contexts_tuples]
    return sorted_contexts


def transform_context(chr_, pos, mut, assembly = hg38):
    ref, alt = tuple(mut.split('/'))
    ref_triplet = assembly(chr_, pos-1, size=3)
    if ref_triplet[1] not in ['C', 'T']:
        ref_triplet = ''.join(list(map(lambda x: cb[x], ref_triplet[::-1])))
        alt = cb[alt]
    return ref_triplet + '>' + alt