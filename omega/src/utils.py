"""This module contains utility functions for the omega package."""
from typing import Callable
from itertools import product

from bgreference import hg38, hg19, mm10, mm39

# region context_store

ASSEMBLY: dict[str, Callable] = {
    'hg38': hg38,
    'hg19': hg19,
    'mm10': mm10,
    'mm39': mm39
}

cb = dict(zip('ACGT', 'TGCA'))


def canonical_channels() -> list[str]:
    """Return the list of canonical mutation contexts in a sorted order."""
    subs: list = [''.join(z) for z in product('CT', 'ACGT') if z[0] != z[1]]
    flanks: list = [''.join(z) for z in product('ACGT', repeat=2)]
    contexts_tuples: list[tuple[str, str]] = [(a, b) for a, b in product(subs, flanks)]
    sorted_contexts_tuples: list[tuple[str, str]] = sorted(contexts_tuples, key=lambda x: (x[0], x[1]))
    sorted_contexts: list[str] = [b[0] + a[0] + b[1] + '>' + a[1] for a, b in sorted_contexts_tuples]

    return sorted_contexts


def transform_context(chr_: str, pos: int, mut: str, assembly: str = "hg38") -> str:
    """Convert a chromosome, position, and mutation to a string context (e.g., ACG>AT)."""
    genome_assembly = ASSEMBLY.get(assembly, hg38)
    _ref, alt = tuple(mut.split('/'))
    ref_triplet = genome_assembly(chr_, pos-1, size=3)
    if ref_triplet[1] not in ['C', 'T']:
        ref_triplet = ''.join(list(map(lambda x: cb[x], ref_triplet[::-1]))) # TODO avoid map
        alt = cb[alt]
    return ref_triplet + '>' + alt


def transform_bracket_context(bracket_context: str) -> str:
    """Convert a bracket context (e.g., A[C>T]G) to a string context (e.g., ACG>AT)."""
    # FIXME: This function is never called
    ref = bracket_context[2]
    alt = bracket_context[4]
    flank1 = bracket_context[0]
    flank2 = bracket_context[-1]
    return flank1 + ref + flank2 + '>' + alt

# endregion context_store