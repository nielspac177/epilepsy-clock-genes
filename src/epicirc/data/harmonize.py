"""Harmonize GWAS summary statistics to a common effect allele.

Correct allele harmonization is where two-sample MR and cross-trait analyses most often go
silently wrong. This module aligns a variant's effect estimate to a reference effect allele,
handling straight matches, allele swaps (flip sign), strand flips, and flagging palindromic
(strand-ambiguous) SNPs so downstream steps can drop or MAF-resolve them.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
Action = Literal["match", "flip_sign", "strand_flip", "strand_flip_and_sign", "ambiguous", "incompatible"]


@dataclass(frozen=True)
class Harmonized:
    beta: float
    effect_allele: str
    other_allele: str
    action: Action
    usable: bool


def complement(allele: str) -> str:
    try:
        return "".join(_COMPLEMENT[b] for b in allele.upper())
    except KeyError as e:  # non-ACGT (indels etc.)
        raise ValueError(f"non-ACGT allele: {allele}") from e


def is_palindromic(a1: str, a2: str) -> bool:
    """A/T or C/G variants are strand-ambiguous."""
    return complement(a1.upper()) == a2.upper()


def harmonize(
    beta: float,
    effect_allele: str,
    other_allele: str,
    ref_effect: str,
    ref_other: str,
    drop_palindromic: bool = True,
) -> Harmonized:
    """Align ``beta`` (for ``effect_allele``) to ``ref_effect``.

    Palindromic SNPs are marked unusable when ``drop_palindromic`` (the safe default), because
    without matched allele frequencies their strand cannot be resolved.
    """
    ea, oa = effect_allele.upper(), other_allele.upper()
    re_, ro = ref_effect.upper(), ref_other.upper()

    if is_palindromic(ea, oa):
        return Harmonized(beta, ea, oa, "ambiguous", usable=not drop_palindromic)

    if {ea, oa} == {re_, ro}:
        if ea == re_:
            return Harmonized(beta, ea, oa, "match", usable=True)
        return Harmonized(-beta, re_, ro, "flip_sign", usable=True)

    # try strand flip
    cea, coa = complement(ea), complement(oa)
    if {cea, coa} == {re_, ro}:
        if cea == re_:
            return Harmonized(beta, re_, ro, "strand_flip", usable=True)
        return Harmonized(-beta, re_, ro, "strand_flip_and_sign", usable=True)

    return Harmonized(beta, ea, oa, "incompatible", usable=False)
