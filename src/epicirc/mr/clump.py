"""Distance-based clumping — select approximately independent instruments without an LD panel.

From the MR playbook: with no PLINK/LD reference available, keep the lead SNP in each locus and
drop any SNP within +/- ``window_kb`` on the same chromosome (greedy, strongest-first). This is an
LD *approximation* (a real LD-clump can retain independent SNPs closer than the window and prune
correlated ones farther apart); results must say so. Adequate for a first-pass two-sample MR until
the LD panel is installed.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Variant:
    snp: str
    chrom: str
    pos: int
    p: float


def distance_clump(variants: list[Variant], window_kb: int = 1000) -> list[Variant]:
    """Greedy distance clumping: strongest SNP first, remove neighbours within the window."""
    window = window_kb * 1000
    kept: list[Variant] = []
    for v in sorted(variants, key=lambda x: x.p):          # strongest (smallest p) first
        if all(not (v.chrom == k.chrom and abs(v.pos - k.pos) <= window) for k in kept):
            kept.append(v)
    return kept
