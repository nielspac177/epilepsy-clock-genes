"""Multiple-testing control: Benjamini-Hochberg FDR and Bonferroni.

Pure-Python (no scipy/statsmodels) so the statistical core is dependency-light and unit-tested.
The pipeline reports BOTH: Bonferroni on the pre-registered primary-contrast family, BH-FDR on
the exploratory subtype x trait x method grid (see config/analysis_params.yaml).
"""
from __future__ import annotations

from typing import Sequence


def bonferroni(pvals: Sequence[float]) -> list[float]:
    """Bonferroni-adjusted p-values (capped at 1.0)."""
    m = len(pvals)
    _check(pvals)
    return [min(1.0, p * m) for p in pvals]


def benjamini_hochberg(pvals: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values (q-values), preserving input order.

    Uses the standard step-up monotonicity enforcement so adjusted values are non-decreasing
    in the original p-value ranking.
    """
    _check(pvals)
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    # walk from largest p-value (rank m) down to smallest (rank 1)
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = pvals[i] * m / rank
        prev = min(prev, val)
        adj[i] = min(1.0, prev)
    return adj


def _check(pvals: Sequence[float]) -> None:
    for p in pvals:
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"p-value out of range: {p}")
