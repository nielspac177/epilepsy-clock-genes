"""Tier-1 engine tests: cache path reproduces geneset_enrich, and LOO/diff behave sanely."""
import random

from epicirc.analysis import geneset_enrich as ge
from epicirc.analysis.tier1 import diff_test, matched_null


def _synthetic():
    """A deterministic top/length/nsnp universe with a planted 'clock' signal."""
    rng = random.Random(7)
    top = {f"G{i}": abs(rng.gauss(3, 2)) for i in range(400)}
    length = {f"G{i}": 1000 + 37 * i for i in range(400)}
    nsnp = {f"G{i}": 5 + (i % 20) for i in range(400)}
    present = [f"G{i}" for i in range(0, 60, 3)]
    for g in present:                 # plant an enrichment
        top[g] += 4.0
    return top, length, nsnp, present


def test_tier1_matched_null_matches_geneset_enrich_bit_for_bit():
    """The refactor must not move the published numbers: identical RNG => identical output."""
    top, length, nsnp, present = _synthetic()
    a = matched_null(top, length, nsnp, present, n_null=1000, seed=42)
    b = ge.matched_null(top, length, nsnp, present, n_null=1000, seed=42)
    for k in ("obs", "null_mean", "emp_p", "obs_lo", "obs_hi", "null_lo", "null_hi"):
        assert a[k] == b[k], f"{k}: {a[k]} != {b[k]}"


def test_leave_one_out_never_raises_and_shrinks_k():
    top, length, nsnp, present = _synthetic()
    full = matched_null(top, length, nsnp, present, n_null=300, seed=1)
    for drop in present:
        sub = [g for g in present if g != drop]
        r = matched_null(top, length, nsnp, sub, n_null=300, seed=1)
        assert r["k"] == full["k"] - 1
        assert 0 < r["emp_p"] <= 1


def test_diff_test_symmetry_and_sign():
    """Same phenotype vs itself => Δ≈0; a stronger A => Δ>0."""
    top, length, nsnp, present = _synthetic()
    same = diff_test(top, length, nsnp, top, length, nsnp, present, n_null=500, seed=3)
    assert abs(same["obs_delta"]) < 1e-9        # identical inputs cancel exactly
    # "focal": remove the planted clock enrichment so B's set sits at baseline (ratio≈1),
    # while A keeps it. Ratio is scale-invariant, so the contrast must be set-specific.
    top_b = dict(top)
    for g in present:
        top_b[g] -= 4.0
    d = diff_test(top, length, nsnp, top_b, length, nsnp, present, n_null=500, seed=3)
    assert d["ratio_a"] > d["ratio_b"]          # A enriched relative to B
    assert d["obs_delta"] > 0
