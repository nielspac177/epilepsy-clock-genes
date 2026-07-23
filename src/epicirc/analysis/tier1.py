"""Tier-1 hardening analyses on the gene-stat cache.

Consumes results/gene_stats/<pheno>.tsv (from gene_stats.py) so a single GWAS scan serves every
analysis. The core matched_null reproduces analysis/geneset_enrich.py bit-identically (same RNG
call order) — validated in tests — then adds:
  * leave-one-out    (is the set-level signal carried by one gene, e.g. PER1?)
  * arbitrary sets   (positive / housekeeping / ion-channel / GO-circadian controls)
  * difference test  (GGE vs focal on the SAME matched-null sets -> Δ ratio, shared-LD absorbed)
"""
from __future__ import annotations

import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path

CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]


def load_cache(path: str):
    """Return top (gene->chi2), length (gene->bp), nsnp (gene->count)."""
    top, length, nsnp = {}, {}, {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            g, c, s, e, L, n, chi, bp = ln.rstrip("\n").split("\t")
            top[g] = float(chi)
            length[g] = int(L)
            nsnp[g] = int(n)
    return top, length, nsnp


def _strata(top, length, nsnp):
    """Decile strata by (length, snp-count) — identical construction to geneset_enrich."""
    genes = [g for g in top if g in length]
    by_len = sorted(genes, key=lambda g: length[g])
    by_snp = sorted(genes, key=lambda g: nsnp.get(g, 0))
    n = len(genes)
    len_dec = {g: (10 * i // n) for i, g in enumerate(by_len)}
    snp_dec = {g: (10 * i // n) for i, g in enumerate(by_snp)}
    strata = defaultdict(list)
    for g in genes:
        strata[(len_dec[g], snp_dec[g])].append(g)
    return strata, len_dec, snp_dec


def matched_null(top, length, nsnp, present, n_null, seed, boot=2000):
    """Competitive covariate-matched null. Reproduces geneset_enrich.matched_null exactly."""
    strata, len_dec, snp_dec = _strata(top, length, nsnp)
    rng = random.Random(seed)
    k = len(present)
    obs = sum(top[g] for g in present) / k
    ge = 0
    null_means = []
    for _ in range(n_null):
        pick = [rng.choice(strata[(len_dec[g], snp_dec[g])]) for g in present]
        m = sum(top[g] for g in pick) / len(pick)
        null_means.append(m)
        if m >= obs:
            ge += 1
    null_mean = sum(null_means) / len(null_means)
    null_means.sort()
    null_lo = null_means[int(0.025 * n_null)]
    null_hi = null_means[int(0.975 * n_null)]
    vals = [top[g] for g in present]
    b = sorted(sum(rng.choice(vals) for _ in range(k)) / k for _ in range(boot))
    obs_lo, obs_hi = b[int(0.025 * boot)], b[int(0.975 * boot)]
    return {"obs": obs, "obs_lo": obs_lo, "obs_hi": obs_hi, "null_mean": null_mean,
            "null_lo": null_lo, "null_hi": null_hi, "emp_p": (ge + 1) / (n_null + 1), "k": k}


def diff_test(top_a, len_a, nsnp_a, top_b, len_b, nsnp_b, present, n_null, seed):
    """Δ = ratio_A - ratio_B on the SAME matched-null gene sets, sampled from A's strata.

    Because each null draw is scored in both phenotypes with identical genes, shared controls and
    shared LD structure cancel in Δ (the difference test is thereby conservative, not liberal).
    """
    # Null universe = genes present in BOTH caches, so every sampled null gene is scorable in
    # both phenotypes (their gene coverage differs slightly due to per-window SNP presence).
    common = {g: top_a[g] for g in top_a if g in top_b}
    len_c = {g: len_a[g] for g in common}
    nsnp_c = {g: nsnp_a[g] for g in common}
    strata, len_dec, snp_dec = _strata(common, len_c, nsnp_c)
    present = [g for g in present if g in top_a and g in top_b]
    k = len(present)
    obs_a = sum(top_a[g] for g in present) / k
    obs_b = sum(top_b[g] for g in present) / k
    rng = random.Random(seed)
    da = db = 0.0  # accumulate null means for ratio denominators
    null_a, null_b, null_d = [], [], []
    for _ in range(n_null):
        pick = [rng.choice(strata[(len_dec[g], snp_dec[g])]) for g in present]
        ma = sum(top_a[g] for g in pick) / k
        mb = sum(top_b[g] for g in pick) / k
        null_a.append(ma); null_b.append(mb)
        da += ma; db += mb
    nm_a, nm_b = da / n_null, db / n_null
    ratio_a, ratio_b = obs_a / nm_a, obs_b / nm_b
    obs_delta = ratio_a - ratio_b
    # null Δ: ratio of each null draw against the mean null (matched sets in both phenotypes)
    ge = 0
    for ma, mb in zip(null_a, null_b):
        d = ma / nm_a - mb / nm_b
        null_d.append(d)
        if abs(d) >= abs(obs_delta):
            ge += 1
    return {"obs_a": obs_a, "obs_b": obs_b, "ratio_a": ratio_a, "ratio_b": ratio_b,
            "obs_delta": obs_delta, "emp_p_2sided": (ge + 1) / (n_null + 1), "k": k}


def cmd_run(a):
    top, length, nsnp = load_cache(a.cache)
    present = _read_set(a.geneset, top)
    r = matched_null(top, length, nsnp, present, a.n_null, a.seed)
    ratio = r["obs"] / r["null_mean"]
    lo, hi = r["obs_lo"] / r["null_mean"], r["obs_hi"] / r["null_mean"]
    _append(a.out, "phenotype\tgeneset\tn_genes\tobs\tobs_lo\tobs_hi\tnull_mean\tnull_lo\tnull_hi\tratio\tratio_lo\tratio_hi\temp_p",
            f"{a.label}\t{a.set_name}\t{r['k']}\t{r['obs']:.3f}\t{r['obs_lo']:.3f}\t{r['obs_hi']:.3f}\t"
            f"{r['null_mean']:.3f}\t{r['null_lo']:.3f}\t{r['null_hi']:.3f}\t{ratio:.3f}\t{lo:.3f}\t{hi:.3f}\t{r['emp_p']:.3e}")
    print(f"[tier1] {a.label} / {a.set_name}: ratio={ratio:.3f} [{lo:.3f},{hi:.3f}] emp_p={r['emp_p']:.3e} (k={r['k']})", file=sys.stderr)


def cmd_loo(a):
    top, length, nsnp = load_cache(a.cache)
    present = _read_set(a.geneset, top)
    full = matched_null(top, length, nsnp, present, a.n_null, a.seed)
    full_ratio = full["obs"] / full["null_mean"]
    hdr = "phenotype\tdropped\tn_genes\tratio\temp_p\tflips_ns_at_0.05"
    _append(a.out, hdr, f"{a.label}\t<none>\t{full['k']}\t{full_ratio:.3f}\t{full['emp_p']:.3e}\t-")
    for drop in present:
        sub = [g for g in present if g != drop]
        r = matched_null(top, length, nsnp, sub, a.n_null, a.seed)
        ratio = r["obs"] / r["null_mean"]
        flips = "YES" if (full["emp_p"] < 0.05 and r["emp_p"] >= 0.05) else "no"
        _append(a.out, hdr, f"{a.label}\t{drop}\t{r['k']}\t{ratio:.3f}\t{r['emp_p']:.3e}\t{flips}")
        print(f"[loo] drop {drop:8s}: ratio={ratio:.3f} emp_p={r['emp_p']:.3e} {flips}", file=sys.stderr)


def cmd_diff(a):
    ta, la, na = load_cache(a.cache_a)
    tb, lb, nb = load_cache(a.cache_b)
    present = _read_set(a.geneset, ta)
    r = diff_test(ta, la, na, tb, lb, nb, present, a.n_null, a.seed)
    _append(a.out, "geneset\tpheno_a\tpheno_b\tn_genes\tratio_a\tratio_b\tdelta\temp_p_2sided",
            f"{a.set_name}\t{a.label_a}\t{a.label_b}\t{r['k']}\t{r['ratio_a']:.3f}\t{r['ratio_b']:.3f}\t{r['obs_delta']:.3f}\t{r['emp_p_2sided']:.3e}")
    print(f"[diff] {a.label_a} vs {a.label_b}: Δ={r['obs_delta']:.3f} (ratio {r['ratio_a']:.3f} vs {r['ratio_b']:.3f}) p2={r['emp_p_2sided']:.3e}", file=sys.stderr)


def _read_set(path, top):
    if path == "CLOCK":
        syms = CLOCK
    else:
        syms = []
        with open(path) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ln.lower().startswith("symbol"):
                    continue
                syms.append(ln.split("\t")[0])
    present = [g for g in syms if g in top]
    missing = [g for g in syms if g not in top]
    if missing:
        print(f"[tier1] {len(missing)} set genes absent from cache: {missing}", file=sys.stderr)
    return present


def _append(path, header, row):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    new = not p.exists()
    with p.open("a") as fh:
        if new:
            fh.write(header + "\n")
        fh.write(row + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    common_seed = dict(type=int, default=20260702)

    r = sub.add_parser("run"); r.add_argument("--cache", required=True); r.add_argument("--geneset", default="CLOCK")
    r.add_argument("--set-name", default="circadian"); r.add_argument("--label", required=True)
    r.add_argument("--out", required=True); r.add_argument("--n-null", type=int, default=10000); r.add_argument("--seed", **common_seed)
    r.set_defaults(func=cmd_run)

    l = sub.add_parser("loo"); l.add_argument("--cache", required=True); l.add_argument("--geneset", default="CLOCK")
    l.add_argument("--label", required=True); l.add_argument("--out", required=True)
    l.add_argument("--n-null", type=int, default=10000); l.add_argument("--seed", **common_seed)
    l.set_defaults(func=cmd_loo)

    d = sub.add_parser("diff"); d.add_argument("--cache-a", required=True); d.add_argument("--cache-b", required=True)
    d.add_argument("--geneset", default="CLOCK"); d.add_argument("--set-name", default="circadian")
    d.add_argument("--label-a", required=True); d.add_argument("--label-b", required=True)
    d.add_argument("--out", required=True); d.add_argument("--n-null", type=int, default=10000); d.add_argument("--seed", **common_seed)
    d.set_defaults(func=cmd_diff)

    a = ap.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
