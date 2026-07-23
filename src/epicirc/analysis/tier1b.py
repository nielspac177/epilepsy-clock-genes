"""Tier-1b hardening: the analyses the adversarial audit demanded before trusting the headline.

  neighbours  : for each clock gene, is the top-SNP inside the gene body or the flank, and which
                OTHER genes' bodies contain that SNP? Then test whether the neighbour is itself
                GGE-specific (its own GGE vs focal top-chi2). Distinguishes "circadian gene" from
                "circadian locus / non-clock neighbour".
  robust      : competitive matched-null using the MEDIAN top-chi2 (outlier-insensitive) and a
                rank-sum statistic, as a cross-check on the mean-of-chi2 headline.
  leavekout   : systematic worst-case leave-k-out (k=1,2,3) — how few genes carry significance.
  multiseed   : the headline enrichment + difference test across many seeds -> median & MC-SD,
                so no single lucky seed is reported.
"""
from __future__ import annotations

import argparse
import random
import statistics
import sys
from itertools import combinations

from epicirc.analysis.tier1 import CLOCK, load_cache, _strata, matched_null, diff_test


# ----------------------------------------------------------------------------- neighbours
def load_gene_windows(genes_path="reference/genes_hg19.tsv"):
    by_chr = {}
    with open(genes_path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                s, e = int(f[2]), int(f[3])
            except ValueError:
                continue
            by_chr.setdefault(f[1], []).append((s, e, f[0]))
    return by_chr


def cmd_neighbours(a):
    """Report body/flank status of each clock top-SNP and the genes overlapping it."""
    gge = {}
    with open(a.cache_gge) as fh:
        next(fh)
        for ln in fh:
            g, c, s, e, L, n, chi, bp = ln.rstrip("\n").split("\t")
            gge[g] = (c, int(s), int(e), float(chi), int(bp))
    top_focal, *_ = load_cache(a.cache_focal)
    by_chr = load_gene_windows(a.genes)
    hdr = ("clock_gene\tchr\ttop_bp\ttop_chi2\tin_body\tdist_to_body_kb\t"
           "overlapping_genes\tgge_specific_neighbour")
    rows = []
    for g in CLOCK:
        if g not in gge:
            continue
        c, s, e, chi, bp = gge[g]
        in_body = s <= bp <= e
        dist = 0 if in_body else round(min(abs(bp - s), abs(bp - e)) / 1000, 1)
        # genes whose BODY contains this SNP (true co-location, not just window)
        overlap = [gn for (gs, ge_, gn) in by_chr.get(c, []) if gs <= bp <= ge_]
        others = [gn for gn in overlap if gn != g]
        # is any non-clock overlapping neighbour itself GGE-specific? (its GGE chi >> focal chi)
        flags = []
        for gn in others:
            gchi = gge.get(gn, (None,)*5)[3]
            fchi = top_focal.get(gn)
            if gchi is not None and fchi is not None:
                flags.append(f"{gn}(GGE {gchi:.1f}/focal {fchi:.1f})")
        rows.append((g, c, bp, f"{chi:.2f}", "YES" if in_body else "no", dist,
                     ";".join(others) if others else "-", " | ".join(flags) if flags else "-"))
    _write(a.out, hdr, rows)
    nbody = sum(1 for r in rows if r[4] == "YES")
    print(f"[neighbours] {nbody}/{len(rows)} clock top-SNPs in-body; {len(rows)-nbody} in flank", file=sys.stderr)
    for r in rows:
        print(f"  {r[0]:9s} chi={r[3]:>6} in_body={r[4]:3s} overlap={r[6]}", file=sys.stderr)


# ----------------------------------------------------------------------------- robust stat
def _matched_null_stat(top, length, nsnp, present, n_null, seed, statfn):
    strata, len_dec, snp_dec = _strata(top, length, nsnp)
    rng = random.Random(seed)
    obs = statfn([top[g] for g in present])
    ge = 0
    for _ in range(n_null):
        pick = [rng.choice(strata[(len_dec[g], snp_dec[g])]) for g in present]
        if statfn([top[g] for g in pick]) >= obs:
            ge += 1
    return obs, (ge + 1) / (n_null + 1)


def cmd_ranktest(a):
    """Scale-free, tail-robust competitive test: mean percentile-rank of clock genes vs matched null.

    Converts each gene's top-chi2 to its genome-wide percentile rank (uniform[0,1] marginally), so
    cross-cohort differences in test-statistic scale / genomic inflation / heavy tails cancel. The
    matched null still controls length and SNP count. This is the correct statistic for replicating
    across cohorts whose raw chi2 distributions differ (e.g. ILAE meta-Z vs FinnGen mixed-model).
    """
    top, length, nsnp = load_cache(a.cache)
    genes = [g for g in top if g in length]
    ranked = sorted(genes, key=lambda g: top[g])
    n = len(ranked)
    pct = {g: (i + 0.5) / n for i, g in enumerate(ranked)}   # midpoint percentile
    present = [g for g in CLOCK if g in top]
    obs, p = _matched_null_stat(pct, length, nsnp, present, a.n_null, a.seed, statistics.mean)
    hdr = "phenotype\tstatistic\tmean_percentile\temp_p\tn_genes"
    _write(a.out, hdr, [(a.label, "rank", f"{obs:.4f}", f"{p:.3e}", len(present))])
    print(f"[ranktest] {a.label}: mean clock percentile={obs:.4f} (0.5=null) emp_p={p:.3e}", file=sys.stderr)


def cmd_robust(a):
    top, length, nsnp = load_cache(a.cache)
    present = [g for g in CLOCK if g in top]
    hdr = "phenotype\tstatistic\tobs\temp_p\tn_genes"
    rows = []
    for name, fn in [("mean", statistics.mean), ("median", statistics.median)]:
        obs, p = _matched_null_stat(top, length, nsnp, present, a.n_null, a.seed, fn)
        rows.append((a.label, name, f"{obs:.3f}", f"{p:.3e}", len(present)))
        print(f"[robust] {a.label} {name}: obs={obs:.3f} emp_p={p:.3e}", file=sys.stderr)
    _write(a.out, hdr, rows)


# ----------------------------------------------------------------------------- leave-k-out
def cmd_leavekout(a):
    top, length, nsnp = load_cache(a.cache)
    present = [g for g in CLOCK if g in top]
    hdr = "phenotype\tk_dropped\tworst_dropped_set\tworst_ratio\tworst_emp_p\tstill_sig_0.05"
    rows = []
    for k in (0, 1, 2, 3):
        worst_p, worst = -1.0, None
        combos = [()] if k == 0 else list(combinations(present, k))
        # for k>=2 the search space is large; cap by pre-ranking genes by contribution
        if k >= 2:
            ranked = sorted(present, key=lambda g: top[g], reverse=True)[:8]  # top-8 contributors
            combos = list(combinations(ranked, k))
        for drop in combos:
            sub = [g for g in present if g not in drop]
            r = matched_null(top, length, nsnp, sub, a.n_null, a.seed)
            ratio = r["obs"] / r["null_mean"]
            if r["emp_p"] > worst_p:
                worst_p, worst = r["emp_p"], (drop, ratio)
        drop, ratio = worst
        rows.append((a.label, k, "+".join(drop) if drop else "<none>", f"{ratio:.3f}",
                     f"{worst_p:.3e}", "YES" if worst_p < 0.05 else "NO"))
        print(f"[leavekout] k={k} worst: drop {drop} -> ratio {ratio:.3f} emp_p {worst_p:.3e} "
              f"{'sig' if worst_p<0.05 else 'NS'}", file=sys.stderr)
    _write(a.out, hdr, rows)


# ----------------------------------------------------------------------------- multi-seed
def cmd_multiseed(a):
    top_g, len_g, ns_g = load_cache(a.cache_gge)
    top_f, len_f, ns_f = load_cache(a.cache_focal)
    present = [g for g in CLOCK if g in top_g]
    seeds = list(range(a.n_seeds))
    enr_p, diff_p, ratios, deltas = [], [], [], []
    for s in seeds:
        r = matched_null(top_g, len_g, ns_g, present, a.n_null, s)
        enr_p.append(r["emp_p"]); ratios.append(r["obs"] / r["null_mean"])
        d = diff_test(top_g, len_g, ns_g, top_f, len_f, ns_f, present, a.n_null, s)
        diff_p.append(d["emp_p_2sided"]); deltas.append(d["obs_delta"])
    hdr = "quantity\tmedian\tmc_sd\tmin\tmax\tn_seeds"
    rows = [
        ("enrich_ratio", f"{statistics.median(ratios):.3f}", f"{statistics.pstdev(ratios):.4f}",
         f"{min(ratios):.3f}", f"{max(ratios):.3f}", a.n_seeds),
        ("enrich_emp_p", f"{statistics.median(enr_p):.2e}", f"{statistics.pstdev(enr_p):.2e}",
         f"{min(enr_p):.2e}", f"{max(enr_p):.2e}", a.n_seeds),
        ("diff_delta", f"{statistics.median(deltas):.3f}", f"{statistics.pstdev(deltas):.4f}",
         f"{min(deltas):.3f}", f"{max(deltas):.3f}", a.n_seeds),
        ("diff_emp_p", f"{statistics.median(diff_p):.2e}", f"{statistics.pstdev(diff_p):.2e}",
         f"{min(diff_p):.2e}", f"{max(diff_p):.2e}", a.n_seeds),
    ]
    _write(a.out, hdr, rows)
    for r in rows:
        print(f"[multiseed] {r[0]}: median={r[1]} mc_sd={r[2]} range=[{r[3]},{r[4]}]", file=sys.stderr)


# ----------------------------------------------------------------------------- io
def _write(path, header, rows):
    from pathlib import Path
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    n = sub.add_parser("neighbours")
    n.add_argument("--cache-gge", required=True); n.add_argument("--cache-focal", required=True)
    n.add_argument("--genes", default="reference/genes_hg19.tsv"); n.add_argument("--out", required=True)
    n.set_defaults(func=cmd_neighbours)

    rk = sub.add_parser("ranktest")
    rk.add_argument("--cache", required=True); rk.add_argument("--label", required=True)
    rk.add_argument("--out", required=True); rk.add_argument("--n-null", type=int, default=10000)
    rk.add_argument("--seed", type=int, default=20260702); rk.set_defaults(func=cmd_ranktest)

    r = sub.add_parser("robust")
    r.add_argument("--cache", required=True); r.add_argument("--label", required=True)
    r.add_argument("--out", required=True); r.add_argument("--n-null", type=int, default=10000)
    r.add_argument("--seed", type=int, default=20260702); r.set_defaults(func=cmd_robust)

    l = sub.add_parser("leavekout")
    l.add_argument("--cache", required=True); l.add_argument("--label", required=True)
    l.add_argument("--out", required=True); l.add_argument("--n-null", type=int, default=10000)
    l.add_argument("--seed", type=int, default=20260702); l.set_defaults(func=cmd_leavekout)

    m = sub.add_parser("multiseed")
    m.add_argument("--cache-gge", required=True); m.add_argument("--cache-focal", required=True)
    m.add_argument("--out", required=True); m.add_argument("--n-null", type=int, default=10000)
    m.add_argument("--n-seeds", type=int, default=50); m.set_defaults(func=cmd_multiseed)

    a = ap.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
