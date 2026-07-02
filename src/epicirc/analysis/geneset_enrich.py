"""Competitive, covariate-matched, LD-robust gene-set enrichment (rigorous Aim 1).

Upgrades the mean-chi2 SNP screen two ways the adversarial review demanded:
  * LD-robust: one statistic per gene = the top-SNP chi-square in the gene window (±flank), so a
    gene contributes a single ~independent value instead of many LD-correlated SNPs.
  * Competitive + covariate-matched null: the circadian set's mean top-chi2 is compared to 10,000
    random gene sets matched on gene length and SNP count (deciles), giving an empirical p that is
    calibrated against the pipeline's own false-positive rate — not a genome-wide mean.

Genome-wide gene universe: GENCODE v19 protein-coding (hg19). Still an approximation to full MAGMA
(no within-gene LD model, top-SNP only), but removes the two biggest confounds of the screen.
"""
from __future__ import annotations

import argparse
import random
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

IL = {"CHR": 0, "BP": 1, "Z": 8}
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]


def load_genes(path: str, flank_kb: int):
    """Return per-chr sorted (win_start, win_end, gene) and gene->length."""
    flank = flank_kb * 1000
    by_chr = defaultdict(list)
    length = {}
    with open(path) as fh:
        for ln in fh:
            g, c, s, e = ln.split("\t")
            s, e = int(s), int(e)
            by_chr[c].append((s - flank, e + flank, g))
            length[g] = e - s
    for c in by_chr:
        by_chr[c].sort()
    return by_chr, length


def gene_top_chi2(ilae_tbl: str, by_chr) -> tuple[dict, dict]:
    """Index SNPs per chr, then per gene take the top-SNP chi2 in its window (LD-robust)."""
    pos_by_chr = defaultdict(list)
    chi_by_chr = defaultdict(list)
    with open(ilae_tbl) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) < 10:
                continue
            c = f[IL["CHR"]]
            if c not in by_chr:
                continue
            try:
                bp = int(f[IL["BP"]]); z = float(f[IL["Z"]])
            except ValueError:
                continue
            pos_by_chr[c].append(bp)
            chi_by_chr[c].append(z * z)
    # sort each chr's SNPs by position (carry chi2 alongside)
    for c in pos_by_chr:
        order = sorted(range(len(pos_by_chr[c])), key=lambda i: pos_by_chr[c][i])
        pos_by_chr[c] = [pos_by_chr[c][i] for i in order]
        chi_by_chr[c] = [chi_by_chr[c][i] for i in order]

    top, nsnp = {}, {}
    for c, windows in by_chr.items():
        pos, chi = pos_by_chr.get(c, []), chi_by_chr.get(c, [])
        if not pos:
            continue
        for ws, we, g in windows:
            lo = bisect_left(pos, ws)
            hi = bisect_right(pos, we)
            if hi > lo:
                seg = chi[lo:hi]
                top[g] = max(seg)
                nsnp[g] = hi - lo
    return top, nsnp


def _decile_key(length_rank, snp_rank, nbins=10):
    return (length_rank * nbins + snp_rank)


def matched_null(top, length, nsnp, clock_present, n_null, seed):
    genes = [g for g in top if g in length]
    # rank genes into deciles by length and by snp count
    by_len = sorted(genes, key=lambda g: length[g])
    by_snp = sorted(genes, key=lambda g: nsnp.get(g, 0))
    n = len(genes)
    len_dec = {g: (10 * i // n) for i, g in enumerate(by_len)}
    snp_dec = {g: (10 * i // n) for i, g in enumerate(by_snp)}
    strata = defaultdict(list)
    for g in genes:
        strata[(len_dec[g], snp_dec[g])].append(g)
    rng = random.Random(seed)
    obs = sum(top[g] for g in clock_present) / len(clock_present)
    ge = 0
    null_means = []
    for _ in range(n_null):
        pick = []
        for g in clock_present:
            pool = strata[(len_dec[g], snp_dec[g])]
            pick.append(rng.choice(pool))
        m = sum(top[g] for g in pick) / len(pick)
        null_means.append(m)
        if m >= obs:
            ge += 1
    null_mean = sum(null_means) / len(null_means)
    return obs, null_mean, (ge + 1) / (n_null + 1)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", default="reference/genes_hg19.tsv")
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank-kb", type=int, default=50)
    ap.add_argument("--n-null", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260702)
    args = ap.parse_args(argv)

    by_chr, length = load_genes(args.genes, args.flank_kb)
    print(f"[enrich2] scanning {args.label} ...", file=sys.stderr)
    top, nsnp = gene_top_chi2(args.outcome, by_chr)
    clock_present = [g for g in CLOCK if g in top]
    obs, null_mean, p = matched_null(top, length, nsnp, clock_present, args.n_null, args.seed)
    ratio = obs / null_mean if null_mean else float("nan")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_header = not out.exists()
    with out.open("a") as fh:
        if write_header:
            fh.write("phenotype\tn_clock_genes\tobs_mean_top_chi2\tmatched_null_mean\tratio\temp_p\n")
        fh.write(f"{args.label}\t{len(clock_present)}\t{obs:.3f}\t{null_mean:.3f}\t{ratio:.3f}\t{p:.2e}\n")
    print(f"[enrich2] {args.label}: obs={obs:.2f} null={null_mean:.2f} ratio={ratio:.2f} "
          f"emp_p={p:.2e} ({len(clock_present)} clock genes)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
