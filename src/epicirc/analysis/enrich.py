"""Tool-free circadian-locus enrichment on real ILAE summary statistics.

A FIRST LOOK, not a substitute for MAGMA/stratified-LDSC (docs/adr/0005). For each phenotype we
compare the mean chi-square of SNPs in/near the 23 core clock genes against the genome-wide mean.
We report the enrichment as a RATIO (mean_chi2_circadian / mean_chi2_genome), which normalizes for
each phenotype's overall power/inflation, plus a z-test. Because SNPs within a gene are in LD, the
all-SNP z is anti-conservative; we therefore also report an LD-robust top-SNP-per-gene summary.

Caveats (stated in any write-up): no LD pruning within loci for the all-SNP test; genome mean
includes other trait loci; this cannot separate circadian biology from co-located genes. It exists
to see whether the real data point the way the hypothesis predicts before the full pipeline runs.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

# real ILAE .tbl columns (whitespace-delimited)
COL = {"CHR": 0, "BP": 1, "SNP": 2, "A1": 3, "A2": 4, "FRQ": 5, "N": 7, "Z": 8, "P": 9}


def load_windows(coord_tsv: str, window_kb: int) -> dict[str, list[tuple[int, int, str]]]:
    w = window_kb * 1000
    out: dict[str, list[tuple[int, int, str]]] = {}
    with open(coord_tsv) as fh:
        next(fh)
        for ln in fh:
            sym, chrom, start, end = ln.split("\t")
            out.setdefault(chrom.strip(), []).append(
                (int(start) - w, int(end) + w, sym))
    return out


def _in_window(chrom: str, bp: int, windows) -> str | None:
    for lo, hi, sym in windows.get(chrom, ()):  # ~<=few genes per chr
        if lo <= bp <= hi:
            return sym
    return None


def enrich_file(path: str, windows, reservoir_stride: int = 37) -> dict:
    g_sum = g_sumsq = 0.0
    g_n = 0
    circ_chi2: list[float] = []
    top_by_gene: dict[str, float] = {}
    top_detail: dict[str, tuple[float, str, str]] = {}  # gene -> (chi2, P, SNP)
    reservoir: list[float] = []  # strided sample for median (lambda_gc)

    with open(path) as fh:
        next(fh)  # header
        for i, ln in enumerate(fh):
            f = ln.split()
            if len(f) < 10:
                continue
            zt = f[COL["Z"]]
            if zt in ("NA", "nan", ""):
                continue
            try:
                z = float(zt)
            except ValueError:
                continue
            chi2 = z * z
            g_sum += chi2
            g_sumsq += chi2 * chi2
            g_n += 1
            if i % reservoir_stride == 0:
                reservoir.append(chi2)
            sym = _in_window(f[COL["CHR"]], int(f[COL["BP"]]), windows)
            if sym is not None:
                circ_chi2.append(chi2)
                if chi2 > top_by_gene.get(sym, -1.0):
                    top_by_gene[sym] = chi2
                    top_detail[sym] = (chi2, f[COL["P"]], f[COL["SNP"]])

    if g_n == 0 or not circ_chi2:
        raise ValueError(f"no usable rows / no circadian SNPs in {path}")

    g_mean = g_sum / g_n
    g_var = max(g_sumsq / g_n - g_mean * g_mean, 1e-12)
    g_sd = math.sqrt(g_var)
    n_c = len(circ_chi2)
    c_mean = sum(circ_chi2) / n_c
    z_enrich = (c_mean - g_mean) / (g_sd / math.sqrt(n_c))
    p_enrich = 0.5 * math.erfc(z_enrich / math.sqrt(2.0))  # one-sided (enrichment)
    reservoir.sort()
    med = reservoir[len(reservoir) // 2] if reservoir else float("nan")
    lambda_gc = med / 0.4549
    top_vals = sorted(top_by_gene.values(), reverse=True)
    top_mean = sum(top_vals) / len(top_vals)

    return {
        "n_genome_snps": g_n, "n_circadian_snps": n_c,
        "genome_mean_chi2": g_mean, "circadian_mean_chi2": c_mean,
        "enrichment_ratio": c_mean / g_mean, "z_enrich": z_enrich, "p_enrich": p_enrich,
        "lambda_gc": lambda_gc, "n_genes_hit": len(top_by_gene),
        "top_snp_per_gene_mean_chi2": top_mean,
        "per_gene": top_detail,
    }


PHENO_FILE = {
    "all_epilepsy": "ILAE3_Caucasian_all_epilepsy_final.tbl",
    "focal_epilepsy": "ILAE3_Caucasian_focal_epilepsy_final.tbl",
    "generalized_gge": "ILAE3_Caucasian_GGE_final.tbl",
    "jme": "ILAE3_JME_final.tbl", "jae": "ILAE3_JAE_final.tbl", "cae": "ILAE3_CAE_final.tbl",
    "gtcs_alone": "ILAE3_GTCS_final.tbl", "focal_hs": "ILAE3_focal_HS_final.tbl",
    "focal_lesion_neg": "ILAE3_focal_lesion_negative_final.tbl",
    "focal_other_lesion": "ILAE3_focal_other_lesion_final.tbl",
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", default="data/raw/final_sumstats")
    ap.add_argument("--coords", default="config/gene_sets/circadian_core_hg19.tsv")
    ap.add_argument("--window-kb", type=int, default=50)
    ap.add_argument("--phenotypes", nargs="*", default=["all_epilepsy", "focal_epilepsy", "generalized_gge"])
    ap.add_argument("--out", default="results/real_circadian_enrichment.tsv")
    args = ap.parse_args(argv)

    windows = load_windows(args.coords, args.window_kb)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["phenotype", "n_circadian_snps", "n_genes_hit", "genome_mean_chi2",
            "circadian_mean_chi2", "enrichment_ratio", "z_enrich", "p_enrich",
            "top_snp_per_gene_mean_chi2", "lambda_gc", "n_genome_snps"]
    rows = []
    for pheno in args.phenotypes:
        path = Path(args.raw_dir) / PHENO_FILE[pheno]
        if not path.exists():
            print(f"[skip] {pheno}: {path} not found", file=sys.stderr)
            continue
        print(f"[enrich] scanning {pheno} ...", file=sys.stderr)
        r = enrich_file(str(path), windows)
        r["phenotype"] = pheno
        # per-gene table (broad-vs-single-locus audit)
        pg = out.with_name(f"real_circadian_pergene_{pheno}.tsv")
        with pg.open("w") as gh:
            gh.write("gene\ttop_chi2\ttop_P\ttop_SNP\tgws\n")
            for sym, (chi2, p, snp) in sorted(r.pop("per_gene").items(), key=lambda x: -x[1][0]):
                gh.write(f"{sym}\t{chi2:.2f}\t{p}\t{snp}\t{int(float(p) < 5e-8)}\n")
        rows.append(r)
        print(f"  {pheno}: ratio={r['enrichment_ratio']:.3f} z={r['z_enrich']:.2f} "
              f"p={r['p_enrich']:.2e} lambda={r['lambda_gc']:.3f} "
              f"({r['n_circadian_snps']} SNPs, {r['n_genes_hit']} genes)", file=sys.stderr)

    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"[enrich] -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
