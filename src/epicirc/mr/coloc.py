"""Formal colocalization (coloc.abf, Giambartolomei 2014 / Wakefield ABF) for clock genes vs GGE.

Tests whether a clock gene's brain cis-eQTL signal (PsychENCODE prefrontal cortex, hg19, N=1387)
and the epilepsy GWAS signal share a single causal variant in the cis region.

Per SNP (matched CHR:POS): eQTL contributes a Wakefield ABF from its p-value + MAF + N (quantitative,
sdY=1); the GWAS contributes an ABF from ILAE beta/se directly. coloc is sign-agnostic (uses Z^2), so
PsychENCODE's missing allele column is not a problem. Posterior PP.H4 near 1 => shared causal variant
(colocalization); high PP.H3 => distinct variants (the cis-MR Wald would then be LD-confounded).
Priors p1=p2=1e-4, p12=1e-5 (coloc defaults). MAF from ILAE A1FREQ.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from epicirc.mr.robust import _probit

N_EQTL = 1387
W_EQTL = 0.15 ** 2      # prior variance, quantitative eQTL
W_GWAS = 0.2 ** 2       # prior variance, case/control log-odds
P1 = P2 = 1e-4
P12 = 1e-5

# PsychENCODE clock file (whitespace): gene_id0 ... SNP_chr8 SNP_start9 ... nominal_pval11 slope12
PE = {"GENE": 0, "CHR": 8, "POS": 9, "P": 11}
IL = {"CHR": 0, "BP": 1, "FRQ": 5, "BETA": 11, "SE": 12}

GENE_ENSG = {"ARNTL": "ENSG00000133794", "PER1": "ENSG00000179094", "NR1D1": "ENSG00000126368",
             "CRY2": "ENSG00000121671", "RORA": "ENSG00000069667", "NPAS2": "ENSG00000170485"}


def _logsumexp(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def _z2_from_p(p: float) -> float:
    p = min(max(p, 1e-300), 1.0)
    arg = min(1.0 - p / 2.0, 1.0 - 1e-15)   # avoid _probit(1.0) domain error for tiny p
    z = _probit(arg)
    return z * z


def _labf(z2: float, v: float, w: float) -> float:
    r = w / (v + w)
    return 0.5 * (math.log(1 - r) + r * z2)


def load_eqtl_gene(pe_file: str, ensg: str) -> dict[tuple[str, int], float]:
    """(chr,pos) -> min nominal p for the gene's cis-eQTLs."""
    out: dict[tuple[str, int], float] = {}
    with open(pe_file) as fh:
        for ln in fh:
            f = ln.split()
            if len(f) <= PE["P"] or not f[PE["GENE"]].startswith(ensg):
                continue
            try:
                key = (f[PE["CHR"]].replace("chr", ""), int(f[PE["POS"]]))
                p = float(f[PE["P"]])
            except ValueError:
                continue
            if key not in out or p < out[key]:
                out[key] = p
    return out


def scan_gwas(ilae_tbl: str, keys: set[tuple[str, int]]) -> dict[tuple[str, int], tuple]:
    found = {}
    with open(ilae_tbl) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) < 13:
                continue
            try:
                key = (f[IL["CHR"]], int(f[IL["BP"]]))
            except ValueError:
                continue
            if key in keys:
                try:
                    found[key] = (float(f[IL["BETA"]]), float(f[IL["SE"]]), float(f[IL["FRQ"]]))
                except ValueError:
                    continue
    return found


def coloc(eqtl: dict, gwas: dict) -> dict:
    shared = sorted(set(eqtl) & set(gwas))
    labf_e, labf_g = [], []
    for k in shared:
        beta_g, se_g, frq = gwas[k]
        maf = min(frq, 1 - frq)
        if maf <= 0 or se_g <= 0:
            continue
        v_e = 1.0 / (2 * N_EQTL * maf * (1 - maf))     # sdY=1 quantitative approx
        labf_e.append(_labf(_z2_from_p(eqtl[k]), v_e, W_EQTL))
        labf_g.append(_labf((beta_g / se_g) ** 2, se_g ** 2, W_GWAS))
    n = len(labf_e)
    if n < 5:
        return {"n_snps": n, "PP.H4": float("nan")}
    l1 = _logsumexp(labf_e)
    l2 = _logsumexp(labf_g)
    l4 = _logsumexp([a + b for a, b in zip(labf_e, labf_g)])
    # H3 = sum_{i!=j} = exp(l1+l2) - exp(l4), in log space
    l3 = math.log(max(math.exp(l1 + l2) - math.exp(l4), 1e-300))
    terms = [0.0, math.log(P1) + l1, math.log(P2) + l2,
             math.log(P1) + math.log(P2) + l3, math.log(P12) + l4]
    denom = _logsumexp(terms)
    pp = [math.exp(t - denom) for t in terms]
    return {"n_snps": n, "PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2],
            "PP.H3": pp[3], "PP.H4": pp[4]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eqtl", default="data/raw/eqtl/psychencode_clock_cis.tsv")
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", nargs="*", default=["ARNTL", "PER1", "NR1D1", "CRY2", "RORA", "NPAS2"])
    args = ap.parse_args(argv)

    rows = []
    for gene in args.genes:
        ensg = GENE_ENSG[gene]
        eqtl = load_eqtl_gene(args.eqtl, ensg)
        gwas = scan_gwas(args.outcome, set(eqtl))
        r = coloc(eqtl, gwas)
        r["gene"] = gene
        rows.append(r)
        h4 = r.get("PP.H4")
        print(f"[coloc] {gene}: PP.H4={h4:.3f} PP.H3={r.get('PP.H3', float('nan')):.3f} "
              f"(n={r['n_snps']})" if h4 == h4 else f"[coloc] {gene}: n={r['n_snps']} (too few)",
              file=sys.stderr)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = ["gene", "n_snps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    with out.open("w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r.get(k), float) else str(r.get(k, "NA"))
                               for k in keys) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
