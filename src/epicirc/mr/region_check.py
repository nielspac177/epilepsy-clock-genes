"""Regional concordance screen — a feasible coloc proxy for the cis-MR hits.

Formal coloc.abf needs full-region eQTL summary stats (eQTLGen's full file is ~16 GB) and SMR-HEIDI
needs an LD panel (unavailable). As an honest intermediate, for a target gene's cis region we ask:
  * where is the GGE association peak (from the FULL ILAE data — unbiased)?
  * where is the eQTL peak (eQTLGen significant SNPs)?
  * is the GGE lead SNP itself a strong eQTL, and vice-versa?
  * how correlated are the two association profiles across shared SNPs?
If the GGE peak coincides with the eQTL peak and the profiles track, colocalization is plausible;
if the GGE signal sits away from the eQTL peak, the cis-MR Wald is likely LD-confounded. This is a
screen, NOT a posterior probability.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

IL = {"CHR": 0, "BP": 1, "SNP": 2, "Z": 8, "P": 9}
EQ = {"SNP": 1, "CHR": 2, "POS": 3, "Z": 6, "SYMBOL": 8}


def load_eqtl_region(eqtl_tsv: str, gene: str) -> dict[int, float]:
    """pos -> |Z| for the gene's significant cis-eQTLs."""
    out: dict[int, float] = {}
    with open(eqtl_tsv) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= EQ["SYMBOL"] or f[EQ["SYMBOL"]] != gene:
                continue
            try:
                out[int(f[EQ["POS"]])] = abs(float(f[EQ["Z"]]))
            except ValueError:
                continue
    return out


def scan_gwas_region(ilae_tbl: str, chrom: str, lo: int, hi: int) -> dict[int, float]:
    """pos -> chi2 for GGE SNPs in the window."""
    out: dict[int, float] = {}
    with open(ilae_tbl) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) < 10 or f[IL["CHR"]] != chrom:
                continue
            try:
                bp = int(f[IL["BP"]])
            except ValueError:
                continue
            if lo <= bp <= hi:
                try:
                    z = float(f[IL["Z"]])
                except ValueError:
                    continue
                out[bp] = z * z
    return out


def _pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    return cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else float("nan")


def check(gene: str, chrom: str, gene_start: int, gene_end: int, eqtl_tsv: str,
          ilae_tbl: str, flank: int = 200_000) -> dict:
    lo, hi = gene_start - flank, gene_end + flank
    eqtl = load_eqtl_region(eqtl_tsv, gene)
    gwas = scan_gwas_region(ilae_tbl, chrom, lo, hi)
    shared = sorted(set(eqtl) & set(gwas))
    res = {"gene": gene, "region": f"chr{chrom}:{lo}-{hi}", "n_eqtl": len(eqtl),
           "n_gwas": len(gwas), "n_shared": len(shared)}
    if gwas:
        gpos = max(gwas, key=gwas.get)
        res["gwas_lead_pos"] = gpos
        res["gwas_lead_chi2"] = round(gwas[gpos], 2)
        res["gwas_lead_is_eqtl"] = gpos in eqtl
    if eqtl:
        epos = max(eqtl, key=eqtl.get)
        res["eqtl_lead_pos"] = epos
        res["eqtl_lead_absZ"] = round(eqtl[epos], 2)
        res["eqtl_lead_gwas_chi2"] = round(gwas.get(epos, float("nan")), 2)
    if "gwas_lead_pos" in res and "eqtl_lead_pos" in res:
        res["lead_distance_kb"] = round(abs(res["gwas_lead_pos"] - res["eqtl_lead_pos"]) / 1000, 1)
    if len(shared) >= 3:
        res["profile_corr"] = round(_pearson([gwas[p] for p in shared],
                                             [eqtl[p] for p in shared]), 3)
    return res


GENES = {  # gene: (chrom, start_hg19, end_hg19)
    "ARNTL": ("11", 13298199, 13408813),
    "NR1D1": ("17", 38249040, 38256978),
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eqtl", default="data/raw/eqtl/clock_cis_eqtls.tsv")
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", nargs="*", default=["ARNTL", "NR1D1"])
    args = ap.parse_args(argv)

    rows = [check(g, *GENES[g], args.eqtl, args.outcome) for g in args.genes]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = ["gene", "region", "n_eqtl", "n_gwas", "n_shared", "gwas_lead_pos", "gwas_lead_chi2",
            "gwas_lead_is_eqtl", "eqtl_lead_pos", "eqtl_lead_absZ", "eqtl_lead_gwas_chi2",
            "lead_distance_kb", "profile_corr"]
    with out.open("w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "NA")) for k in keys) + "\n")
    for r in rows:
        print(f"[region] {r['gene']}: GWAS lead chi2={r.get('gwas_lead_chi2')} "
              f"(is_eqtl={r.get('gwas_lead_is_eqtl')}), eQTL lead GGE chi2={r.get('eqtl_lead_gwas_chi2')}, "
              f"lead dist={r.get('lead_distance_kb')}kb, profile r={r.get('profile_corr')}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
