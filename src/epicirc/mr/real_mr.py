"""Real two-sample MR: circadian/sleep exposure -> ILAE epilepsy outcome.

Implements the playbook flow with the tools we actually have (no PLINK/LD panel):
  1. read exposure sumstats, keep genome-wide-significant SNPs (P < --pval);
  2. distance-clump to ~independent instruments (epicirc.mr.clump);
  3. match to the ILAE .tbl outcome by CHR:POS (NOT rsID) — one streamed pass over the big file;
  4. harmonize outcome effect to the exposure effect allele, dropping palindromes/incompatibles;
  5. IVW + MR-Egger (epicirc.mr.two_sample).

BUILD SAFETY: the ILAE files are GRCh37/hg19. The exposure MUST be GRCh37 too (a GRCh38 harmonised
file will match ~0 SNPs). Pass GRCh37 columns via --exp-* flags. Sample overlap is low here (ILAE
has no UK Biobank; UKB-based sleep exposures overlap only the controls minimally) so estimates are
conservative — but still report it. Always run a positive-control exposure through this same driver
before trusting a null.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

from epicirc.data.harmonize import harmonize
from epicirc.mr.clump import Variant, distance_clump
from epicirc.mr.two_sample import ivw, mr_egger

# ILAE .tbl column indices (whitespace-delimited)
IL = {"CHR": 0, "BP": 1, "SNP": 2, "A1": 3, "A2": 4, "BETA": 11, "SE": 12}


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_exposure(path: str, cols: dict, pval: float) -> dict[tuple[str, int], dict]:
    """Return {(chr,pos): {ea,oa,beta,p}} for exposure SNPs passing P<pval."""
    out: dict[tuple[str, int], dict] = {}
    with _open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                p = float(row[cols["p"]])
                if p >= pval:
                    continue
                chrom = str(row[cols["chr"]]).replace("chr", "")
                pos = int(float(row[cols["pos"]]))
                out[(chrom, pos)] = {
                    "ea": row[cols["ea"]].upper(), "oa": row[cols["oa"]].upper(),
                    "beta": float(row[cols["beta"]]), "p": p}
            except (KeyError, ValueError):
                continue
    return out


def match_outcome(ilae_tbl: str, keys: set[tuple[str, int]]) -> dict[tuple[str, int], dict]:
    """One streamed pass over the ILAE .tbl, pulling only the instrument CHR:POS."""
    found: dict[tuple[str, int], dict] = {}
    with open(ilae_tbl) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) < 13:
                continue
            key = (f[IL["CHR"]], int(f[IL["BP"]]))
            if key in keys:
                try:
                    found[key] = {"a1": f[IL["A1"]].upper(), "a2": f[IL["A2"]].upper(),
                                  "beta": float(f[IL["BETA"]]), "se": float(f[IL["SE"]])}
                except ValueError:
                    continue
    return found


def run(exposure: str, cols: dict, ilae_tbl: str, pval: float, clump_kb: int):
    exp = read_exposure(exposure, cols, pval)
    variants = [Variant(f"{c}:{p}", c, p, exp[(c, p)]["p"]) for (c, p) in exp]
    clumped = distance_clump(variants, window_kb=clump_kb)
    keys = {(v.chrom, v.pos) for v in clumped}
    out = match_outcome(ilae_tbl, keys)

    bx, by, sey, used = [], [], [], 0
    for (c, p), o in out.items():
        e = exp[(c, p)]
        h = harmonize(o["beta"], o["a1"], o["a2"], e["ea"], e["oa"])
        if not h.usable:
            continue
        bx.append(e["beta"]); by.append(h.beta); sey.append(o["se"]); used += 1
    return bx, by, sey, len(clumped), used


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exposure", required=True, help="exposure sumstats (GRCh37!)")
    ap.add_argument("--outcome", required=True, help="ILAE .tbl outcome file")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pval", type=float, default=5e-8)
    ap.add_argument("--clump-kb", type=int, default=1000)
    # exposure column names (defaults for a typical sleep-trait file; override per dataset)
    ap.add_argument("--exp-chr", default="CHR")
    ap.add_argument("--exp-pos", default="POS")
    ap.add_argument("--exp-ea", default="Effect_Allele")
    ap.add_argument("--exp-oa", default="Other_Allele")
    ap.add_argument("--exp-beta", default="Effect")
    ap.add_argument("--exp-p", default="Pval")
    args = ap.parse_args(argv)

    cols = {"chr": args.exp_chr, "pos": args.exp_pos, "ea": args.exp_ea, "oa": args.exp_oa,
            "beta": args.exp_beta, "p": args.exp_p}
    bx, by, sey, n_clumped, n_used = run(args.exposure, cols, args.outcome, args.pval, args.clump_kb)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("method\tn_instruments\tbeta\tse\tp\textra\n")
        if n_used < 3:
            fh.write(f"insufficient_instruments\t{n_used}\tNA\tNA\tNA\t"
                     f"{{'clumped':{n_clumped}}}\n")
            print(f"[real-mr] {n_used}/{n_clumped} instruments matched — too few", file=sys.stderr)
            return 0
        for est in (ivw(bx, by, sey), mr_egger(bx, by, sey)):
            fh.write(f"{est.method}\t{est.n_snps}\t{est.beta:.5f}\t{est.se:.5f}\t"
                     f"{est.p:.3e}\t{est.extra}\n")
    print(f"[real-mr] {n_used}/{n_clumped} instruments; IVW+Egger -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
