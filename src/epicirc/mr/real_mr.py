"""Real two-sample MR: circadian/sleep exposure -> ILAE epilepsy outcome.

Playbook flow with the tools we have (no PLINK/LD panel):
  1. read exposure sumstats, keep P < --pval;
  2. distance-clump to ~independent instruments (epicirc.mr.clump);
  3. match to the ILAE .tbl outcome by CHR:POS (NOT rsID) in one streamed pass;
  4. harmonize outcome effect to the exposure effect allele, dropping palindromes/incompatibles;
  5. run IVW + MR-Egger + weighted median + MR-PRESSO-style + Steiger (epicirc.mr.robust).

BUILD SAFETY: ILAE is GRCh37/hg19; the exposure MUST be GRCh37 (a GRCh38 file matches ~0 SNPs).
Sample overlap is low (ILAE has no UK Biobank) -> conservative. Emits a per-instrument table so
sensitivity analyses can be re-run without touching the big files.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

from epicirc.data.harmonize import harmonize
from epicirc.mr.clump import Variant, distance_clump
from epicirc.mr.robust import mr_presso_lite, steiger_filter, weighted_median
from epicirc.mr.two_sample import ivw, mr_egger

# ILAE .tbl column indices (whitespace-delimited)
IL = {"CHR": 0, "BP": 1, "SNP": 2, "A1": 3, "A2": 4, "N": 7, "BETA": 11, "SE": 12}


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_exposure(path: str, cols: dict, pval: float) -> dict[tuple[str, int], dict]:
    out: dict[tuple[str, int], dict] = {}
    with _open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                p = float(row[cols["p"]])
                if p >= pval:
                    continue
                chrom = str(row[cols["chr"]]).replace("chr", "")
                pos = int(float(row[cols["pos"]]))
                out[(chrom, pos)] = {
                    "ea": row[cols["ea"]].upper(), "oa": row[cols["oa"]].upper(),
                    "beta": float(row[cols["beta"]]), "se": float(row[cols["se"]]),
                    "eaf": float(row[cols["eaf"]]), "p": p}
            except (KeyError, ValueError):
                continue
    return out


def match_outcome(ilae_tbl: str, keys: set[tuple[str, int]]) -> dict[tuple[str, int], dict]:
    found: dict[tuple[str, int], dict] = {}
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
                    found[key] = {"a1": f[IL["A1"]].upper(), "a2": f[IL["A2"]].upper(),
                                  "beta": float(f[IL["BETA"]]), "se": float(f[IL["SE"]]),
                                  "n": float(f[IL["N"]])}
                except ValueError:
                    continue
    return found


def build_instruments(exposure: str, cols: dict, ilae_tbl: str, pval: float, clump_kb: int):
    exp = read_exposure(exposure, cols, pval)
    variants = [Variant(f"{c}:{p}", c, p, exp[(c, p)]["p"]) for (c, p) in exp]
    clumped = distance_clump(variants, window_kb=clump_kb)
    keys = {(v.chrom, v.pos) for v in clumped}
    out = match_outcome(ilae_tbl, keys)
    rows = []
    for (c, p), o in out.items():
        e = exp[(c, p)]
        h = harmonize(o["beta"], o["a1"], o["a2"], e["ea"], e["oa"])
        if not h.usable:
            continue
        rows.append({"snp": f"{c}:{p}", "chr": c, "pos": p, "beta_exp": e["beta"],
                     "se_exp": e["se"], "eaf": e["eaf"], "beta_out": h.beta,
                     "se_out": o["se"], "n_out": o["n"]})
    return rows, len(clumped)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exposure", required=True)
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pval", type=float, default=5e-8)
    ap.add_argument("--clump-kb", type=int, default=1000)
    ap.add_argument("--exp-chr", default="CHR")
    ap.add_argument("--exp-pos", default="BP")
    ap.add_argument("--exp-ea", default="ALLELE1")
    ap.add_argument("--exp-oa", default="ALLELE0")
    ap.add_argument("--exp-beta", default="LOGOR")
    ap.add_argument("--exp-se", default="LOGOR_SE")
    ap.add_argument("--exp-eaf", default="A1FREQ")
    ap.add_argument("--exp-p", default="P_BOLT_LMM")
    ap.add_argument("--n-exp", type=float, default=400000.0)
    args = ap.parse_args(argv)

    cols = {"chr": args.exp_chr, "pos": args.exp_pos, "ea": args.exp_ea, "oa": args.exp_oa,
            "beta": args.exp_beta, "se": args.exp_se, "eaf": args.exp_eaf, "p": args.exp_p}
    rows, n_clumped = build_instruments(args.exposure, cols, args.outcome, args.pval, args.clump_kb)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    # per-instrument table (for re-running sensitivity without the big files)
    inst = out.with_name(out.stem + ".instruments.tsv")
    with inst.open("w") as fh:
        fh.write("snp\tchr\tpos\tbeta_exp\tse_exp\teaf\tbeta_out\tse_out\tn_out\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in
                     ["snp", "chr", "pos", "beta_exp", "se_exp", "eaf", "beta_out", "se_out", "n_out"]) + "\n")

    with out.open("w") as fh:
        fh.write("method\tn_snps\tbeta\tse\tp\textra\n")
        if len(rows) < 3:
            fh.write(f"insufficient_instruments\t{len(rows)}\tNA\tNA\tNA\t{{'clumped':{n_clumped}}}\n")
            print(f"[real-mr] {len(rows)}/{n_clumped} instruments — too few", file=sys.stderr)
            return 0
        bx = [r["beta_exp"] for r in rows]; sx = [r["se_exp"] for r in rows]
        by = [r["beta_out"] for r in rows]; sy = [r["se_out"] for r in rows]
        eaf = [r["eaf"] for r in rows]; nout = rows[0]["n_out"]

        e_ivw = ivw(bx, by, sy)
        e_egg = mr_egger(bx, by, sy)
        e_wm = weighted_median(bx, by, sy)
        e_pr = mr_presso_lite(bx, by, sy)
        e_st = steiger_filter(bx, sx, eaf, args.n_exp, by, sy, nout)
        for est in (e_ivw, e_egg):
            fh.write(f"{est.method}\t{est.n_snps}\t{est.beta:.5f}\t{est.se:.5f}\t{est.p:.3e}\t{est.extra}\n")
        for est in (e_wm, e_pr, e_st):
            fh.write(f"{est.method}\t{est.n_snps}\t{est.beta:.5f}\t{est.se:.5f}\t{est.p:.3e}\t{est.extra}\n")
    print(f"[real-mr] {len(rows)}/{n_clumped} instruments; 5 methods -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
