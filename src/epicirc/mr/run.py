"""Two-sample MR (Snakemake `mr` rule): exposure -> outcome.

Selects instruments from the exposure (genome-wide or top-k fallback), harmonizes alleles to the
outcome, drops palindromic/incompatible SNPs, and runs the tested IVW + MR-Egger estimators.
Full sensitivity (weighted median, MR-PRESSO, CAUSE, MRlap) is delegated to R in real runs
(epicirc.mr.two_sample.run_twosamplemr_r).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from epicirc.data.harmonize import harmonize
from epicirc.io import read_sumstats
from epicirc.mr.two_sample import ivw, mr_egger


def select_and_harmonize(exposure: str, outcome: str, pval: float = 5e-8, top_k: int = 40):
    exp = {r["SNP"]: r for r in read_sumstats(exposure)}
    out = {r["SNP"]: r for r in read_sumstats(outcome)}
    shared = [s for s in exp if s in out]
    ranked = sorted(shared, key=lambda s: float(exp[s]["P"]))
    hits = [s for s in ranked if float(exp[s]["P"]) < pval] or ranked[:top_k]

    bx, by, sey = [], [], []
    for s in hits:
        e, o = exp[s], out[s]
        h = harmonize(float(o["BETA"]), o["A1"], o["A2"], e["A1"], e["A2"])
        if not h.usable:
            continue
        bx.append(float(e["BETA"]))
        by.append(h.beta)
        sey.append(float(o["SE"]))
    return bx, by, sey


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exposure", required=True)
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    bx, by, sey = select_and_harmonize(args.exposure, args.outcome)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("method\tn_snps\tbeta\tse\tp\textra\n")
        if len(bx) < 3:
            fh.write(f"insufficient_instruments\t{len(bx)}\tNA\tNA\tNA\t{{}}\n")
            print(f"[mr] only {len(bx)} instruments; skipping -> {out}", file=sys.stderr)
            return 0
        for est in (ivw(bx, by, sey), mr_egger(bx, by, sey)):
            fh.write(f"{est.method}\t{est.n_snps}\t{est.beta:.5f}\t{est.se:.5f}\t"
                     f"{est.p:.3e}\t{est.extra}\n")
    print(f"[mr] {len(bx)} instruments; IVW+Egger -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
