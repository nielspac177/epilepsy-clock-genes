"""Genetic correlation (Snakemake `ldsc_rg` rule).

Real mode calls LDSC (epicirc.ldsc.rg.run_rg) in the legacy Python-2 env. Synthetic mode estimates
rg as the Pearson correlation of per-SNP effects over shared SNPs (Fisher-z p-value) and emits a
log in LDSC's format so epicirc.ldsc.rg.parse_rg_log can read it unchanged. The cross-trait
intercept (gcov_int) is emitted as the shared-control rho proxy for the MAGMA heterogeneity arm.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from epicirc.io import read_sumstats


def _fisher_p(r: float, n: int) -> float:
    if n <= 3 or abs(r) >= 1:
        return 1.0
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(n - 3)
    return 2.0 * (0.5 * math.erfc(abs(z) / math.sqrt(2.0)))


def _load_betas(path: str) -> dict[str, float]:
    return {row["SNP"]: float(row["BETA"]) for row in read_sumstats(path)}


def _h2_block(betas: dict[str, float], which: int) -> tuple[str, float, float]:
    zs = list(betas.values())
    n = len(zs)
    mean_chi2 = sum(z * z for z in zs) / n if n else 1.0
    h2 = max(mean_chi2 - 1.0, 1e-4)
    se = max(mean_chi2 * math.sqrt(2.0 / max(n, 2)), 1e-6)
    block = (f"Heritability of phenotype {which}\n"
             f"Total Observed scale h2: {h2:.4f} ({se:.4f})\n")
    return block, h2, se


def rg_from_files(a: str, b: str) -> dict:
    ba, bb = _load_betas(a), _load_betas(b)
    shared = sorted(set(ba) & set(bb))
    n = len(shared)
    xa = [ba[s] for s in shared]
    xb = [bb[s] for s in shared]
    ma, mb = sum(xa) / n, sum(xb) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(xa, xb)) / n
    va = sum((x - ma) ** 2 for x in xa) / n
    vb = sum((y - mb) ** 2 for y in xb) / n
    r = cov / math.sqrt(va * vb) if va > 0 and vb > 0 else 0.0
    return {"rg": r, "se": (1 - r * r) / math.sqrt(max(n - 3, 1)), "p": _fisher_p(r, n),
            "gcov_int": round(cov * 5, 4), "gcov_int_se": 0.005, "n": n}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--out", required=True)   # basename; writes {out}.log
    args = ap.parse_args(argv)

    res = rg_from_files(args.a, args.b)
    ba, bb = _load_betas(args.a), _load_betas(args.b)
    blk1, *_ = _h2_block(ba, 1)
    blk2, *_ = _h2_block(bb, 2)

    log = Path(args.out).with_suffix(".log")
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as fh:
        fh.write(blk1 + "\n" + blk2 + "\n")
        fh.write("Summary of Genetic Correlation Results\n")
        fh.write("p1 p2 rg se z p gcov_int gcov_int_se\n")
        fh.write(f"{Path(args.a).name} {Path(args.b).name} "
                 f"{res['rg']:.4f} {res['se']:.4f} {res['rg']/max(res['se'],1e-9):.4f} "
                 f"{res['p']:.3e} {res['gcov_int']:.4f} {res['gcov_int_se']:.4f}\n")
    print(f"[ldsc-rg] rg={res['rg']:.3f} p={res['p']:.2e} (n={res['n']}) -> {log}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
