"""Aim 3 (Snakemake `aim3_severity` rule): within-epilepsy pharmacoresistance.

Reports circadian gene-set involvement in the drug-resistant-epilepsy GWAS. NOTE (docs/adr/0006):
DRE is type-confounded; the crude estimate below MUST be type-conditioned (within-type or mtCOJO
on focal/GGE liability) in real runs before any interpretation. This synthetic step emits the
crude value plus an explicit flag that conditioning is pending.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from epicirc.magma.geneset import parse_gsa_out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/analysis_params.yaml")
    ap.add_argument("--magma-dir", default="results/magma")
    ap.add_argument("--out", default="results/aim3_severity.tsv")
    args = ap.parse_args(argv)

    gsa = Path(args.magma_dir) / "circadian_core__drug_resistant_epilepsy.gsa.out"
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("outcome\testimand\tbeta\tse\tp\ttype_conditioned\tnote\n")
        if gsa.exists():
            core = next((r for r in parse_gsa_out(gsa) if r.set_name == "circadian_core"), None)
            if core:
                fh.write(f"drug_resistant_epilepsy\twithin_epilepsy_pharmacoresistance\t"
                         f"{core.beta:.4f}\t{core.beta_se:.4f}\t{core.p:.3e}\tno\t"
                         f"CRUDE — type-conditioning pending (see ADR-0006)\n")
                print(f"[aim3] DRE circadian beta={core.beta:.3f} p={core.p:.2e} (crude)",
                      file=sys.stderr)
                return 0
        fh.write("drug_resistant_epilepsy\twithin_epilepsy_pharmacoresistance\tNA\tNA\tNA\tno\t"
                 "missing_inputs\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
