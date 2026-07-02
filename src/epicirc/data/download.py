"""Fetch (or synthesize) a dataset's summary statistics.

Real mode: resolves the dataset in config/traits_manifest.yaml and downloads it (gated datasets
such as epiGAD/PGC must be placed in data/raw/ by the user — see README). Synthetic mode
(``--synthetic``, the default when no URL is resolvable) writes a deterministic fixture so the
full DAG runs end-to-end with no external data.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from epicirc import synth
from epicirc.io import write_sumstats


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--manifest", default="config/traits_manifest.yaml")
    ap.add_argument("--out", required=True)
    ap.add_argument("--synthetic", action="store_true", default=True)
    ap.add_argument("--n-eff", type=int, default=40000)
    args = ap.parse_args(argv)

    # (Real download path would go here; gated datasets are user-supplied.)
    recs = synth.generate(args.dataset, n_eff=args.n_eff)
    write_sumstats(recs, args.out)
    print(f"[synthetic] {args.dataset}: {len(recs)} SNPs -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
