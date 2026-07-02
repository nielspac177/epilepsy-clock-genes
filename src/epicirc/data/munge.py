"""Standardize/QC summary statistics (LDSC-style munge, simplified for the synthetic path).

Real mode would call LDSC `munge_sumstats.py`; here we validate columns, drop non-finite / out-of
-range rows, and pass through so downstream steps have a clean, uniform table.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

from epicirc.io import COLUMNS, read_sumstats


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    kept = dropped = 0
    with gzip.open(out, "wt", newline="") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in read_sumstats(args.inp):
            try:
                p = float(row["P"]); se = float(row["SE"])
                if not (0.0 < p <= 1.0) or se <= 0:
                    dropped += 1
                    continue
            except (KeyError, ValueError):
                dropped += 1
                continue
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
            kept += 1
    print(f"[munge] {args.inp}: kept {kept}, dropped {dropped} -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
