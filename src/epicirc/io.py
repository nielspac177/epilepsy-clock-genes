"""Minimal gzipped-TSV read/write for summary statistics (stdlib only)."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Iterable, Iterator

COLUMNS = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N", "GENE", "CIRCADIAN"]


def write_sumstats(records: Iterable, out: str | Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUMNS)
        for r in records:
            w.writerow([r.snp, r.chrom, r.bp, r.a1, r.a2, r.frq, r.beta, r.se, r.p, r.n,
                        r.gene, int(r.circadian)])
    return out


def read_sumstats(path: str | Path) -> Iterator[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row
