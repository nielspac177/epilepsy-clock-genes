"""Gene-level aggregation (Snakemake `magma_gene` rule).

Real mode calls MAGMA (see epicirc.magma.geneset.gene_analysis) once a macOS binary + reference
panel are installed. Synthetic mode aggregates per-SNP z to a gene z (Stouffer, LD ignored) so the
gene-set step has input. Output columns: GENE CHR NSNPS ZSTAT CIRCADIAN.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

from epicirc.io import read_sumstats


def aggregate_genes(sumstats_path: str) -> list[tuple[str, int, int, float, int]]:
    z_by_gene: dict[str, list[float]] = defaultdict(list)
    circ_by_gene: dict[str, int] = {}
    chr_by_gene: dict[str, int] = {}
    for row in read_sumstats(sumstats_path):
        gene = row["GENE"]
        z = float(row["BETA"]) / float(row["SE"])
        z_by_gene[gene].append(z)
        circ_by_gene[gene] = int(row["CIRCADIAN"])
        chr_by_gene[gene] = int(row["CHR"])
    out = []
    for gene, zs in z_by_gene.items():
        m = len(zs)
        gene_z = sum(zs) / math.sqrt(m) if m else 0.0
        out.append((gene, chr_by_gene[gene], m, gene_z, circ_by_gene[gene]))
    return sorted(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sumstats", required=True)
    ap.add_argument("--epi", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    genes = aggregate_genes(args.sumstats)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("GENE\tCHR\tNSNPS\tZSTAT\tCIRCADIAN\n")
        for gene, chrom, m, z, circ in genes:
            fh.write(f"{gene}\t{chrom}\t{m}\t{z:.6f}\t{circ}\n")
    print(f"[magma-gene] {len(genes)} genes -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
