"""Competitive gene-set analysis (Snakemake `magma_geneset` rule).

Competitive test: regress gene ZSTAT on set-membership. The beta is the standardized mean
difference (genes in set vs out); p is one-sided (enrichment). Output mimics MAGMA .gsa.out:
VARIABLE NGENES BETA BETA_STD SE P. Set membership resolved from the .genes.raw CIRCADIAN column
(for circadian_core) or a gene-list file / config set name.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from epicirc.geneset.curate import load_geneset


def _norm_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def _read_genes(genes_raw: str) -> list[tuple[str, float, int]]:
    rows = []
    with open(genes_raw) as fh:
        next(fh)  # header
        for ln in fh:
            gene, _chr, _n, z, circ = ln.split("\t")
            rows.append((gene, float(z), int(circ)))
    return rows


def _membership(genes, set_name: str, set_file: str | None) -> set[str]:
    if set_file and Path(set_file).exists() and set_file.endswith(".tsv"):
        return set(load_geneset(set_file).symbols)
    if set_name == "circadian_core":
        return {g for g, _z, circ in genes if circ == 1}
    if set_name == "negative_controls":
        # a size-matched background slice (deterministic): first N non-circadian genes
        n = sum(1 for _g, _z, c in genes if c == 1)
        bg = [g for g, _z, c in genes if c == 0]
        return set(bg[:n])
    cfg = Path("config/gene_sets") / f"{set_name}.tsv"
    if cfg.exists():
        return set(load_geneset(cfg).symbols)
    return set()


def competitive_test(genes, in_set: set[str]) -> tuple[int, float, float, float]:
    z_in = [z for g, z, _c in genes if g in in_set]
    z_out = [z for g, z, _c in genes if g not in in_set]
    n_in, n_out = len(z_in), len(z_out)
    if n_in < 2 or n_out < 2:
        raise ValueError("need >=2 genes in and out of set")
    m_in, m_out = sum(z_in) / n_in, sum(z_out) / n_out
    var_in = sum((z - m_in) ** 2 for z in z_in) / (n_in - 1)
    var_out = sum((z - m_out) ** 2 for z in z_out) / (n_out - 1)
    se = math.sqrt(var_in / n_in + var_out / n_out)
    beta = m_in - m_out
    pooled_sd = math.sqrt((var_in * (n_in - 1) + var_out * (n_out - 1)) / (n_in + n_out - 2))
    beta_std = beta / pooled_sd if pooled_sd > 0 else 0.0
    z = beta / se if se > 0 else 0.0
    p_one_sided = _norm_sf(z)          # enrichment (positive beta) is the alternative
    return n_in, beta, beta_std, se, p_one_sided


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--set", dest="set_arg", required=True, help="set name or .tsv gene list")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    set_name = Path(args.set_arg).stem
    genes = _read_genes(args.genes)
    in_set = _membership(genes, set_name, args.set_arg if args.set_arg.endswith(".tsv") else None)
    n_in, beta, beta_std, se, p = competitive_test(genes, in_set)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("VARIABLE\tTYPE\tNGENES\tBETA\tBETA_STD\tSE\tP\n")
        fh.write(f"{set_name}\tSET\t{n_in}\t{beta:.6f}\t{beta_std:.6f}\t{se:.6f}\t{p:.3e}\n")
    print(f"[magma-set] {set_name}: NGENES={n_in} BETA={beta:.3f} P={p:.2e} -> {out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
