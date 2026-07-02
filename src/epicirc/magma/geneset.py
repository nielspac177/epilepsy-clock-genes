"""Wrapper around MAGMA for gene-level and competitive gene-set analysis (Aim 1a).

Two steps: (1) SNP -> gene aggregation using the reference LD panel; (2) competitive gene-set
test for the pre-registered circadian set plus the negative/positive controls, so every reported
p-value sits next to its empirical null.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeneSetResult:
    set_name: str
    n_genes: int
    beta: float
    beta_se: float
    p: float


def gene_analysis(sumstats_pval_file: str | Path, n: int, ref_bfile: str | Path,
                  gene_loc: str | Path, out: str | Path, magma: str = "magma") -> Path:
    """Run MAGMA SNP-wise gene analysis; returns the .genes.raw path for the set step."""
    out = Path(out)
    subprocess.run(
        [magma, "--bfile", str(ref_bfile), "--gene-annot", f"{gene_loc}",
         "--pval", str(sumstats_pval_file), f"N={n}", "--out", str(out)],
        check=True,
    )
    return out.with_suffix(".genes.raw")


def geneset_analysis(genes_raw: str | Path, set_file: str | Path, out: str | Path,
                     magma: str = "magma") -> Path:
    """Competitive gene-set test. ``set_file`` lists each set as 'SETNAME gene1 gene2 ...'."""
    out = Path(out)
    subprocess.run(
        [magma, "--gene-results", str(genes_raw), "--set-annot", str(set_file),
         "--out", str(out)],
        check=True,
    )
    return out.with_suffix(".gsa.out")


def parse_gsa_out(path: str | Path) -> list[GeneSetResult]:
    """Parse MAGMA .gsa.out (whitespace table with columns VARIABLE NGENES BETA BETA_STD SE P)."""
    results: list[GeneSetResult] = []
    for ln in Path(path).read_text().splitlines():
        if not ln.strip() or ln.startswith("#") or ln.startswith("VARIABLE"):
            continue
        parts = ln.split()
        # VARIABLE TYPE NGENES BETA BETA_STD SE P  (column layout can vary by MAGMA version)
        try:
            name = parts[0]
            ngenes = int(parts[2])
            beta = float(parts[3])
            se = float(parts[-2])
            p = float(parts[-1])
        except (IndexError, ValueError):
            continue
        results.append(GeneSetResult(name, ngenes, beta, se, p))
    if not results:
        raise ValueError(f"no gene-set rows parsed from {path}")
    return results
