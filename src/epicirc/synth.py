"""Deterministic synthetic GWAS generator for end-to-end pipeline validation.

Produces SNP-level summary statistics with a *shared latent structure* across datasets so that
(a) genetic correlation and (b) a focal-vs-GGE circadian difference are demonstrable without any
real data or external tools. Everything is stdlib + seeded, so the whole DAG is reproducible.

Design: each SNP carries two latent factors shared across ALL datasets — a generic effect ``u``
and a circadian factor ``c`` (non-zero only for SNPs in circadian genes). A dataset's per-SNP
effect is ``u + delta*c (if circadian) + dataset noise``. Datasets with a larger circadian
``delta`` therefore correlate more strongly with circadian-exposure datasets and show larger
circadian gene-set enrichment. This is a *demonstration fixture*, not a simulation study.
"""
from __future__ import annotations

import math
import random
import zlib
from dataclasses import dataclass

N_GENES_BACKGROUND = 200
SNPS_PER_GENE = 10

# circadian genes injected as signal-carrying (matches config/gene_sets/circadian_core.tsv)
CIRCADIAN_GENES = (
    "ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
    "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP",
    "NFIL3", "TIMELESS", "BHLHE40", "BHLHE41", "CIART",
)

# circadian effect size (delta) by dataset — encodes the ground-truth hypothesis GGE > focal
CIRCADIAN_DELTA = {
    "all_epilepsy": 0.22,
    "focal_epilepsy": 0.10,
    "generalized_gge": 0.38,
    "jme": 0.40, "cae": 0.35, "jae": 0.33, "gtcs_alone": 0.34,
    "focal_hs": 0.11, "focal_lesion_neg": 0.10, "focal_other_lesion": 0.09,
    "drug_resistant_epilepsy": 0.16,
    # circadian exposures carry strong circadian signal (they *are* the clock phenotype)
    "chronotype_morning": 0.60, "insomnia": 0.45, "sleep_duration": 0.40, "daytime_napping": 0.42,
    # comorbid/confounder traits: modest circadian loading
    "mdd": 0.20, "bip": 0.18, "adhd": 0.17, "bmi": 0.12, "t2d": 0.10,
}
DEFAULT_DELTA = 0.05

_U_SEED = 10_000
_C_SEED = 5_000_000


@dataclass(frozen=True)
class SnpRecord:
    snp: str
    chrom: int
    bp: int
    a1: str
    a2: str
    frq: float
    beta: float
    se: float
    p: float
    n: int
    gene: str
    circadian: bool


def _latent_u(idx: int) -> float:
    return random.Random(_U_SEED + idx).gauss(0.0, 0.03)


def _latent_c(idx: int) -> float:
    return random.Random(_C_SEED + idx).gauss(0.0, 0.03)


def _norm_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def _gene_layout() -> list[tuple[str, bool]]:
    genes: list[tuple[str, bool]] = [(g, True) for g in CIRCADIAN_GENES]
    genes += [(f"BKGD{i:04d}", False) for i in range(N_GENES_BACKGROUND)]
    return genes


def generate(name: str, n_eff: int = 40000, delta: float | None = None) -> list[SnpRecord]:
    """Generate deterministic sumstats for dataset ``name``."""
    d = CIRCADIAN_DELTA.get(name, DEFAULT_DELTA) if delta is None else delta
    noise_rng = random.Random(zlib.crc32(name.encode()))   # stable across processes
    alleles = [("A", "G"), ("C", "T"), ("A", "C"), ("T", "G")]
    recs: list[SnpRecord] = []
    idx = 0
    for chrom, (gene, is_circ) in enumerate(_gene_layout(), start=1):
        for k in range(SNPS_PER_GENE):
            u, c = _latent_u(idx), _latent_c(idx)
            frq = 0.05 + noise_rng.random() * 0.45
            se = 1.0 / math.sqrt(2.0 * frq * (1 - frq) * max(n_eff, 1))
            eff = u + (d * c if is_circ else 0.0) + noise_rng.gauss(0.0, se)
            z = eff / se
            a1, a2 = alleles[idx % len(alleles)]
            recs.append(SnpRecord(
                snp=f"rs{idx:06d}", chrom=(chrom % 22) + 1, bp=(k + 1) * 1000,
                a1=a1, a2=a2, frq=round(frq, 4), beta=round(eff, 6), se=round(se, 6),
                p=max(2.0 * _norm_sf(abs(z)), 1e-300), n=n_eff, gene=gene, circadian=is_circ,
            ))
            idx += 1
    return recs
