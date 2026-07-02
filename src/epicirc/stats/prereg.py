"""Generate pre-registration planning tables BEFORE any outcome is analysed (docs/adr/0003).

The key artifact is ``estimable_cells.tsv``: for each epilepsy phenotype it reports the effective
sample size and an *expected* LDSC heritability z-score, then flags whether the cell is estimable
(z > gate). This fixes the multiple-testing family to estimable cells only, so underpowered
subtypes are declared "not estimable — insufficient power" up front rather than mined as nulls.

**Approximation.** LDSC's h2 standard error scales, to first order, as SE(h2_obs) ~ M / N_eff,
so the expected z-score is

    z_h2  ~  N_eff * h2_liability / sqrt(2 * M)

with M the number of LD-independent SNPs (~1.1e6 HapMap3). This is a *planning* estimate only;
the actual z is computed by LDSC at runtime. Calibration check (h2_liab as configured below):
GGE (N_eff~26k) -> z~5, focal (~50k) -> z~7 (estimable); JME borderline; CAE/JAE -> z<2
(not estimable) — matching the adversarial power analysis.
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import yaml

M_LD_INDEP = 1.1e6            # ~HapMap3 LD-independent SNPs used by LDSC
DEFAULT_H2_GATE = 4.0

# Assumed liability-scale SNP-h2 by phenotype class (planning priors; cited in the manuscript).
# Focal/all-epilepsy lower; GGE and generalized syndromes higher (well-documented).
H2_LIABILITY_PRIOR = {
    "all_epilepsy": 0.18,
    "focal_epilepsy": 0.18,
    "focal_hs": 0.20,
    "focal_lesion_neg": 0.20,
    "focal_other_lesion": 0.18,
    "generalized_gge": 0.32,
    "jme": 0.40,
    "cae": 0.40,
    "jae": 0.40,
    "gtcs_alone": 0.35,
}
DEFAULT_H2_PRIOR = 0.20


@dataclass(frozen=True)
class CellPlan:
    phenotype: str
    n_cases: int
    n_controls: int
    n_eff: float
    h2_liab_prior: float
    expected_h2_z: float
    estimable: bool


def effective_n(n_cases: int, n_controls: int) -> float:
    """Effective sample size for a case/control study: 4 / (1/Ncase + 1/Ncontrol)."""
    if n_cases <= 0 or n_controls <= 0:
        raise ValueError("case/control counts must be positive")
    return 4.0 / (1.0 / n_cases + 1.0 / n_controls)


def expected_h2_zscore(n_eff: float, h2_liab: float, m: float = M_LD_INDEP) -> float:
    """Planning approximation for the LDSC heritability z-score (see module docstring)."""
    if n_eff <= 0 or h2_liab <= 0:
        raise ValueError("n_eff and h2_liab must be positive")
    return n_eff * h2_liab / math.sqrt(2.0 * m)


def plan_cells(manifest: dict, gate: float = DEFAULT_H2_GATE) -> list[CellPlan]:
    cells: list[CellPlan] = []
    for name, meta in manifest.get("epilepsy", {}).items():
        nc, nk = int(meta["n_cases"]), int(meta["n_controls"])
        neff = effective_n(nc, nk)
        h2 = H2_LIABILITY_PRIOR.get(name, DEFAULT_H2_PRIOR)
        z = expected_h2_zscore(neff, h2)
        cells.append(CellPlan(name, nc, nk, neff, h2, z, z > gate))
    return sorted(cells, key=lambda c: c.expected_h2_z, reverse=True)


def write_estimable_cells(cells: list[CellPlan], out: str | Path) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    header = "phenotype\tn_cases\tn_controls\tn_eff\th2_liab_prior\texpected_h2_z\testimable\n"
    lines = [
        f"{c.phenotype}\t{c.n_cases}\t{c.n_controls}\t{c.n_eff:.0f}\t"
        f"{c.h2_liab_prior:.2f}\t{c.expected_h2_z:.2f}\t{'yes' if c.estimable else 'no'}\n"
        for c in cells
    ]
    out.write_text(header + "".join(lines))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default="config/traits_manifest.yaml")
    ap.add_argument("--out", default="results/estimable_cells.tsv")
    ap.add_argument("--gate", type=float, default=DEFAULT_H2_GATE)
    args = ap.parse_args(argv)

    manifest = yaml.safe_load(Path(args.manifest).read_text())
    cells = plan_cells(manifest, gate=args.gate)
    write_estimable_cells(cells, args.out)

    est = [c.phenotype for c in cells if c.estimable]
    not_est = [c.phenotype for c in cells if not c.estimable]
    print(f"Pre-registration estimable-cell plan -> {args.out} (gate z>{args.gate})")
    print(f"  ESTIMABLE ({len(est)}): {', '.join(est) or 'none'}")
    print(f"  NOT estimable ({len(not_est)}): {', '.join(not_est) or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
