"""Two-sample Mendelian randomization estimators.

Closed-form IVW and MR-Egger are implemented in pure Python (testable, no R needed). Heavier
sensitivity analyses (weighted median, MR-PRESSO global test, Steiger filtering with proper SEs)
are delegated to the TwoSampleMR R package via ``run_twosamplemr_r`` so we use the validated
reference implementation rather than re-deriving it.

Conventions: ``bx``/``sex`` are instrument-exposure effects/SEs, ``by``/``sey`` outcome effects/SEs.
Instruments must already be harmonized (see epicirc.data.harmonize) and LD-clumped.
"""
from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class MREstimate:
    method: str
    beta: float
    se: float
    z: float
    p: float
    n_snps: int
    extra: dict


def _norm_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def _p_from_z(z: float) -> float:
    return 2.0 * _norm_sf(abs(z))


def ivw(bx: Sequence[float], by: Sequence[float], sey: Sequence[float]) -> MREstimate:
    """Fixed-effect inverse-variance-weighted estimator.

    Assumes no measurement error in bx (NOME); the standard first-line MR estimator.
    """
    _validate(bx, by, sey)
    num = sum(x * y / s**2 for x, y, s in zip(bx, by, sey))
    den = sum(x**2 / s**2 for x, s in zip(bx, sey))
    if den <= 0:
        raise ValueError("degenerate weights in IVW")
    beta = num / den
    se = math.sqrt(1.0 / den)
    z = beta / se
    # Cochran's Q for instrument heterogeneity (pleiotropy indicator)
    q = sum((y - beta * x) ** 2 / s**2 for x, y, s in zip(bx, by, sey))
    return MREstimate("ivw", beta, se, z, _p_from_z(z), len(bx), {"cochran_q": q, "q_df": len(bx) - 1})


def mr_egger(bx: Sequence[float], by: Sequence[float], sey: Sequence[float]) -> MREstimate:
    """MR-Egger weighted regression; intercept estimates directional pleiotropy.

    Returns the slope as the causal estimate and reports the intercept (and its p-value) in
    ``extra`` — a non-zero intercept flags horizontal pleiotropy, the key MR threat the
    adversarial review emphasized.
    """
    _validate(bx, by, sey)
    n = len(bx)
    if n < 3:
        raise ValueError("MR-Egger needs >=3 instruments")
    w = [1.0 / s**2 for s in sey]
    sw = sum(w)
    mx = sum(wi * x for wi, x in zip(w, bx)) / sw
    my = sum(wi * y for wi, y in zip(w, by)) / sw
    sxx = sum(wi * (x - mx) ** 2 for wi, x in zip(w, bx))
    sxy = sum(wi * (x - mx) * (y - my) for wi, x, y in zip(w, bx, by))
    slope = sxy / sxx
    intercept = my - slope * mx
    # residual variance for SEs
    resid = [y - intercept - slope * x for x, y in zip(bx, by)]
    dof = n - 2
    sigma2 = sum(wi * r**2 for wi, r in zip(w, resid)) / dof
    se_slope = math.sqrt(sigma2 / sxx)
    se_int = math.sqrt(sigma2 * (1.0 / sw + mx**2 / sxx))
    z = slope / se_slope
    z_int = intercept / se_int
    return MREstimate(
        "mr_egger", slope, se_slope, z, _p_from_z(z), n,
        {"intercept": intercept, "intercept_se": se_int,
         "intercept_p": _p_from_z(z_int), "pleiotropy_flag": _p_from_z(z_int) < 0.05},
    )


def run_twosamplemr_r(harmonized_tsv: str | Path, out_dir: str | Path, rscript: str = "Rscript") -> Path:
    """Delegate full MR (weighted median, MR-PRESSO, Steiger) to the R script.

    Requires R with the TwoSampleMR + MRPRESSO packages installed (see environment.yml). Kept as
    a subprocess boundary so the reference implementation is authoritative.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("twosamplemr_runner.R")
    subprocess.run([rscript, str(script), str(harmonized_tsv), str(out_dir)], check=True)
    return out_dir / "mr_results.tsv"


def _validate(bx, by, sey) -> None:
    if not (len(bx) == len(by) == len(sey)):
        raise ValueError("bx, by, sey must be equal length")
    if len(bx) < 2:
        raise ValueError("need >=2 instruments")
    if any(s <= 0 for s in sey):
        raise ValueError("outcome SEs must be positive")
