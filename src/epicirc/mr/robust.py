"""Robust two-sample MR estimators: weighted median, MR-PRESSO-style outlier removal, Steiger.

These stress-test an IVW estimate against invalid instruments. Pure Python; the weighted median
and outlier logic follow Bowden 2016 / Verbanck 2018, with documented simplifications (the
PRESSO global test uses the analytic Cochran-Q outlier contribution rather than the full
simulation; Steiger uses a continuous variance-explained approximation for the binary outcome).
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass


@dataclass(frozen=True)
class RobustEstimate:
    method: str
    beta: float
    se: float
    p: float
    n_snps: int
    extra: dict


def _p_from_z(z: float) -> float:
    return 2.0 * (0.5 * math.erfc(abs(z) / math.sqrt(2.0)))


def _wm_point(betas: list[float], weights: list[float]) -> float:
    """Weighted-median point estimate of a list of ratio estimates (Bowden 2016)."""
    order = sorted(range(len(betas)), key=lambda i: betas[i])
    b = [betas[i] for i in order]
    w = [weights[i] for i in order]
    sw = sum(w)
    cum = 0.0
    p = []
    for wi in w:
        p.append((cum + wi / 2.0) / sw)
        cum += wi
    # find crossing of 0.5
    for k in range(len(p)):
        if p[k] >= 0.5:
            if k == 0:
                return b[0]
            frac = (0.5 - p[k - 1]) / (p[k] - p[k - 1])
            return b[k - 1] + (b[k] - b[k - 1]) * frac
    return b[-1]


def weighted_median(bx, by, byse, n_boot: int = 1000, seed: int = 20260702) -> RobustEstimate:
    ratios = [y / x for x, y in zip(bx, by)]      # by/bx (Wald ratio per instrument)
    weights = [(x ** 2) / (s ** 2) for x, s in zip(bx, byse)]
    est = _wm_point(ratios, weights)
    rng = random.Random(seed)
    boots = []
    se_ratio = [s / abs(x) for x, s in zip(bx, byse)]
    for _ in range(n_boot):
        rb = [rng.gauss(r, se) for r, se in zip(ratios, se_ratio)]
        boots.append(_wm_point(rb, weights))
    mean = sum(boots) / len(boots)
    se = math.sqrt(sum((x - mean) ** 2 for x in boots) / (len(boots) - 1))
    z = est / se if se > 0 else 0.0
    return RobustEstimate("weighted_median", est, se, _p_from_z(z), len(bx), {})


def _ivw_point(bx, by, byse):
    den = sum(x * x / s ** 2 for x, s in zip(bx, byse))
    num = sum(x * y / s ** 2 for x, y, s in zip(bx, by, byse))
    beta = num / den
    se = math.sqrt(1.0 / den)
    return beta, se


def mr_presso_lite(bx, by, byse, alpha: float = 0.05) -> RobustEstimate:
    """Flag instruments whose Q-contribution is a significant outlier; return outlier-corrected IVW."""
    beta, _ = _ivw_point(bx, by, byse)
    n = len(bx)
    q_contrib = [(y - beta * x) ** 2 / s ** 2 for x, y, s in zip(bx, by, byse)]
    thresh = _chi2_1df_quantile(1 - alpha / n)  # Bonferroni across instruments
    outliers = [i for i, q in enumerate(q_contrib) if q > thresh]
    keep = [i for i in range(n) if i not in set(outliers)]
    if len(keep) < 3:
        keep = list(range(n))
        outliers = []
    bxk = [bx[i] for i in keep]; byk = [by[i] for i in keep]; bsk = [byse[i] for i in keep]
    b2, se2 = _ivw_point(bxk, byk, bsk)
    z = b2 / se2
    global_q = sum(q_contrib)
    return RobustEstimate("mr_presso_lite", b2, se2, _p_from_z(z), len(keep),
                          {"n_outliers": len(outliers), "global_Q": global_q, "Q_df": n - 1})


def steiger_filter(bx, bxse, eaf, n_exp, by, byse, n_out) -> RobustEstimate:
    """Keep instruments where exposure variance-explained > outcome variance-explained, then IVW.

    r2 ~ 2*eaf*(1-eaf)*beta^2 (continuous approx; applied to both traits). Guards against reverse
    causation / mismeasured direction.
    """
    keep = []
    for i in range(len(bx)):
        f = eaf[i]
        r2_exp = 2 * f * (1 - f) * bx[i] ** 2
        r2_out = 2 * f * (1 - f) * by[i] ** 2
        if r2_exp >= r2_out:
            keep.append(i)
    if len(keep) < 3:
        keep = list(range(len(bx)))
    bxk = [bx[i] for i in keep]; byk = [by[i] for i in keep]; bsk = [byse[i] for i in keep]
    b, se = _ivw_point(bxk, byk, bsk)
    z = b / se
    return RobustEstimate("steiger_ivw", b, se, _p_from_z(z), len(keep),
                          {"n_removed": len(bx) - len(keep)})


def _chi2_1df_quantile(p: float) -> float:
    """Inverse CDF of chi-square(1df) = (probit((1+p)/2))^2, via a rational approximation."""
    z = _probit((1 + p) / 2)
    return z * z


def _probit(p: float) -> float:
    """Acklam's inverse-normal approximation."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
