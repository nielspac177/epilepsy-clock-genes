"""Focal-vs-GGE heterogeneity test for a gene-set / genetic-correlation effect.

The core methodological subtlety (confirmed against the adversarial design review): the ILAE
focal and GGE GWAS are two case/control studies that **share the same controls**. Shared controls
induce a *positive* correlation ``rho`` between the two effect estimates. For the **difference**
b_focal - b_gge the shared control frequency largely cancels, so the true variance

    Var(b_focal - b_gge) = se_focal^2 + se_gge^2 - 2 * rho * se_focal * se_gge

is *smaller* than the independence sum se_focal^2 + se_gge^2. Ignoring rho (rho=0) therefore
makes the type-difference test **conservative** — it double-counts the control variance and
loses power to detect genuine focal-vs-GGE differences. It does not create false positives, so
rho=0 is a *safe* default; supplying the true rho recovers power.

``rho`` is not guessed. It is estimated empirically from the **cross-trait LDSC intercept**
between the focal and GGE GWAS (which captures the shared-sample correlation), passed in by the
pipeline. With rho unknown we fall back to rho=0 with a loud warning (conservative, power-losing).
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass


@dataclass(frozen=True)
class HeterogeneityResult:
    delta: float          # b_a - b_b
    se_delta: float
    z: float
    p_two_sided: float
    rho_used: float
    rho_assumed: bool     # True if rho defaulted to 0 (optimistic)


def _norm_sf(x: float) -> float:
    """Upper-tail standard-normal survival function (no scipy dependency)."""
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def heterogeneity_test(
    b_a: float,
    se_a: float,
    b_b: float,
    se_b: float,
    rho: float | None = None,
) -> HeterogeneityResult:
    """Test whether effect ``a`` (e.g. focal) differs from effect ``b`` (e.g. GGE).

    Parameters
    ----------
    b_a, se_a, b_b, se_b:
        Effect estimates and standard errors for the two epilepsy types.
    rho:
        Correlation between the two estimates induced by shared controls, from the cross-trait
        LDSC intercept. If ``None``, defaults to 0 with a warning (anti-conservative).
    """
    if se_a <= 0 or se_b <= 0:
        raise ValueError("standard errors must be positive")

    rho_assumed = rho is None
    if rho_assumed:
        warnings.warn(
            "rho not supplied; assuming independence (rho=0). For a shared-control case-case "
            "difference this is conservative (power-losing), not anti-conservative. Pass the "
            "cross-trait LDSC intercept to recover power.",
            stacklevel=2,
        )
        rho = 0.0
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")

    var_delta = se_a**2 + se_b**2 - 2.0 * rho * se_a * se_b
    if var_delta <= 0:
        raise ValueError("non-positive variance for the difference; check rho/SEs")

    delta = b_a - b_b
    se_delta = math.sqrt(var_delta)
    z = delta / se_delta
    p = 2.0 * _norm_sf(abs(z))
    return HeterogeneityResult(
        delta=delta, se_delta=se_delta, z=z, p_two_sided=p,
        rho_used=rho, rho_assumed=rho_assumed,
    )
