"""Cold–hot diverging colormap matching the Fox/Horn LNM figures, plus a resolver.

The tremor (Goede) and gait (Luo) network figures render connectivity with the
neuroimaging-standard *cold–hot* scheme: negative connectivity in blue→cyan, positive
in red→yellow, with the sub-threshold band left transparent so the anatomical underlay
shows through. ``lnm_cold_hot`` reproduces that look self-containedly (matplotlib only);
:func:`resolve_cmap` additionally accepts nilearn's ``cold_hot`` / ``roy_big_bl`` (the
Connectome-Workbench colormap) and any matplotlib name.

Accessibility note: cold–hot avoids the red–green axis, so it is robust to the common
(protan/deutan) color-vision deficiencies, and each half-ramp is luminance-monotonic so
magnitude survives in grayscale (a mid-blue and a mid-red differ, but the **poles*
blue↔yellow are the tritan-weak pair — for a strictly colorblind-first figure pass
``cmap="RdBu_r"``). See docs/figure_style_notes.md.
"""
from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.colors import Colormap, LinearSegmentedColormap, ListedColormap

# Anchor colors on a 0→1 axis where 0 = most-negative, 0.5 = zero, 1 = most-positive.
# Cold half: deep blue → cyan.  Hot half: red → orange → yellow.  Near-white at center
# (hidden by thresholding in practice, but keeps an unthresholded map readable).
_COLD_HOT_ANCHORS = [
    (0.00, (0.13, 0.20, 0.51)),   # deep blue  (strong negative)
    (0.15, (0.16, 0.36, 0.78)),   # blue
    (0.32, (0.27, 0.60, 0.94)),   # medium blue
    (0.46, (0.62, 0.85, 0.99)),   # light cyan (just below zero)
    (0.50, (0.96, 0.96, 0.96)),   # near-white (zero)
    (0.54, (0.90, 0.35, 0.20)),   # red        (just above zero)
    (0.68, (0.96, 0.52, 0.09)),   # orange
    (0.85, (0.99, 0.74, 0.09)),   # amber
    (1.00, (1.00, 0.95, 0.28)),   # yellow     (strong positive)
]

lnm_cold_hot: LinearSegmentedColormap = LinearSegmentedColormap.from_list(
    "lnm_cold_hot",
    [(pos, rgb) for pos, rgb in _COLD_HOT_ANCHORS],
    N=256,
)
lnm_cold_hot.set_bad(alpha=0.0)  # masked (sub-threshold) voxels are transparent

# Register once so "lnm_cold_hot" works everywhere a name is accepted.
try:
    mpl.colormaps.register(lnm_cold_hot, force=True)
    mpl.colormaps.register(lnm_cold_hot.reversed(), force=True)
except (ValueError, AttributeError):  # older mpl or already registered
    pass


def resolve_cmap(cmap: str | Colormap | None) -> Colormap:
    """Return a matplotlib ``Colormap`` from a name, a ``Colormap``, or ``None``.

    Names are resolved in order: our custom map, matplotlib's registry, then nilearn's
    neuroimaging colormaps (``cold_hot``, ``roy_big_bl``, ``blue_red``, …). ``None`` and
    the sentinel ``"lnm_cold_hot"`` both give the paper-matched default.
    """
    if cmap is None:
        return lnm_cold_hot
    if isinstance(cmap, Colormap):
        return cmap
    if cmap == "lnm_cold_hot":
        return lnm_cold_hot
    try:
        return mpl.colormaps[cmap]
    except KeyError:
        pass
    try:  # nilearn ships the Fox-lab cold_hot and Workbench roy_big_bl
        from nilearn.plotting import cm as ncm

        cand = getattr(ncm, cmap, None)
        if isinstance(cand, Colormap):
            return cand
    except Exception:  # nilearn optional at colormap-resolution time
        pass
    raise ValueError(f"Unknown colormap: {cmap!r}")


def recenter_diverging(
    cmap: Colormap, vmin: float, vmax: float, vcenter: float = 0.0, n: int = 256
) -> Colormap:
    """Warp a diverging ``cmap`` so ``vcenter`` lands on its midpoint under linear scaling.

    A diverging colormap fixes its neutral color at position 0.5. Under a plain linear
    ``Normalize(vmin, vmax)`` with **asymmetric** limits, data value ``vcenter`` (0 for a
    signed map) maps to ``(vcenter-vmin)/(vmax-vmin) ≠ 0.5``, so zero no longer renders as
    the neutral color and the sign→hue contract breaks. Resampling the colormap so its two
    halves occupy ``[vmin, vcenter]`` and ``[vcenter, vmax]`` restores it — and, unlike
    ``TwoSlopeNorm``, works with backends (nilearn ``plot_surf_stat_map``, SUITPy) that only
    accept a colormap plus a linear ``[vmin, vmax]`` scale. A symmetric (or non-straddling)
    range is returned unchanged.
    """
    if not (vmin < vcenter < vmax):
        return cmap
    p = (vcenter - vmin) / (vmax - vmin)
    if abs(p - 0.5) < 1e-6:
        return cmap
    xs = np.linspace(0.0, 1.0, n)
    src = np.where(xs <= p, 0.5 * xs / p, 0.5 + 0.5 * (xs - p) / (1.0 - p))
    out = ListedColormap(cmap(src), name=f"{cmap.name}_recentered")
    out.set_bad(cmap.get_bad())
    return out


def auto_threshold(
    nifti_or_values,
    *,
    percentile: float = 75.0,
) -> float:
    """A magnitude threshold at the given percentile of ``|non-zero|`` voxel values.

    Convenience for suppressing the near-zero haze on unthresholded stat maps. Callers
    that already have a principled threshold (an FWE t-cut, |R| gate) should pass it
    explicitly; this is only a sensible default for exploratory rendering.
    """
    from ._base import map_values

    if isinstance(nifti_or_values, np.ndarray):
        d = nifti_or_values
    else:
        d = map_values(nifti_or_values)
    nz = np.abs(d[np.isfinite(d) & (d != 0)])
    return float(np.percentile(nz, percentile)) if nz.size else 0.0
