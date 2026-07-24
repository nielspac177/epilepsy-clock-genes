"""Shared plumbing for the LNM renderers: I/O, symmetric limits, and a result type.

The three public renderers (``surface``, ``cerebellum``, ``combine``) all speak the
same small vocabulary defined here so a map drawn on the cortex and the same map drawn
on the cerebellum share one colormap, one symmetric scale, and one threshold — which is
the whole point of rendering a volumetric map "in the Fox/Horn LNM style".
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.colors import Normalize

# ---------------------------------------------------------------------------
# Publication style — sans-serif, thin spines, matches pipeline/06 conventions.
# Applied lazily (not on import) so importing the package never mutates a caller's
# global rcParams unexpectedly; call apply_style() from the entry points.
# ---------------------------------------------------------------------------
_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # editable text in Illustrator (TrueType, not Type-3)
    "ps.fonttype": 42,
}


def apply_style() -> None:
    """Set the shared publication rcParams (idempotent)."""
    mpl.rcParams.update(_STYLE)


def load_map(nifti: str | Path | nib.Nifti1Image) -> nib.Nifti1Image:
    """Accept a path or an already-loaded image; return a NIfTI image."""
    if isinstance(nifti, (str, Path)):
        return nib.load(str(nifti))
    return nifti


def map_values(nifti: str | Path | nib.Nifti1Image) -> np.ndarray:
    """Finite voxel values of a map with NaN/inf coerced to 0 (as displayed)."""
    d = load_map(nifti).get_fdata()
    return np.nan_to_num(d, nan=0.0, posinf=0.0, neginf=0.0)


def symmetric_limits(
    nifti: str | Path | nib.Nifti1Image | np.ndarray,
    *,
    percentile: float = 98.0,
    round_to: float | None = None,
    minimum: float = 1e-6,
) -> tuple[float, float]:
    """Robust symmetric display limits ``(-v, +v)`` for a signed statistical map.

    ``v`` is the ``percentile`` of the absolute non-zero voxel values, so a few extreme
    voxels do not blow out the scale. Symmetric limits keep the diverging colormap
    centered on zero — the honest choice for a signed map where +2 and −2 must read as
    equal-magnitude opposites.

    Rounding of ``v`` to a clean colorbar limit: ``round_to=None`` (default) rounds up to a
    step **scaled to the magnitude** of ``v`` (a t-map ~3.9 → 4.0; an r-map ~0.28 → 0.30),
    so it never flattens a small-valued map; a positive ``round_to`` forces that fixed step;
    ``0`` disables rounding.
    """
    d = nifti if isinstance(nifti, np.ndarray) else map_values(nifti)
    nz = np.abs(d[np.isfinite(d) & (d != 0)])
    v = float(np.percentile(nz, percentile)) if nz.size else minimum
    if v > 0:
        if round_to is None:                       # adaptive: ~1/2 order of magnitude
            step = 0.5 * 10.0 ** np.floor(np.log10(v))
            v = float(np.ceil(v / step) * step)
        elif round_to > 0:
            v = float(np.ceil(v / round_to) * round_to)
    return -max(v, minimum), max(v, minimum)


def make_norm(vmin: float, vmax: float) -> Normalize:
    return Normalize(vmin=vmin, vmax=vmax)


def normalize_threshold(threshold) -> tuple[float, float] | None:
    """Coerce a threshold into a ``(neg, pos)`` band, or ``None``.

    A scalar ``t`` becomes the symmetric band ``(-|t|, +|t|)``; a 2-tuple is honored as
    an asymmetric ``(neg, pos)`` band (the papers threshold each side independently, e.g.
    positive from ``+0.7`` but negative from ``-1.0``). Voxels inside the band are hidden.
    """
    if threshold is None:
        return None
    if np.isscalar(threshold):
        t = abs(float(threshold))
        return (-t, t)
    lo, hi = float(threshold[0]), float(threshold[1])
    return (min(lo, hi), max(lo, hi))


def _resolve_symmetric(
    nifti,
    vmin: float | None,
    vmax: float | None,
    *,
    percentile: float,
    round_to: float | None,
) -> tuple[float, float]:
    """Fill in missing limits from the data, forcing symmetry about zero."""
    if vmin is not None and vmax is not None:
        return float(vmin), float(vmax)
    lo, hi = symmetric_limits(nifti, percentile=percentile, round_to=round_to)
    if vmax is not None:
        return -abs(float(vmax)), abs(float(vmax))
    if vmin is not None:
        return -abs(float(vmin)), abs(float(vmin))
    return lo, hi


@dataclass
class Scale:
    """The resolved display scale shared by every renderer (one source of truth)."""

    cmap: mpl.colors.Colormap
    vmin: float
    vmax: float
    norm: Normalize
    threshold_band: tuple[float, float] | None


def resolve_scale(
    nifti,
    *,
    cmap,
    vmin: float | None,
    vmax: float | None,
    threshold,
    percentile: float,
    round_to: float | None,
) -> Scale:
    """Resolve colormap + symmetric-or-asymmetric limits + threshold band, once.

    Centralizes the choice so the cortex, the cerebellum, and the colorbar in a combined
    panel are guaranteed identical. When the resolved limits straddle zero **asymmetrically**
    the diverging colormap is recentered on zero (see
    :func:`~pipeline.render.colormaps.recenter_diverging`) so a signed value's sign always
    maps to the correct hue.
    """
    from .colormaps import recenter_diverging, resolve_cmap  # lazy: avoid import cycle

    cmap = resolve_cmap(cmap)
    vmin, vmax = _resolve_symmetric(nifti, vmin, vmax, percentile=percentile, round_to=round_to)
    if vmin < 0 < vmax and abs(vmin) != abs(vmax):
        cmap = recenter_diverging(cmap, vmin, vmax, 0.0)
    return Scale(cmap, vmin, vmax, make_norm(vmin, vmax), normalize_threshold(threshold))


@dataclass
class RenderResult:
    """A rendered figure plus the scale that produced it.

    Carrying ``cmap``/``norm``/``threshold`` alongside the figure is what lets
    :func:`combine_surface_flatmap` draw one shared colorbar for a cortex panel and a
    cerebellum panel that were rendered by two different toolkits (nilearn, SUITPy).
    """

    figure: plt.Figure
    axes: list = field(default_factory=list)
    cmap: mpl.colors.Colormap | None = None
    norm: Normalize | None = None
    vmin: float | None = None
    vmax: float | None = None
    threshold: float | None = None

    @property
    def mappable(self) -> mpl.cm.ScalarMappable:
        sm = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
        sm.set_array([])
        return sm

    def save(
        self,
        stem: str | Path,
        *,
        formats: tuple[str, ...] = ("png", "pdf"),
        dpi: int = 300,
        transparent: bool = False,
    ) -> list[Path]:
        """Write ``stem.png``/``stem.pdf`` (300 dpi by default). Returns the paths."""
        stem = Path(stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        out = []
        for ext in formats:
            p = stem.with_suffix(f".{ext}")
            self.figure.savefig(
                p, dpi=dpi, bbox_inches="tight", transparent=transparent,
                facecolor=self.figure.get_facecolor(),
            )
            out.append(p)
        return out

    def close(self) -> None:
        plt.close(self.figure)
