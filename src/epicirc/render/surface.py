"""Render an MNI volumetric map on the inflated fsaverage cortical surface.

``plot_map_surface`` reproduces the cortical half of the Fox/Horn LNM figure: the map is
projected from the volume onto the fsaverage surface (sampling between the white and pial
meshes) and shown on the inflated surface as a 2×2 grid of L/R × lateral/medial views,
shaded by sulcal depth, with a matched cold–hot diverging colormap and threshold.
"""
from __future__ import annotations

import inspect
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from nilearn import datasets, plotting, surface

from ._base import RenderResult, apply_style, load_map, resolve_scale

__all__ = ["plot_map_surface", "project_to_fsaverage"]

# nilearn kwargs vary by version (e.g. `darkness` removed in 0.13); introspect what's accepted.
try:
    _PLOT_SURF_PARAMS = set(inspect.signature(plotting.plot_surf_stat_map).parameters)
except (ValueError, TypeError):
    _PLOT_SURF_PARAMS = set()

_VIEW_ANGLES = {"lateral", "medial", "dorsal", "ventral", "anterior", "posterior"}
_HEMIS = {"left", "right"}


def _reject_unknown(tokens, allowed: set, kind: str) -> None:
    bad = [t for t in tokens if t not in allowed]
    if bad:
        raise ValueError(f"Unknown {kind}(s) {bad}; choose from {sorted(allowed)}")


@lru_cache(maxsize=4)
def _fsaverage(mesh: str):
    """Fetch (and cache) an fsaverage surface bundle. 'fsaverage' = high-res (164k)."""
    return datasets.fetch_surf_fsaverage(mesh)


def project_to_fsaverage(
    nifti: str | Path | nib.Nifti1Image,
    *,
    mesh: str = "fsaverage",
    hemis: tuple[str, ...] = ("left", "right"),
    interpolation: str = "linear",
    radius: float = 3.0,
) -> dict[str, np.ndarray]:
    """Project a volume to fsaverage vertices per hemisphere (white↔pial sampling).

    Returns ``{"left": tex, "right": tex}``. Sampling between the white and pial meshes
    (rather than at a single depth) is the nilearn-recommended way to avoid missing
    gyral-crown or sulcal-fundus signal.
    """
    _reject_unknown(hemis, _HEMIS, "hemi")
    img = load_map(nifti)
    fs = _fsaverage(mesh)
    out = {}
    for h in hemis:
        tex = surface.vol_to_surf(
            img, fs[f"pial_{h}"], inner_mesh=fs[f"white_{h}"],
            interpolation=interpolation, radius=radius,
        )
        out[h] = np.nan_to_num(tex)
    return out


def plot_map_surface(
    nifti: str | Path | nib.Nifti1Image,
    *,
    cmap="lnm_cold_hot",
    vmin: float | None = None,
    vmax: float | None = None,
    threshold: float | None = None,
    mesh: str = "fsaverage",
    surf_kind: str = "infl",
    hemis: tuple[str, ...] = ("left", "right"),
    views: tuple[str, ...] = ("lateral", "medial"),
    bg_on_data: bool = True,
    darkness: float = 0.5,
    alpha: float = 1.0,
    percentile: float = 98.0,
    round_to: float | None = None,
    grid_labels: bool = True,
    colorbar: bool = True,
    cbar_label: str | None = None,
    title: str | None = None,
    figure: plt.Figure | None = None,
    axes: np.ndarray | list | None = None,
    figsize: tuple[float, float] | None = None,
) -> RenderResult:
    """Render ``nifti`` on the inflated cortex as a ``len(views)`` × ``len(hemis)`` grid.

    Parameters mirror the paper style: ``cmap`` cold–hot diverging, symmetric ``vmin``/
    ``vmax`` (auto from a robust percentile when omitted), and a magnitude ``threshold``
    that leaves the sub-threshold band transparent. Returns a :class:`RenderResult`
    carrying the figure and the exact scale (so a shared colorbar can be reused).

    Pass ``figure``+``axes`` (a grid of 3-D axes shaped ``len(views)``×``len(hemis)``) to
    draw into an existing composite; otherwise a new figure is created. Set
    ``grid_labels=False`` to suppress the per-grid hemisphere/view labels when a composite
    labels the shared geometry once itself (e.g. a montage of several panels).
    """
    apply_style()
    _reject_unknown(views, _VIEW_ANGLES, "view")
    _reject_unknown(hemis, _HEMIS, "hemi")
    sc = resolve_scale(nifti, cmap=cmap, vmin=vmin, vmax=vmax, threshold=threshold,
                       percentile=percentile, round_to=round_to)
    cmap, vmin, vmax, norm, thr_band = sc.cmap, sc.vmin, sc.vmax, sc.norm, sc.threshold_band

    textures = project_to_fsaverage(nifti, mesh=mesh, hemis=hemis)
    fs = _fsaverage(mesh)
    nrow, ncol = len(views), len(hemis)

    owns_figure = axes is None
    if axes is None:
        figsize = figsize or (2.6 * ncol, 2.4 * nrow)
        figure = figure or plt.figure(figsize=figsize)
        axes = np.empty((nrow, ncol), dtype=object)
        for i in range(nrow):
            for j in range(ncol):
                axes[i, j] = figure.add_subplot(nrow, ncol, i * ncol + j + 1, projection="3d")
    else:
        axes = np.asarray(axes, dtype=object).reshape(nrow, ncol)
        figure = figure or axes[0, 0].get_figure()

    for i, view in enumerate(views):
        for j, hemi in enumerate(hemis):
            ax = axes[i, j]
            tex = textures[hemi]
            nthr = None
            if thr_band is not None:
                lo, hi = thr_band
                if abs(lo) == abs(hi):          # symmetric -> nilearn's native path
                    nthr = abs(hi)
                else:                            # asymmetric -> hide the band ourselves
                    tex = tex.copy()
                    tex[(tex > lo) & (tex < hi)] = np.nan
            surf_kwargs = dict(
                bg_map=fs[f"sulc_{hemi}"],
                hemi=hemi,
                view=view,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                threshold=nthr,
                symmetric_cbar=abs(vmin) == abs(vmax),
                bg_on_data=bg_on_data,
                alpha=alpha,
                colorbar=False,
                axes=ax,
                figure=figure,
            )
            # `darkness` was deprecated in nilearn 0.11 and removed in 0.13; pass it only if supported.
            if "darkness" in _PLOT_SURF_PARAMS:
                surf_kwargs["darkness"] = darkness
            plotting.plot_surf_stat_map(
                fs[f"{surf_kind}_{hemi}"],
                tex,
                **surf_kwargs,
            )
            ax.set_title("")
    _zoom_3d(axes)
    if owns_figure:
        figure.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.10,
                               wspace=-0.05, hspace=-0.05)
    if grid_labels:
        _grid_labels(figure, axes, hemis, views)

    result = RenderResult(figure, list(axes.ravel()), cmap, norm, vmin, vmax, threshold)
    if colorbar:
        _add_colorbar(figure, result, cbar_label)
    if title:
        figure.suptitle(title, fontsize=10, fontweight="bold", y=0.985)

    return result


def _zoom_3d(axes: np.ndarray) -> None:
    """Zoom each 3-D surface axis to remove nilearn's default whitespace margins."""
    for ax in axes.ravel():
        try:
            ax.set_box_aspect(None, zoom=1.35)
        except Exception:
            pass
        ax.margins(0)


def _grid_labels(figure, axes: np.ndarray, hemis, views) -> None:
    """Label columns with hemisphere and rows with view, from the axes' geometry.

    Uses the subplot cell rectangles (robust for 3-D axes and for a sub-grid embedded in
    a composite figure), so the same call labels a standalone cortex grid and the cortex
    block inside :func:`combine_surface_flatmap`.
    """
    nrow, ncol = axes.shape
    pos = np.vectorize(lambda a: a.get_position(), otypes=[object])(axes)
    top = max(pos[0, j].y1 for j in range(ncol))
    left = min(pos[i, 0].x0 for i in range(nrow))
    for j, h in enumerate(hemis):
        xc = float(np.mean([pos[i, j].x0 + pos[i, j].width / 2 for i in range(nrow)]))
        figure.text(xc, min(top + 0.012, 0.995), h.capitalize(),
                    ha="center", va="bottom", fontsize=9)
    if nrow > 1:
        for i, v in enumerate(views):
            yc = float(np.mean([pos[i, j].y0 + pos[i, j].height / 2 for j in range(ncol)]))
            figure.text(max(left - 0.006, 0.004), yc, v,
                        ha="left", va="center", rotation=90, fontsize=8)


def _add_colorbar(figure, result: RenderResult, label: str | None) -> None:
    cax = figure.add_axes([0.30, 0.045, 0.40, 0.022])
    cb = figure.colorbar(result.mappable, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_ticks([result.vmin, 0, result.vmax])
    if label:
        cb.set_label(label, fontsize=8)
