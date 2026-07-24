"""Compose the cortical surface and the cerebellar flatmap into one LNM-style panel.

``combine_surface_flatmap`` is the deliverable the downstream figure sections (S2–S6)
call: it resolves the colormap, symmetric limits, and threshold **once** and hands the
identical scale to both the cortex (:func:`plot_map_surface`) and the SUIT flatmap
(:func:`plot_map_flatmap`), then draws a single shared colorbar. One volumetric map in,
one 300-dpi PDF/PNG out, matched exactly as the tremor (Goede) and gait (Luo) figures do.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

from ._base import RenderResult, apply_style, resolve_scale
from .cerebellum import plot_map_flatmap
from .surface import plot_map_surface

__all__ = ["combine_surface_flatmap"]


def combine_surface_flatmap(
    nifti: str | Path | nib.Nifti1Image,
    *,
    cmap="lnm_cold_hot",
    vmin: float | None = None,
    vmax: float | None = None,
    threshold: float | None = None,
    percentile: float = 98.0,
    round_to: float | None = None,
    hemis: tuple[str, ...] = ("left", "right"),
    views: tuple[str, ...] = ("lateral", "medial"),
    mesh: str = "fsaverage",
    title: str | None = None,
    cbar_label: str = "connectivity (t)",
    flatmap_label: str = "cerebellum (SUIT flatmap)",
    output: str | Path | None = None,
    dpi: int = 300,
    figsize: tuple[float, float] = (7.4, 4.2),
    cortex_frac: float = 0.62,
) -> RenderResult:
    """Render ``nifti`` as cortex (2×2 inflated views) + cerebellar flatmap inset.

    All display parameters are resolved once and shared, so the two anatomies are strictly
    comparable. When ``output`` is given, writes ``output.png`` and ``output.pdf`` at
    ``dpi``. Returns a :class:`RenderResult` for the whole panel.
    """
    apply_style()
    # Resolve the shared scale ONCE (symmetric-or-asymmetric limits, recentered cmap for
    # the colorbar). Children receive the resolved limits + the ORIGINAL cmap and recenter
    # it themselves against the same limits — so all three stay identical without double-
    # warping an already-recentered colormap.
    sc = resolve_scale(nifti, cmap=cmap, vmin=vmin, vmax=vmax, threshold=threshold,
                       percentile=percentile, round_to=round_to)
    vmin, vmax = sc.vmin, sc.vmax

    fig = plt.figure(figsize=figsize)
    # Reserve vertical room for multi-line titles: the "Left"/"Right" column labels sit just
    # above the cortex grid (top + 0.012), so a 2+ line suptitle (anchored at y=0.98, growing
    # downward) would land on them. Push the grid — and thus those labels — down one title-line
    # worth per extra line. One-line titles are unchanged.
    n_title_lines = (title.count("\n") + 1) if title else 1
    grid_top = 0.89 - 0.05 * (n_title_lines - 1)
    outer = fig.add_gridspec(
        1, 2, width_ratios=[cortex_frac, 1 - cortex_frac],
        left=0.01, right=0.99, top=grid_top, bottom=0.14, wspace=0.02,
    )

    # --- cortex: 2×2 grid of 3-D axes in the left cell -----------------------
    nrow, ncol = len(views), len(hemis)
    cgs = outer[0, 0].subgridspec(nrow, ncol, wspace=-0.02, hspace=-0.08)
    cortex_axes = np.empty((nrow, ncol), dtype=object)
    for i in range(nrow):
        for j in range(ncol):
            cortex_axes[i, j] = fig.add_subplot(cgs[i, j], projection="3d")
    plot_map_surface(
        nifti, cmap=cmap, vmin=vmin, vmax=vmax, threshold=threshold, mesh=mesh,
        hemis=hemis, views=views, colorbar=False, figure=fig, axes=cortex_axes,
    )

    # --- cerebellum: flatmap centered in the right cell ----------------------
    fgs = outer[0, 1].subgridspec(3, 1, height_ratios=[0.12, 1.0, 0.12])
    flat_ax = fig.add_subplot(fgs[1, 0])
    plot_map_flatmap(
        nifti, cmap=cmap, vmin=vmin, vmax=vmax, threshold=threshold,
        colorbar=False, title=flatmap_label, figure=fig, axes=flat_ax,
    )

    result = RenderResult(fig, list(cortex_axes.ravel()) + [flat_ax], sc.cmap, sc.norm,
                          vmin, vmax, threshold)
    _shared_colorbar(fig, result, cbar_label)
    if title:
        fig.suptitle(title, fontsize=11, fontweight="bold", y=0.98)

    if output is not None:
        result.save(output, dpi=dpi)
    return result


def _shared_colorbar(fig, result: RenderResult, label: str) -> None:
    cax = fig.add_axes([0.30, 0.075, 0.40, 0.024])
    cb = fig.colorbar(result.mappable, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_ticks([result.vmin, 0, result.vmax])
    cb.set_label(label, fontsize=8)
