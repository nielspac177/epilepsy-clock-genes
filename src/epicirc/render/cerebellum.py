"""Render an MNI volumetric map on the Diedrichsen SUIT cerebellar flatmap (SUITPy).

The cerebellum is where MRgFUS-VIM networks live, and the LNM papers show it as a SUIT
flatmap rather than slices. ``plot_map_flatmap`` samples the volume onto the flatmap
surface (``space='FSL'`` uses FSL-MNI152 surfaces, so an MNI152 connectivity map maps
directly, no SUIT reslicing) and draws it with the same cold–hot colormap, symmetric
scale, and threshold as the cortical panel.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from SUITPy import flatmap as suit_flatmap

from ._base import RenderResult, apply_style, load_map, resolve_scale

__all__ = ["plot_map_flatmap", "map_to_flatmap"]


def map_to_flatmap(
    nifti: str | Path | nib.Nifti1Image,
    *,
    space: str = "FSL",
    ignore_zeros: bool = True,
    depths=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
) -> np.ndarray:
    """Sample a volume onto the SUIT flatmap vertices (length ~28935).

    ``space='FSL'`` (alias of SUITPy's ``'MNI'``) selects the FSL-MNI152 white/pial
    cerebellar surfaces, so a standard 2 mm MNI152 connectivity/stat map is sampled in
    its own space without an intermediate SUIT normalization step.
    """
    data = suit_flatmap.vol_to_surf(
        load_map(nifti), space=space, ignore_zeros=ignore_zeros, depths=list(depths)
    )
    data = np.asarray(data, dtype=float)
    if data.ndim > 1:
        data = data[:, 0]
    return data


def plot_map_flatmap(
    nifti: str | Path | nib.Nifti1Image,
    *,
    cmap="lnm_cold_hot",
    vmin: float | None = None,
    vmax: float | None = None,
    threshold: float | None = None,
    space: str = "FSL",
    percentile: float = 98.0,
    round_to: float | None = None,
    borders: bool = True,
    bordercolor: str = "0.35",
    bordersize: float = 1.2,
    underscale: tuple[float, float] = (-0.6, 0.6),
    colorbar: bool = True,
    cbar_label: str | None = None,
    title: str | None = None,
    figure: plt.Figure | None = None,
    axes: plt.Axes | None = None,
    figsize: tuple[float, float] = (3.2, 3.0),
) -> RenderResult:
    """Render ``nifti`` on the SUIT flatmap with the paper cold–hot style.

    Symmetric ``vmin``/``vmax`` (auto from a robust percentile when omitted) and a
    magnitude ``threshold`` matched to the cortical panel; sub-threshold vertices are
    left transparent so the gray cerebellar underlay shows through. Returns a
    :class:`RenderResult`. Pass ``figure``+``axes`` to draw into a composite.
    """
    apply_style()
    sc = resolve_scale(nifti, cmap=cmap, vmin=vmin, vmax=vmax, threshold=threshold,
                       percentile=percentile, round_to=round_to)
    cmap, vmin, vmax, norm = sc.cmap, sc.vmin, sc.vmax, sc.norm

    data = map_to_flatmap(nifti, space=space)

    if axes is None:
        figure = figure or plt.figure(figsize=figsize)
        axes = figure.add_subplot(111)
    else:
        figure = figure or axes.get_figure()
    plt.sca(axes)

    thr = list(sc.threshold_band) if sc.threshold_band is not None else None
    suit_flatmap.plot(
        data.copy(),            # SUITPy thresholds in place; keep our array clean
        overlay_type="func",
        cmap=cmap,
        cscale=[vmin, vmax],
        threshold=thr,
        borders="borders.txt" if borders else None,
        bordercolor=bordercolor,
        bordersize=bordersize,
        underscale=list(underscale),
        new_figure=False,
        colorbar=False,
        backgroundcolor="w",
    )
    axes.set_aspect("equal")
    axes.axis("off")
    if title:
        axes.set_title(title, fontsize=9, pad=2)

    result = RenderResult(figure, [axes], cmap, norm, vmin, vmax, threshold)
    if colorbar:
        _add_flatmap_colorbar(figure, axes, result, cbar_label)
    return result


def _add_flatmap_colorbar(figure, axes, result: RenderResult, label: str | None) -> None:
    pos = axes.get_position()
    cax = figure.add_axes([pos.x0 + pos.width * 0.30, pos.y0 - 0.02,
                           pos.width * 0.40, 0.02])
    cb = figure.colorbar(result.mappable, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_ticks([result.vmin, 0, result.vmax])
    if label:
        cb.set_label(label, fontsize=8)
