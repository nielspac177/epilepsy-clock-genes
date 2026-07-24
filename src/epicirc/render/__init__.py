"""LNM-style rendering for MNI volumetric maps (Fox/Horn figure standard).

Render any signed statistical map on the inflated fsaverage cortex **and** the Diedrichsen
SUIT cerebellar flatmap, with a matched cold–hot diverging colormap, symmetric limits, and
threshold — the visual grammar of the tremor (Goede) and gait (Luo) network figures.

    from pipeline.render import combine_surface_flatmap
    combine_surface_flatmap(
        "maps/imbalance_network/1year/specificity_t.nii.gz",
        threshold=2.3, title="Imbalance network (1 year)",
        output="results/figures/fig_imbalance_1year",
    )

Lower-level entry points (:func:`plot_map_surface`, :func:`plot_map_flatmap`) render each
anatomy alone; :func:`combine_surface_flatmap` composes both with one shared scale/colorbar.
See docs/figure_style_notes.md for the paper choices this reproduces.
"""
from __future__ import annotations

from ._base import RenderResult, symmetric_limits
from .cerebellum import map_to_flatmap, plot_map_flatmap
from .colormaps import auto_threshold, lnm_cold_hot, resolve_cmap
from .combine import combine_surface_flatmap
from .surface import plot_map_surface, project_to_fsaverage

__all__ = [
    "plot_map_surface",
    "plot_map_flatmap",
    "combine_surface_flatmap",
    "project_to_fsaverage",
    "map_to_flatmap",
    "lnm_cold_hot",
    "resolve_cmap",
    "auto_threshold",
    "symmetric_limits",
    "RenderResult",
]
