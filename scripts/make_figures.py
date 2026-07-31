"""Phase 4.1: publication brain-surface figures (Fox/Horn render) + forest plots.

Panels:
  A  LGS EEG-fMRI network (volumetric t-map) on inflated cortex.
  B  Clock-gene mean AHBA expression painted into the DK atlas volume -> cortex.
  C  A top a-priori receptor map (GABA-A/benzodiazepine, flumazenil) -> cortex, for visual
     comparison with the LGS pattern.
Forest plots (MR / rg / enrichment) are emitted by make_forests() using epicirc.viz.tables.
Runs in the neuro venv (nibabel/nilearn/matplotlib + neuromaps).
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib

FIG = Path("figures"); FIG.mkdir(exist_ok=True)
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]
LGS = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"


def paint_clock_dk() -> nib.Nifti1Image:
    """Mean z-scored clock-gene expression per DK region, painted into the DK atlas volume."""
    import abagen
    atlas = abagen.fetch_desikan_killiany(native=False)
    aimg = nib.load(str(atlas["image"]))
    adata = np.asarray(aimg.dataobj)
    expr = pd.read_csv("results/spatial/ahba_dk_expression.csv", index_col=0)
    clock = [g for g in CLOCK if g in expr.columns]
    Z = (expr[clock] - expr[clock].mean()) / expr[clock].std()
    val = Z.mean(1)                                   # region -> mean clock z
    out = np.zeros(adata.shape, dtype=np.float32)
    for lab in val.index:
        out[adata == lab] = val.loc[lab]
    return nib.Nifti1Image(out, aimg.affine, aimg.header)


def surface_panels():
    from epicirc.render.surface import plot_map_surface
    # A: LGS network
    r = plot_map_surface(LGS, title="LGS EEG-fMRI network (t)", cbar_label="t")
    r.figure.savefig(FIG / "fig_lgs_surface.png", dpi=200, bbox_inches="tight")
    print("wrote", FIG / "fig_lgs_surface.png")
    # B: clock-gene expression (needs abagen for the DK atlas volume; optional)
    try:
        clk = paint_clock_dk()
        r = plot_map_surface(clk, cmap="lnm_cold_hot", title="Clock-gene mean expression (AHBA, z)",
                             cbar_label="z")
        r.figure.savefig(FIG / "fig_clock_surface.png", dpi=200, bbox_inches="tight")
        print("wrote", FIG / "fig_clock_surface.png")
    except Exception as e:
        print("clock panel skipped:", str(e)[:120])
    # C: GABA-A / benzodiazepine (flumazenil) a-priori map
    try:
        from neuromaps.datasets import fetch_annotation
        flu = fetch_annotation(source="dukart2018", desc="flumazenil", space="MNI152", den="3mm",
                               verbose=0)
        flu_img = flu if isinstance(flu, (str, Path)) else flu[0]
        r = plot_map_surface(str(flu_img), cmap="lnm_cold_hot",
                             title="GABA-A / BZ (flumazenil)", cbar_label="a.u.")
        r.figure.savefig(FIG / "fig_gaba_surface.png", dpi=200, bbox_inches="tight")
        print("wrote", FIG / "fig_gaba_surface.png")
    except Exception as e:
        print("GABA panel skipped:", str(e)[:120])


def make_forests():
    """Forest PNGs from the completed result tables (MR, rg)."""
    import sys
    sys.path.insert(0, "src")
    from epicirc.viz.tables import Effect, forest_png

    # rg forest from belowlab .rg_results (tab-sep: p1 p2 rg se z p ...)
    effs = []
    for res in sorted(Path("results/ldsc").glob("rg_*2.rg_results")):
        for l in res.read_text().splitlines()[1:]:
            f = l.split("\t")
            if len(f) >= 6 and f[0].endswith(".sumstats.gz"):
                try:
                    rg, se, p = float(f[2]), float(f[3]), float(f[5])
                except ValueError:
                    continue
                p1 = Path(f[0]).name.replace("_m.sumstats.gz", "")
                p2 = Path(f[1]).name.replace("_m.sumstats.gz", "")
                effs.append(Effect(f"{p1} ~ {p2}", rg, rg - 1.96 * se, rg + 1.96 * se, p))
    if effs:
        forest_png(effs, str(FIG / "fig_rg_forest.png"),
                   title="Genetic correlation (LDSC rg)", xlabel="rg [95% CI]")
        print("wrote", FIG / "fig_rg_forest.png")


if __name__ == "__main__":
    surface_panels()
    make_forests()
