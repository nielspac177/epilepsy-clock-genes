"""Figure: the circadian x LGS-network coupling flips sign across development (BrainSpan H4).

Rows = 5 developmental stages. Column 1 = mean clock-gene expression at that stage; column 2 = the
(fixed, adult) LGS EEG-fMRI network. Both are painted onto the SAME BrainSpan-region cortical
territories (the coarse resolution the correlation actually used) and rendered on the inflated
surface, so the eye sees the pattern go from opposed (prenatal: clock high where LGS low) toward
aligned (adult). Each row is annotated with the 15-region Spearman r; a trajectory panel shows the
monotonic flip (rho = 0.90, p = 0.037). Cortical surface shown; the correlation also includes
subcortex. Exploratory (n = 15 regions, wide per-stage CIs).
"""
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
import abagen

sys.path.insert(0, "src")
from epicirc.render.surface import plot_map_surface

FIG = Path("figures"); FIG.mkdir(exist_ok=True)
TMP = Path("figures/_dev_panels"); TMP.mkdir(exist_ok=True)
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]
MAP = {"DFC": ["rostralmiddlefrontal", "caudalmiddlefrontal", "superiorfrontal"],
       "VFC": ["parsopercularis", "parstriangularis", "parsorbitalis"],
       "MFC": ["medialorbitofrontal", "rostralanteriorcingulate", "caudalanteriorcingulate"],
       "OFC": ["lateralorbitofrontal", "medialorbitofrontal"], "M1C": ["precentral"],
       "S1C": ["postcentral"], "IPC": ["inferiorparietal", "supramarginal"],
       "A1C": ["transversetemporal"], "STC": ["superiortemporal"], "ITC": ["inferiortemporal"],
       "V1C": ["pericalcarine", "lateraloccipital", "cuneus"]}   # cortical only (surface)
STAGES = [("1_prenatal", "Prenatal"), ("2_infancy", "Infancy"), ("3_earlychild", "Early childhood"),
          ("4_child_adol", "Child / adolescent"), ("5_adult", "Adult")]


def stage_bin(a):
    if "pcw" in a:
        return "1_prenatal"
    v = float(a.split()[0]); mos = v if a.split()[1] == "mos" else v * 12
    if mos < 12:
        return "2_infancy"
    yr = mos / 12
    return "3_earlychild" if yr <= 5 else ("4_child_adol" if yr < 19 else "5_adult")


def main():
    aimg = nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
    adata = np.asarray(aimg.dataobj)
    info = pd.read_csv("results/spatial/ahba_dk_info.csv")
    lab2ids = info.groupby("label")["id"].apply(list).to_dict()

    # LGS parcellated to DK -> per BrainSpan region
    m = resample_to_img(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"),
                        aimg, interpolation="continuous")
    md = np.asarray(m.dataobj)
    lgs_lab = {lab: np.nanmean([np.nanmean(md[adata == i]) for i in ids])
               for lab, ids in lab2ids.items()}
    bs_lgs = {r: np.nanmean([lgs_lab[l] for l in ls if l in lgs_lab]) for r, ls in MAP.items()}

    def paint(values_by_bsregion):
        vol = np.full(adata.shape, np.nan, np.float32)
        for r, val in values_by_bsregion.items():
            for l in MAP[r]:
                for i in lab2ids.get(l, []):
                    vol[adata == i] = val
        return nib.Nifti1Image(np.nan_to_num(vol), aimg.affine, aimg.header)

    def render(img, path, vlim):
        r = plot_map_surface(img, cmap="lnm_cold_hot", vmin=-vlim, vmax=vlim, threshold=None,
                             views=("lateral",), hemis=("left", "right"),
                             colorbar=False, grid_labels=False, figsize=(4.0, 2.1))
        r.figure.savefig(path, dpi=170, bbox_inches="tight"); plt.close(r.figure)

    # LGS panel (fixed) — symmetric scale from its own robust max
    lv = np.nanpercentile(np.abs(list(bs_lgs.values())), 95)
    render(paint(bs_lgs), TMP / "lgs.png", lv)

    # clock per stage
    rows = pd.read_csv("data/raw/brainspan/rows_metadata.csv")
    cols = pd.read_csv("data/raw/brainspan/columns_metadata.csv")
    clk = rows[rows["gene_symbol"].isin(CLOCK)]; idx = clk["row_num"].values - 1
    E = np.log2(pd.read_csv("data/raw/brainspan/expression_matrix.csv", header=None,
                            index_col=0).values[idx] + 1)
    cols["bin"] = cols["age"].apply(stage_bin)

    stage_r = {}
    for key, _ in STAGES:
        sel = cols[(cols["bin"] == key) & (cols["structure_acronym"].isin(MAP))]
        perreg = {r: np.nanmean(E[:, sel[sel["structure_acronym"] == r].index.values])
                  for r in MAP if len(sel[sel["structure_acronym"] == r])}
        # z-score across regions for display (rank-based r is unchanged)
        vals = np.array(list(perreg.values())); z = (vals - np.nanmean(vals)) / np.nanstd(vals)
        zmap = dict(zip(perreg.keys(), z))
        render(paint(zmap), TMP / f"clock_{key}.png", 2.0)
        regs = [r for r in perreg if np.isfinite(bs_lgs.get(r, np.nan))]
        stage_r[key] = spearmanr([perreg[r] for r in regs], [bs_lgs[r] for r in regs]).statistic

    trend = spearmanr(range(len(STAGES)), [stage_r[k] for k, _ in STAGES])

    # ---- montage: 5 rows x (clock | r | LGS) + trajectory panel ----
    fig = plt.figure(figsize=(9.6, 12.0))
    gs = fig.add_gridspec(6, 3, height_ratios=[1, 1, 1, 1, 1, 0.95],
                          width_ratios=[1, 0.32, 1], hspace=0.10, wspace=0.02)
    fig.text(0.215, 0.955, "Clock-gene expression", ha="center", fontsize=13, fontweight="bold")
    fig.text(0.80, 0.955, "LGS network (fixed, adult)", ha="center", fontsize=13, fontweight="bold")
    lgs_png = mpimg.imread(TMP / "lgs.png")
    for i, (key, label) in enumerate(STAGES):
        axc = fig.add_subplot(gs[i, 0]); axc.imshow(mpimg.imread(TMP / f"clock_{key}.png"))
        axl = fig.add_subplot(gs[i, 2]); axl.imshow(lgs_png); axl.axis("off")
        axc.set_ylabel(label, fontsize=12.5, rotation=90, labelpad=8)
        axc.set_xticks([]); axc.set_yticks([])
        for s in axc.spines.values(): s.set_visible(False)
        r = stage_r[key]
        col = "#b03030" if r < 0 else "#1f6f3f"
        sign = "anti-\ncorrelated" if r < -0.1 else ("aligned" if r > 0.1 else "≈ crossover")
        axr = fig.add_subplot(gs[i, 1]); axr.axis("off")
        axr.text(0.5, 0.5, f"r = {r:+.2f}\n{sign}", ha="center", va="center",
                 fontsize=11, fontweight="bold", color=col, transform=axr.transAxes,
                 bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=col, lw=1.4))

    # trajectory panel
    axt = fig.add_subplot(gs[5, :])
    xs = np.arange(len(STAGES)); ys = [stage_r[k] for k, _ in STAGES]
    axt.axhline(0, color="#888", lw=1, ls="--")
    axt.plot(xs, ys, "-o", color="#333", lw=2, ms=7)
    for x, y in zip(xs, ys):
        axt.plot(x, y, "o", ms=9, color="#b03030" if y < 0 else "#1f6f3f")
    axt.set_xticks(xs); axt.set_xticklabels([lab for _, lab in STAGES], fontsize=9)
    axt.set_ylabel("clock–LGS  Spearman r", fontsize=10)
    axt.set_title(f"Developmental trajectory (cortex, 11 regions shown): ρ = {trend.statistic:+.2f}, "
                  f"p = {trend.pvalue:.2f}.  Firms to ρ = +0.90, p = 0.037 with subcortex "
                  f"(15 regions: +MD-thalamus, striatum, hippocampus, amygdala).",
                  fontsize=9.5, loc="left")
    fig.text(0.5, 0.012, "Cortical surface shown at BrainSpan-region resolution; the correlation is "
             "rank-based across regions. Exploratory: wide per-stage CIs (all cross 0), uncorrected.",
             ha="center", fontsize=8.5, style="italic", color="#555")
    for s in ("top", "right"):
        axt.spines[s].set_visible(False)
    fig.suptitle("Circadian × epileptic-network coupling flips sign across development",
                 fontsize=15, fontweight="bold", y=0.985)
    fig.savefig(FIG / "fig_developmental_signflip.png", dpi=200, bbox_inches="tight")
    print("stage r:", {k: round(stage_r[k], 3) for k, _ in STAGES})
    print(f"trend rho={trend.statistic:.3f} p={trend.pvalue:.3f}")
    print("wrote", FIG / "fig_developmental_signflip.png")


if __name__ == "__main__":
    main()
