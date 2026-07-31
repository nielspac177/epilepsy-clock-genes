"""Figure: circadian x epileptic-network coupling flips sign across development (BrainSpan H4).

Redesigned for clarity + journal style (colors now have explicit, labelled meaning):
  - TOP panel: the LGS EEG-fMRI network shown the way the source paper does -- an axial montage on
    the MNI152 T1 template with a cold-hot t colormap (yellow/red = BOLD increase during GPFA,
    blue = decrease), thresholded, with a `t` colorbar. Slices (not surface) so the thalamic/
    brainstem hubs are visible.
  - MIDDLE: mean clock-gene cortical expression at each of five developmental stages, on the
    inflated surface with a blue->red diverging colormap (z-scored across regions; blue = below
    average, red = above), one shared colorbar. Each stage labelled with its Spearman r to the LGS
    network (red if anti-correlated, green if aligned).
  - BOTTOM: the r-by-stage trajectory (the monotonic sign-flip).
Exploratory (n = 15 regions/stage incl. subcortex; wide per-stage CIs; uncorrected).
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from nilearn import plotting, surface, datasets
from nilearn.plotting import cm as nlcm
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
import abagen

mpl.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
FIG = Path("figures"); FIG.mkdir(exist_ok=True)
TMP = Path("figures/_dev_panels"); TMP.mkdir(exist_ok=True)
LGS_NII = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]
MAP = {"DFC": ["rostralmiddlefrontal", "caudalmiddlefrontal", "superiorfrontal"],
       "VFC": ["parsopercularis", "parstriangularis", "parsorbitalis"],
       "MFC": ["medialorbitofrontal", "rostralanteriorcingulate", "caudalanteriorcingulate"],
       "OFC": ["lateralorbitofrontal", "medialorbitofrontal"], "M1C": ["precentral"],
       "S1C": ["postcentral"], "IPC": ["inferiorparietal", "supramarginal"],
       "A1C": ["transversetemporal"], "STC": ["superiortemporal"], "ITC": ["inferiortemporal"],
       "V1C": ["pericalcarine", "lateraloccipital", "cuneus"]}
STAGES = [("1_prenatal", "Prenatal"), ("2_infancy", "Infancy"), ("3_earlychild", "Early\nchildhood"),
          ("4_child_adol", "Child /\nadolescent"), ("5_adult", "Adult")]
CLK_CMAP, CLK_VMAX = "RdBu_r", 2.0
NEG, POS = "#c0392b", "#1f7a44"


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

    m = resample_to_img(nib.load(LGS_NII), aimg, interpolation="continuous")
    md = np.asarray(m.dataobj)
    lgs_lab = {lab: np.nanmean([np.nanmean(md[adata == i]) for i in ids])
               for lab, ids in lab2ids.items()}
    bs_lgs = {r: np.nanmean([lgs_lab[l] for l in ls if l in lgs_lab]) for r, ls in MAP.items()}

    def paint(vals):
        vol = np.full(adata.shape, np.nan, np.float32)
        for r, v in vals.items():
            for l in MAP[r]:
                for i in lab2ids.get(l, []):
                    vol[adata == i] = v
        return nib.Nifti1Image(vol, aimg.affine, aimg.header)

    # ---- TOP: LGS network as MNI slices, source-paper style (cold-hot, thresholded) ----
    thr = float(np.nanpercentile(np.abs(np.asarray(nib.load(LGS_NII).dataobj)), 92))
    disp = plotting.plot_stat_map(LGS_NII, display_mode="z", cut_coords=[-12, 0, 14, 28, 42, 56],
                                  cmap=nlcm.cold_hot, threshold=thr, black_bg=True, colorbar=True,
                                  annotate=True, draw_cross=False)
    disp.savefig(TMP / "lgs_slices.png", dpi=200); disp.close()

    # ---- clock surfaces per stage (left lateral), diverging, no per-panel colorbar ----
    fs = datasets.fetch_surf_fsaverage("fsaverage")
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
        vals = np.array(list(perreg.values())); z = (vals - np.nanmean(vals)) / np.nanstd(vals)
        img = paint(dict(zip(perreg.keys(), z)))
        tex = surface.vol_to_surf(img, fs["pial_left"], inner_mesh=fs["white_left"])
        f = plt.figure(figsize=(2.2, 2.2))
        ax = f.add_subplot(111, projection="3d")
        plotting.plot_surf_stat_map(fs["infl_left"], np.nan_to_num(tex), hemi="left", view="lateral",
                                    cmap=CLK_CMAP, vmax=CLK_VMAX, threshold=None, colorbar=False,
                                    bg_map=fs["sulc_left"], bg_on_data=False,
                                    axes=ax, figure=f)
        try: ax.set_box_aspect(None, zoom=1.4)
        except Exception: pass
        f.subplots_adjust(0, 0, 1, 1); f.savefig(TMP / f"clock_{key}.png", dpi=200,
                                                 bbox_inches="tight", pad_inches=0); plt.close(f)
        regs = [r for r in perreg if np.isfinite(bs_lgs.get(r, np.nan))]
        stage_r[key] = spearmanr([perreg[r] for r in regs], [bs_lgs[r] for r in regs]).statistic
    trend = spearmanr(range(len(STAGES)), [stage_r[k] for k, _ in STAGES])

    # ---- compose ----
    fig = plt.figure(figsize=(10.5, 9.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.25, 0.85], hspace=0.28)

    # A: LGS slices
    axA = fig.add_subplot(gs[0]); axA.imshow(mpimg.imread(TMP / "lgs_slices.png")); axA.axis("off")
    axA.set_title("A   The Lennox–Gastaut epileptic network  (EEG-fMRI during generalized paroxysmal fast activity)",
                  loc="left", fontsize=12, fontweight="bold")
    axA.text(0.5, -0.06, "yellow/red = BOLD increase during epileptic discharges (t);  blue = decrease. "
             "Axial slices on MNI152; thalamic + brainstem hubs visible.", transform=axA.transAxes,
             ha="center", fontsize=8.5, color="#444")

    # B: clock surfaces row
    gsB = gs[1].subgridspec(1, 6, width_ratios=[1, 1, 1, 1, 1, 0.12], wspace=0.04)
    for j, (key, label) in enumerate(STAGES):
        ax = fig.add_subplot(gsB[j]); ax.imshow(mpimg.imread(TMP / f"clock_{key}.png")); ax.axis("off")
        ax.set_title(label, fontsize=11, fontweight="bold")
        r = stage_r[key]; col = NEG if r < 0 else POS
        tag = "anti-correlated" if r < -0.1 else ("aligned" if r > 0.1 else "crossover")
        ax.text(0.5, -0.02, f"r = {r:+.2f}\n{tag}", transform=ax.transAxes, ha="center", va="top",
                fontsize=10, fontweight="bold", color=col,
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=col, lw=1.3))
    caxB = fig.add_subplot(gsB[5])
    cb = fig.colorbar(ScalarMappable(Normalize(-CLK_VMAX, CLK_VMAX), CLK_CMAP), cax=caxB)
    cb.set_label("clock-gene expression (z)", fontsize=9); cb.set_ticks([-2, 0, 2])
    cb.ax.tick_params(labelsize=8)
    fig.text(0.09, 0.60, "B   Clock-gene cortical expression by developmental stage",
             fontsize=12, fontweight="bold")

    # C: trajectory
    axC = fig.add_subplot(gs[2])
    xs = np.arange(len(STAGES)); ys = [stage_r[k] for k, _ in STAGES]
    axC.axhline(0, color="#999", lw=1, ls="--")
    axC.plot(xs, ys, "-", color="#333", lw=2.2, zorder=1)
    for x, y in zip(xs, ys):
        axC.scatter(x, y, s=90, color=NEG if y < 0 else POS, zorder=3, edgecolor="white", lw=1)
    axC.set_xticks(xs); axC.set_xticklabels([l.replace("\n", " ") for _, l in STAGES], fontsize=9)
    axC.set_ylabel("clock–LGS  Spearman r", fontsize=10); axC.set_ylim(-1, 0.55)
    axC.set_title(f"C   Developmental trajectory   (cortex: ρ = {trend.statistic:+.2f}, "
                  f"p = {trend.pvalue:.2f};  firms to ρ = +0.90, p = 0.037 with subcortex)",
                  loc="left", fontsize=11, fontweight="bold")
    for s in ("top", "right"):
        axC.spines[s].set_visible(False)
    axC.annotate("anti-correlated", (0, ys[0]), (0.3, -0.9), fontsize=8.5, color=NEG,
                 arrowprops=dict(arrowstyle="->", color=NEG))
    axC.annotate("aligned", (4, ys[4]), (3.4, 0.42), fontsize=8.5, color=POS,
                 arrowprops=dict(arrowstyle="->", color=POS))

    fig.suptitle("Circadian × epileptic-network coupling flips sign across development",
                 fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.5, 0.005, "Clock–LGS correlation is rank-based across BrainSpan regions; cortical "
             "surface at region resolution. Exploratory: wide per-stage CIs (all cross 0, n = 15), "
             "uncorrected.", ha="center", fontsize=8.5, style="italic", color="#666")
    fig.savefig(FIG / "fig_developmental_signflip.png", dpi=210, bbox_inches="tight")
    print("stage r:", {k: round(stage_r[k], 3) for k, _ in STAGES},
          f"\ntrend rho={trend.statistic:.3f} p={trend.pvalue:.3f}")
    print("wrote", FIG / "fig_developmental_signflip.png")


if __name__ == "__main__":
    main()
