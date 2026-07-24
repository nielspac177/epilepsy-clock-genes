"""Phase D: does clock-gene expression spatially co-localize with the LGS EEG-fMRI network?

Parcellates the unthresholded LGS t-map into the DK atlas, then correlates it (across 83 regions)
with each clock gene's AHBA expression and with the set-level clock map. Two nulls are enforced:
  (1) SPATIAL — brainsmash variogram surrogates of the LGS map (preserves spatial autocorrelation).
  (2) RANDOM-GENE — 10,000 random gene sets of the same size (specificity vs any brain-expressed set).
Only signals surviving BOTH nulls are meaningful. No causal language is drawn from a spatial r.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
import abagen
from brainsmash.mapgen.base import Base

OUT = Path("results/spatial"); OUT.mkdir(parents=True, exist_ok=True)
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]
LGS_MAP = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"  # unthresholded t-scores


def parcellate(map_path, atlas_img):
    m = resample_to_img(nib.load(map_path), atlas_img, interpolation="continuous")
    md = np.asarray(m.dataobj); ad = np.asarray(atlas_img.dataobj)
    labels = [int(x) for x in np.unique(ad) if x != 0]
    return {lab: float(np.nanmean(md[ad == lab])) for lab in labels}


def main():
    atlas = abagen.fetch_desikan_killiany(native=False)
    atlas_img = nib.load(str(atlas["image"]))
    expr = pd.read_csv(OUT / "ahba_dk_expression.csv", index_col=0)          # 83 x genes
    cent = pd.read_csv(OUT / "ahba_dk_centroids.csv").set_index("id")

    lgs = parcellate(LGS_MAP, atlas_img)                                     # {label: mean t}
    lgs = pd.Series(lgs).reindex(expr.index)                                 # align to expr regions

    # regions usable in both, and present clock genes
    ok = expr.notna().all(1) & lgs.notna()
    reg = expr.index[ok]
    lgs_v = lgs.loc[reg].to_numpy()
    E = expr.loc[reg]
    clock = [g for g in CLOCK if g in E.columns]
    print(f"regions used: {len(reg)}/83 | clock genes: {len(clock)}")

    # distance matrix (mm) among used regions -> brainsmash variogram surrogates of LGS
    xyz = cent.loc[reg, ["x", "y", "z"]].to_numpy()
    D = cdist(xyz, xyz)
    np.random.seed(20260724)   # seed the brainsmash null for reproducible p_spatial (audit fix)
    gen = Base(x=lgs_v, D=D, resample=True)
    surr = gen(n=1000)                                                       # 1000 x nregion

    def spatial_p(vec):
        r = spearmanr(vec, lgs_v).statistic
        rn = np.array([spearmanr(vec, s).statistic for s in surr])
        return r, (np.sum(np.abs(rn) >= abs(r)) + 1) / (len(rn) + 1)

    # set-level clock map = mean of z-scored clock-gene expression
    Z = (E - E.mean()) / E.std()
    clock_map = Z[clock].mean(1).to_numpy()
    r_set, p_spatial = spatial_p(clock_map)

    # random-gene null for the SET (same size, real LGS map)
    rng = np.random.default_rng(20260723)
    allg = E.columns.to_numpy()
    r_rand = []
    for _ in range(10000):
        pick = rng.choice(allg, len(clock), replace=False)
        cm = Z[pick].mean(1).to_numpy()
        r_rand.append(spearmanr(cm, lgs_v).statistic)
    r_rand = np.array(r_rand)
    p_gene = (np.sum(np.abs(r_rand) >= abs(r_set)) + 1) / (len(r_rand) + 1)

    # per-gene (spatial null only)
    rows = [("CLOCK_SET", r_set, p_spatial, p_gene)]
    for g in clock:
        r, p = spatial_p(E[g].to_numpy())
        rows.append((g, r, p, np.nan))
    df = pd.DataFrame(rows, columns=["gene", "spearman_r", "p_spatial", "p_randomgene"])
    df.to_csv(OUT / "lgs_clock_spatial.tsv", sep="\t", index=False)

    print(f"\nSET-level: r={r_set:.3f}  p_spatial={p_spatial:.3f}  p_randomgene={p_gene:.3f}")
    print(df.sort_values("p_spatial").to_string(index=False))


if __name__ == "__main__":
    main()
