"""Phase 2.6: parcellation robustness of the LGS receptor findings in Schaefer-400.

The headline LGS x receptor correlations were vertex-wise on fsLR-32k. Here they are recomputed in a
completely different parcellation (Schaefer-400, volumetric MNI) with a brainsmash variogram null on
the parcel centroids (spatial-autocorrelation-preserving) -- a different atlas AND a different null
family. Receptor claims that reproduce here are not DK/vertex or spin-test artifacts.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import nibabel as nib
from nilearn import datasets as nldata
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from brainsmash.mapgen.base import Base
from neuromaps.datasets import fetch_annotation

OUT = Path("results/neuromaps_battery"); OUT.mkdir(parents=True, exist_ok=True)
LGS = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"
SEED = 20260731

# MNI-available maps to re-test (source, desc, space, den, label)
MAPS = [
    ("dubois2015", "abp688", "MNI152", "1mm", "mGluR5 (ABP688)"),
    ("dukart2018", "flumazenil", "MNI152", "3mm", "GABAa/BZ (flumazenil)"),
    ("lukow2022", "ro154513", "MNI152", "2mm", "GABAa/BZ (Ro15-4513)"),
    ("castrillon2023", "cmrglc", "MNI152", "3mm", "glucose metabolism"),
]


def parcellate(img, atlas_img, adata, labels):
    m = resample_to_img(img, atlas_img, interpolation="continuous", force_resample=True,
                        copy_header=True)
    md = np.asarray(m.dataobj)
    return np.array([np.nanmean(md[adata == lab]) for lab in labels])


def main():
    sch = nldata.fetch_atlas_schaefer_2018(n_rois=400, resolution_mm=2)
    atlas_img = nib.load(sch["maps"]); adata = np.asarray(atlas_img.dataobj)
    labels = [int(x) for x in np.unique(adata) if x != 0]
    print(f"Schaefer-400: {len(labels)} parcels", flush=True)

    # parcel centroids (voxel->mm) for the variogram null
    aff = atlas_img.affine
    cents = []
    for lab in labels:
        vox = np.argwhere(adata == lab).mean(0)
        cents.append(nib.affines.apply_affine(aff, vox))
    D = cdist(np.array(cents), np.array(cents))

    lgs_v = parcellate(nib.load(LGS), atlas_img, adata, labels)
    ok = np.isfinite(lgs_v)
    np.random.seed(SEED)
    surr = Base(x=lgs_v[ok], D=D[np.ix_(ok, ok)], resample=True)(n=1000)

    rows = []
    for src, desc, space, den, lab in MAPS:
        a = fetch_annotation(source=src, desc=desc, space=space, den=den, verbose=0)
        aimg = nib.load(a if isinstance(a, str) else a[0])
        mv = parcellate(aimg, atlas_img, adata, labels)
        both = ok & np.isfinite(mv)
        r = spearmanr(lgs_v[both], mv[both]).statistic
        rn = np.array([spearmanr(s, mv[ok]).statistic for s in surr])
        p = (np.sum(np.abs(rn) >= abs(r)) + 1) / (len(rn) + 1)
        rows.append((lab, r, p))
        print(f"  {lab:24s} r={r:+.3f}  p_variogram={'<0.001' if p<=0.001 else f'{p:.3f}'}", flush=True)

    with (OUT / "schaefer400_robust.tsv").open("w") as fh:
        fh.write("map\tspearman_r_schaefer400\tp_variogram\n")
        for lab, r, p in rows:
            fh.write(f"{lab}\t{r:.4f}\t{p:.4f}\n")
    print("SCHAEFER_DONE", flush=True)


if __name__ == "__main__":
    main()
