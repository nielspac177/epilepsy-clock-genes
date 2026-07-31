"""Phase 2.4: is the LGS-receptor spatial overlap SPECIFIC, or just the shared principal gradient?

The battery shows LGS correlates with 36/68 maps and most strongly with the principal FC gradient
(margulies fcgradient01, r=0.70). mGluR5/GABA-A receptor maps also track that gradient, so a raw
spin-significant correlation does NOT prove a receptor-specific relationship. Here we compute the
PARTIAL spatial correlation of LGS with each a-priori/highlight map CONTROLLING for the principal
gradient (and, in a second model, for myelin too), with a spin-test null on LGS. A receptor claim
survives only if the partial correlation remains spin-significant after removing the gradient.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
from neuromaps import transforms, nulls
from neuromaps.datasets import fetch_annotation
from neuromaps.images import load_data

OUT = Path("results/neuromaps_battery"); OUT.mkdir(parents=True, exist_ok=True)
LGS = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"
N_PERM = 1000
SEED = 20260731

# maps to test (all brought to fsLR-32k); (source, desc, space, den, label)
MAPS = [
    ("margulies2016", "fcgradient01", "fsLR", "32k", "principal gradient"),
    ("hcps1200", "myelinmap", "fsLR", "32k", "myelin (T1w/T2w)"),
    ("dubois2015", "abp688", "MNI152", "1mm", "mGluR5 (ABP688)"),
    ("lukow2022", "ro154513", "MNI152", "2mm", "GABAa/BZ (Ro15-4513)"),
    ("dukart2018", "flumazenil", "MNI152", "3mm", "GABAa/BZ (flumazenil)"),
    ("castrillon2023", "cmrglc", "MNI152", "3mm", "glucose metabolism"),
]


def to32k(src, desc, space, den):
    a = fetch_annotation(source=src, desc=desc, space=space, den=den, verbose=0)
    if space == "MNI152":
        return load_data(transforms.mni152_to_fslr(a, "32k"))
    if den != "32k":
        return load_data(transforms.fslr_to_fslr(a, "32k"))
    return load_data(a)


def _resid_rank(y, X):
    """Residuals of rank(y) on [1, rank(X columns)] (Spearman partialling)."""
    yr = rankdata(y)
    Xr = np.column_stack([np.ones(len(yr))] + [rankdata(X[:, j]) for j in range(X.shape[1])])
    beta, *_ = np.linalg.lstsq(Xr, yr, rcond=None)
    return yr - Xr @ beta


def partial_spearman(a, b, covars):
    """partial Spearman(a, b | covars) via residual correlation of rank-residuals."""
    ra = _resid_rank(a, covars)
    rb = _resid_rank(b, covars)
    return float(np.corrcoef(ra, rb)[0, 1])


def main():
    lgs = load_data(transforms.mni152_to_fslr(LGS, "32k"))
    maps = {}
    for src, desc, space, den, lab in MAPS:
        maps[lab] = to32k(src, desc, space, den)
        print(f"loaded {lab}", flush=True)
    grad = maps["principal gradient"]
    mye = maps["myelin (T1w/T2w)"]

    rot = nulls.alexander_bloch(lgs, atlas="fsLR", density="32k", n_perm=N_PERM, seed=SEED)
    # valid vertices (finite across everything)
    stack = np.vstack([lgs, grad, mye] + [maps[l] for _, _, _, _, l in MAPS])
    ok = np.all(np.isfinite(stack), axis=0) & (lgs != 0)

    rows = []
    receptor_maps = [(l) for _, _, _, _, l in MAPS if l not in ("principal gradient",)]
    for lab in receptor_maps:
        m = maps[lab][ok]
        g = grad[ok]; my = mye[ok]
        # raw, +gradient, +gradient+myelin
        raw = partial_spearman(lgs[ok], m, np.empty((ok.sum(), 0)))
        pg = partial_spearman(lgs[ok], m, g[:, None])
        pgm = partial_spearman(lgs[ok], m, np.column_stack([g, my]))
        # spin null for the gradient-controlled partial r
        null = []
        for perm in rot.T:
            ls = perm[ok] if perm.shape[0] == lgs.shape[0] else perm  # rot rows are vertices
            null.append(partial_spearman(ls, m, g[:, None]))
        null = np.array(null)
        p_spin = (np.sum(np.abs(null) >= abs(pg)) + 1) / (len(null) + 1)
        rows.append((lab, raw, pg, pgm, p_spin))
        print(f"  {lab:24s} raw={raw:+.3f}  |grad={pg:+.3f}  |grad+mye={pgm:+.3f}  "
              f"p_spin(|grad)={'<0.001' if p_spin<=0.001 else f'{p_spin:.3f}'}", flush=True)

    with (OUT / "lgs_partial_gradient.tsv").open("w") as fh:
        fh.write("map\tr_raw\tr_partial_grad\tr_partial_grad_myelin\tp_spin_partial_grad\n")
        for lab, raw, pg, pgm, p in rows:
            fh.write(f"{lab}\t{raw:.4f}\t{pg:.4f}\t{pgm:.4f}\t{p:.4f}\n")
    print("wrote", OUT / "lgs_partial_gradient.tsv")
    print("DOMINANCE_DONE")


if __name__ == "__main__":
    main()
