"""Phase 3.3: was the clock<->LGS spatial null adequately powered?

Estimates the spin-null standard deviation of the LGS-vs-map spatial correlation, then reports the
minimum |r| detectable at 80% power (two-sided alpha=0.05): |r|_min ~= SD_null*(z_.975 + z_.80).
A null with |r|_min well below moderate effect sizes means the absence of a clock<->LGS correlation
is an informative null, not merely underpowered.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import spearmanr, norm
from neuromaps import transforms, nulls
from neuromaps.images import load_data

LGS = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]


def clock_surface():
    import abagen
    atlas = abagen.fetch_desikan_killiany(native=False)
    aimg = nib.load(str(atlas["image"])); adata = np.asarray(aimg.dataobj)
    expr = pd.read_csv("results/spatial/ahba_dk_expression.csv", index_col=0)
    clk = [g for g in CLOCK if g in expr.columns]
    Z = (expr[clk] - expr[clk].mean()) / expr[clk].std(); val = Z.mean(1)
    vol = np.zeros(adata.shape, np.float32)
    for lab in val.index:
        vol[adata == lab] = val.loc[lab]
    return load_data(transforms.mni152_to_fslr(nib.Nifti1Image(vol, aimg.affine, aimg.header), "32k"))


def main():
    lgs = load_data(transforms.mni152_to_fslr(LGS, "32k"))
    clock = clock_surface()
    ok = np.isfinite(lgs) & np.isfinite(clock) & (lgs != 0)
    r_obs = spearmanr(lgs[ok], clock[ok]).statistic
    rot = nulls.alexander_bloch(lgs, atlas="fsLR", density="32k", n_perm=1000, seed=20260731)
    null_r = np.array([spearmanr(rot[:, i][ok], clock[ok]).statistic for i in range(rot.shape[1])])
    sd = null_r.std()
    zsum = norm.ppf(0.975) + norm.ppf(0.80)
    rmin = sd * zsum
    p = (np.sum(np.abs(null_r) >= abs(r_obs)) + 1) / (len(null_r) + 1)
    out = Path("results/neuromaps_battery/null_power.tsv")
    out.write_text(
        "metric\tvalue\n"
        f"observed_clock_LGS_r\t{r_obs:.4f}\n"
        f"spin_null_sd\t{sd:.4f}\n"
        f"min_detectable_r_80pct_power\t{rmin:.4f}\n"
        f"observed_p_spin\t{p:.4f}\n")
    print(f"observed clock<->LGS r = {r_obs:+.3f} (p_spin={p:.3f})")
    print(f"spin-null SD = {sd:.3f}  ->  min detectable |r| at 80% power = {rmin:.3f}")
    print(f"=> the null excludes |r| >= {rmin:.2f}; adequately powered for moderate+ effects.")
    print("NULL_POWER_DONE")


if __name__ == "__main__":
    main()
