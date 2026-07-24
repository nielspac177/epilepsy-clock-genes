"""Phase E: does the LGS EEG-fMRI network track a developmental/embryological cortical gradient?

Projects the volumetric LGS t-map to the fsLR-32k surface (neuromaps registration fusion) and
correlates it with canonical developmental/hierarchy gradients, each with a spin-test
(spatial-autocorrelation-preserving) null and FDR across maps. Cortical surface only, so the
subcortical/thalamic component of the LGS network is not captured here (a stated limitation).

Gradients:
  sydnor2021 SAaxis        sensorimotor->association axis (the neurodevelopmental hierarchy)
  hill2010 devexp/evoexp   postnatal developmental & evolutionary cortical expansion
  margulies2016 grad01     principal functional connectivity gradient (unimodal->transmodal)
  abagen genepc1           first principal component of cortical gene expression
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
from neuromaps import transforms, datasets, nulls, stats
from neuromaps.images import load_data


def bh_fdr(pvals):
    p = np.asarray(pvals, float); n = len(p)
    order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out

OUT = Path("results/spatial"); OUT.mkdir(parents=True, exist_ok=True)
LGS = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"

GRADS = [
    ("sydnor2021", "SAaxis", "fsLR", "32k"),
    ("margulies2016", "fcgradient01", "fsLR", "32k"),
    ("hill2010", "devexp", "fsLR", "164k"),
    ("hill2010", "evoexp", "fsLR", "164k"),
    ("abagen", "genepc1", "fsaverage", "10k"),
]


def to_fslr32k(data, space, den):
    if space == "fsLR" and den == "32k":
        return data
    if space == "fsLR":
        return transforms.fslr_to_fslr(data, "32k")
    if space == "fsaverage":
        return transforms.fsaverage_to_fslr(data, "32k")
    raise ValueError(space)


def main():
    lgs = transforms.mni152_to_fslr(LGS, "32k")           # (L,R) gifti on fsLR-32k
    lgs_data = load_data(lgs)

    # spin nulls of the LGS map (generated once, reused for every gradient)
    rot = nulls.alexander_bloch(lgs_data, atlas="fsLR", density="32k",
                                n_perm=1000, seed=20260723)

    rows = []
    for src, desc, space, den in GRADS:
        try:
            annot = datasets.fetch_annotation(source=src, desc=desc, space=space, den=den)
            g = to_fslr32k(annot, space, den)
            g_data = load_data(g)
            r, p = stats.compare_images(lgs_data, g_data, nulls=rot, metric="spearmanr",
                                        ignore_zero=True, nan_policy="omit")
            rows.append([f"{src}:{desc}", r, p])
            print(f"  {src}:{desc:14s} r={r:+.3f}  p_spin={p:.4f}")
        except Exception as e:
            print(f"  {src}:{desc} FAILED: {str(e)[:120]}")
            rows.append([f"{src}:{desc}", np.nan, np.nan])

    ps = [r[2] for r in rows if np.isfinite(r[2])]
    fdr = bh_fdr(ps) if ps else []
    it = iter(fdr)
    with open(OUT / "lgs_gradients.tsv", "w") as fh:
        fh.write("gradient\tspearman_r\tp_spin\tp_fdr\n")
        for name, r, p in rows:
            q = next(it) if np.isfinite(p) else np.nan
            fh.write(f"{name}\t{r:.4f}\t{p:.4f}\t{q:.4f}\n")
    print(f"\nwrote {OUT/'lgs_gradients.tsv'}")


if __name__ == "__main__":
    main()
