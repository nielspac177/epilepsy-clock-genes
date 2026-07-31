"""Phase 2.5: does clock-gene cortical expression align with ANY canonical brain map?

The clock<->LGS spatial bridge is null. This asks the constructive complement: paint mean clock-gene
AHBA expression onto the fsLR-32k surface and correlate it (spin nulls, FDR) against the neuromaps
battery (fsLR-32k native + MNI152->32k maps). If circadian expression aligns with the principal
gradient / a receptor / metabolism, that is a positive finding; if it aligns with nothing, that
reinforces that circadian genetics has no macroscale cortical-topographic signature.
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from neuromaps import transforms, nulls, stats
from neuromaps.datasets import available_annotations, fetch_annotation
from neuromaps.images import load_data

OUT = Path("results/neuromaps_battery"); OUT.mkdir(parents=True, exist_ok=True)
N_PERM = 1000
SEED = 20260731
CLOCK = ["ARNTL", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1",
         "NR1D2", "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
         "BHLHE40", "BHLHE41", "CIART"]


def clock_surface_32k():
    import abagen
    atlas = abagen.fetch_desikan_killiany(native=False)
    aimg = nib.load(str(atlas["image"])); adata = np.asarray(aimg.dataobj)
    expr = pd.read_csv("results/spatial/ahba_dk_expression.csv", index_col=0)
    clk = [g for g in CLOCK if g in expr.columns]
    Z = (expr[clk] - expr[clk].mean()) / expr[clk].std()
    val = Z.mean(1)
    vol = np.zeros(adata.shape, dtype=np.float32)
    for lab in val.index:
        vol[adata == lab] = val.loc[lab]
    img = nib.Nifti1Image(vol, aimg.affine, aimg.header)
    return load_data(transforms.mni152_to_fslr(img, "32k"))


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def main():
    clock = clock_surface_32k()
    print("clock map on fsLR-32k ready", flush=True)
    rot = nulls.alexander_bloch(clock, atlas="fsLR", density="32k", n_perm=N_PERM, seed=SEED)

    # maps reachable to fsLR-32k without workbench: fsLR/32k native or MNI152
    best = {}
    for src, desc, space, den in available_annotations():
        if (space == "fsLR" and den == "32k") or space == "MNI152":
            best.setdefault((src, desc), (space, den))
    print(f"{len(best)} maps testable at fsLR-32k", flush=True)

    rows, log = [], []
    for i, ((src, desc), (space, den)) in enumerate(sorted(best.items())):
        try:
            a = fetch_annotation(source=src, desc=desc, space=space, den=den, verbose=0)
            g = load_data(transforms.mni152_to_fslr(a, "32k")) if space == "MNI152" else load_data(a)
            r, p = stats.compare_images(clock, g, nulls=rot, metric="spearmanr",
                                        ignore_zero=True, nan_policy="omit")
            rows.append([f"{src}:{desc}", r, p])
            print(f"  [{i+1}/{len(best)}] {src}:{desc:16s} r={r:+.3f} p={p:.4f}", flush=True)
        except Exception as e:
            log.append(f"{src}:{desc} FAILED: {str(e)[:120]}")

    ps = [r[2] for r in rows]
    q = bh_fdr(ps)
    for row, qq in zip(rows, q):
        row.append(qq)
    with (OUT / "clock_battery.tsv").open("w") as fh:
        fh.write("map\tspearman_r\tp_spin\tp_fdr\n")
        for m, r, p, qq in sorted(rows, key=lambda x: x[2]):
            fh.write(f"{m}\t{r:.4f}\t{p:.4f}\t{qq:.4f}\n")
    nsig = sum(1 for r in rows if r[3] < 0.05)
    print(f"\n{nsig}/{len(rows)} maps FDR<0.05 with clock expression", flush=True)
    print("CLOCK_BATTERY_DONE", flush=True)


if __name__ == "__main__":
    main()
