"""Phase 2.1/2.2: spatial fingerprint of the LGS EEG-fMRI network across the full neuromaps battery.

Robust design: for each annotation, bring the LGS t-map INTO the annotation's native space by
registration fusion (mni152_to_fslr / mni152_to_fsaverage) -- NO Connectome Workbench needed for
MNI152/fsLR-32k/fsLR-164k/fsaverage maps. fsLR-4k (MEG) and civet still need wb_command; those are
attempted only if wb_command is on PATH, else logged. One spin-null set per (space,density) is
cached and reused. A-priori epilepsy-relevant maps (GABA-A/BZ flumazenil & Ro15-4513, mGluR5
abp688, MEG delta/theta) are tagged so they form a small pre-registered primary family separate
from the exploratory battery, with FDR applied within each family.
"""
import warnings; warnings.filterwarnings("ignore")
import shutil
import sys
from pathlib import Path

import numpy as np
from neuromaps import transforms, nulls, stats
from neuromaps.datasets import available_annotations, fetch_annotation
from neuromaps.images import load_data

OUT = Path("results/neuromaps_battery"); OUT.mkdir(parents=True, exist_ok=True)
LGS = "EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"
N_PERM = 1000
SEED = 20260730
HAVE_WB = shutil.which("wb_command") is not None

APRIORI = {
    ("dukart2018", "flumazenil"): "GABAa/BZ (flumazenil)",
    ("norgaard2021", "flumazenil"): "GABAa/BZ (flumazenil)",
    ("lukow2022", "ro154513"): "GABAa/BZ (Ro15-4513)",
    ("dubois2015", "abp688"): "mGluR5 (ABP688)",
    ("rosaneto", "abp688"): "mGluR5 (ABP688)",
    ("smart2019", "abp688"): "mGluR5 (ABP688)",
    ("hcps1200", "megdelta"): "MEG delta power",
    ("hcps1200", "megtheta"): "MEG theta power",
}

# preference order of (space,density) per (source,desc): workbench-free first
def space_rank(space, den):
    order = {("MNI152", "*"): 0, ("fsLR", "32k"): 1, ("fsLR", "164k"): 2,
             ("fsaverage", "164k"): 3, ("fsaverage", "10k"): 4}
    if space == "MNI152":
        return 0
    r = order.get((space, den))
    if r is not None:
        return r
    if space == "fsLR" and den == "4k":
        return 10 if HAVE_WB else 99      # MEG: only with workbench
    return 98                              # civet / others: skip


_lgs_cache, _null_cache = {}, {}


def lgs_in(space, den):
    """LGS projected into (space, den); MNI152 target -> fsLR 32k."""
    key = (space, den)
    if key in _lgs_cache:
        return _lgs_cache[key]
    if space == "MNI152":
        data = load_data(transforms.mni152_to_fslr(LGS, "32k")); tgt = ("fsLR", "32k")
    elif space == "fsLR" and den == "4k":
        data = load_data(transforms.fslr_to_fslr(transforms.mni152_to_fslr(LGS, "32k"), "4k"))
        tgt = ("fsLR", "4k")
    elif space == "fsLR":
        data = load_data(transforms.mni152_to_fslr(LGS, den)); tgt = ("fsLR", den)
    elif space == "fsaverage":
        data = load_data(transforms.mni152_to_fsaverage(LGS, den)); tgt = ("fsaverage", den)
    else:
        raise ValueError(space)
    _lgs_cache[key] = (data, tgt)
    return _lgs_cache[key]


def spin_null(tgt):
    if tgt in _null_cache:
        return _null_cache[tgt]
    lgs_data = next(v[0] for v in _lgs_cache.values() if v[1] == tgt)
    _null_cache[tgt] = nulls.alexander_bloch(lgs_data, atlas=tgt[0], density=tgt[1],
                                             n_perm=N_PERM, seed=SEED)
    return _null_cache[tgt]


def annot_in_common(src, desc, space, den):
    annot = fetch_annotation(source=src, desc=desc, space=space, den=den, verbose=0)
    if space == "MNI152":
        return load_data(transforms.mni152_to_fslr(annot, "32k"))
    return load_data(annot)


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def main():
    print(f"[battery] wb_command available: {HAVE_WB}", flush=True)
    # group by (source,desc); choose best workbench-free space
    options = {}
    for src, desc, space, den in available_annotations():
        options.setdefault((src, desc), []).append((space, den))
    todo = []
    for (src, desc), spaces in options.items():
        ranked = sorted(spaces, key=lambda sd: space_rank(*sd))
        best = ranked[0]
        if space_rank(*best) >= 98:
            continue  # unreachable without workbench (skip)
        todo.append((src, desc, best[0], best[1]))
    todo.sort()
    print(f"[battery] {len(todo)} annotations reachable "
          f"({sum(1 for t in todo if (t[0],t[1]) in APRIORI)} a-priori)", flush=True)

    rows, log = [], []
    for i, (src, desc, space, den) in enumerate(todo):
        tag = APRIORI.get((src, desc))
        try:
            lgs_data, tgt = lgs_in(space, den)
            rot = spin_null(tgt)
            g_data = annot_in_common(src, desc, space, den)
            r, p = stats.compare_images(lgs_data, g_data, nulls=rot, metric="spearmanr",
                                        ignore_zero=True, nan_policy="omit")
            rows.append([f"{src}:{desc}", f"{tgt[0]}/{tgt[1]}", r, p,
                         "apriori" if tag else "exploratory", tag or ""])
            print(f"  [{i+1}/{len(todo)}] {src}:{desc:16s} r={r:+.3f} p={p:.4f}"
                  f"{'  <'+tag if tag else ''}", flush=True)
        except Exception as e:
            log.append(f"{src}:{desc} ({space}/{den}) FAILED: {str(e)[:150]}")
            print(f"  [{i+1}/{len(todo)}] {src}:{desc} FAILED: {str(e)[:70]}", flush=True)

    def add_fdr(subset):
        ps = [row[3] for row in subset]
        q = bh_fdr(ps) if ps else []
        for row, qq in zip(subset, q):
            row.append(qq)
        for row in subset[len(q):]:
            row.append(float("nan"))

    expl = [r for r in rows if r[4] == "exploratory"]
    apri = [r for r in rows if r[4] == "apriori"]
    add_fdr(expl); add_fdr(apri)

    with (OUT / "lgs_battery.tsv").open("w") as fh:
        fh.write("map\tspace\tspearman_r\tp_spin\tfamily\tapriori_label\tp_fdr\n")
        for m, sp, r, p, fam, lab, q in sorted(rows, key=lambda x: x[3]):
            fh.write(f"{m}\t{sp}\t{r:.4f}\t{p:.4f}\t{fam}\t{lab}\t{q:.4f}\n")
    with (OUT / "battery_failures.log").open("w") as fh:
        fh.write("\n".join(log) + ("\n" if log else ""))
    print(f"\n[battery] wrote {len(rows)} results, {len(log)} failures", flush=True)
    print("BATTERY_DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
