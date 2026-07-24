"""Phase D/E foundation: build the Allen Human Brain Atlas region x gene expression matrix.

Uses abagen with the Desikan-Killiany atlas (MNI152 volumetric; 83 cortical+subcortical regions,
so the thalamocortical LGS network is covered) and abagen/Arnatkeviciute-2019 default processing.
Left hemisphere is retained for all regions (donor coverage), and we also emit region centroids in
MNI space for the distance-based spatial-autocorrelation null (brainsmash).

Outputs (results/spatial/):
  ahba_dk_expression.csv   region x gene (normalized microarray expression)
  ahba_dk_info.csv         atlas region labels (id, label, hemisphere, structure)
  ahba_dk_centroids.csv    region id -> MNI x,y,z centroid (mm), for the variogram null
"""
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import abagen

OUT = Path("results/spatial")
OUT.mkdir(parents=True, exist_ok=True)


def main():
    atlas = abagen.fetch_desikan_killiany(native=False)   # MNI152 volumetric DK
    print("atlas image:", atlas["image"])

    # NOTE: 2 of 6 AHBA donors have a stale Allen download URL (file 178238266 -> 404) in abagen
    # 0.1.3, so we use the 4 available donors, which INCLUDE both bilateral donors (9861, 10021).
    # A 4-donor AHBA is a documented limitation (reduced subcortical coverage in particular).
    DONORS = ["9861", "10021", "12876", "14380"]
    expr = abagen.get_expression_data(
        atlas["image"], atlas["info"],
        donors=DONORS,
        lr_mirror="bidirectional",   # mirror samples across hemispheres to improve coverage
        missing="centroids",
        norm_matched=False,
        verbose=1,
    )
    expr.to_csv(OUT / "ahba_dk_expression.csv")
    print(f"expression: {expr.shape[0]} regions x {expr.shape[1]} genes")

    info = pd.read_csv(atlas["info"])
    info.to_csv(OUT / "ahba_dk_info.csv", index=False)

    # region centroids in MNI mm (label image -> voxel COM -> world coords)
    img = nib.load(str(atlas["image"]))
    data = np.asarray(img.dataobj)
    aff = img.affine
    rows = []
    for lab in np.unique(data):
        if lab == 0:
            continue
        ijk = np.argwhere(data == lab).mean(0)
        xyz = aff[:3, :3] @ ijk + aff[:3, 3]
        rows.append((int(lab), *xyz))
    pd.DataFrame(rows, columns=["id", "x", "y", "z"]).to_csv(OUT / "ahba_dk_centroids.csv", index=False)
    print(f"centroids: {len(rows)} regions -> {OUT/'ahba_dk_centroids.csv'}")

    # sanity: are the clock genes present?
    CLOCK = ["ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2",
             "RORA", "RORB", "RORC", "CSNK1D", "CSNK1E", "FBXL3", "DBP", "NFIL3", "TIMELESS",
             "BHLHE40", "BHLHE41", "CIART", "ARNTL2"]
    present = [g for g in CLOCK if g in expr.columns]
    print(f"clock genes present in AHBA: {len(present)}/23 -> {present}")


if __name__ == "__main__":
    main()
