"""Harmonize GGE z-scores to the 1000G LD panel allele order for SuSiE-RSS fine-mapping.

plink2 --r-unphased counts the .bim A1 (col5) allele, so the sign of each GGE z must be flipped
whenever the GGE effect allele equals the panel's A2 rather than A1. SNPs whose alleles don't match
the panel (or are strand-ambiguous A/T, C/G) are dropped. Emits an aligned z vector, a subsetted LD
matrix, and a SNP table (id, bp) for annotation.
"""
import numpy as np

AMB = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def main():
    # panel alleles in matrix order
    order, a1, a2 = [], {}, {}
    with open("results/finemap/locus_1000g.bim") as fh:
        for ln in fh:
            f = ln.split()
            order.append(f[1]); a1[f[1]] = f[4].upper(); a2[f[1]] = f[5].upper()
    # GGE stats
    gge = {}
    with open("results/finemap/gge_per1_locus.tsv") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            gge[f[0]] = (f[3].upper(), f[4].upper(), float(f[6]), int(f[2]), float(f[8]))  # A1,A2,Z,BP,NEFF

    keep_idx, z, snps, bps, neffs, dropped = [], [], [], [], [], 0
    for i, snp in enumerate(order):
        if snp not in gge:
            dropped += 1; continue
        gA1, gA2, gz, bp, neff = gge[snp]
        pa1, pa2 = a1[snp], a2[snp]
        if (gA1, gA2) in AMB:            # strand-ambiguous -> drop
            dropped += 1; continue
        if {gA1, gA2} != {pa1, pa2}:     # allele mismatch -> drop
            dropped += 1; continue
        sign = 1.0 if gA1 == pa1 else -1.0
        keep_idx.append(i); z.append(sign * gz); snps.append(snp); bps.append(bp); neffs.append(neff)

    keep_idx = np.array(keep_idx)
    print(f"kept {len(keep_idx)} / {len(order)} SNPs (dropped {dropped}: missing/mismatch/ambiguous)")

    # subset the LD matrix
    R = np.loadtxt("results/finemap/locus_ld.unphased.vcor1")
    R = R[np.ix_(keep_idx, keep_idx)]
    np.fill_diagonal(R, 1.0)
    np.savetxt("results/finemap/ld_aligned.txt", R, fmt="%.5f")
    np.savetxt("results/finemap/z_aligned.txt", np.array(z), fmt="%.6f")
    with open("results/finemap/snps_aligned.tsv", "w") as fh:
        fh.write("snp\tbp\tz\tneff\n")
        for s, b, zz, n in zip(snps, bps, z, neffs):
            fh.write(f"{s}\t{b}\t{zz:.4f}\t{n:.1f}\n")
    print(f"median Neff = {np.median(neffs):.0f}; lead |z| = {np.max(np.abs(z)):.2f}")


if __name__ == "__main__":
    main()
