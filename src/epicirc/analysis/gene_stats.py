"""Cache the LD-robust gene-level statistic (top-SNP chi2) for a whole GWAS once.

Scanning a 430 MB ILAE .tbl is ~90 s; every Tier-1 analysis (leave-one-out, difference test,
control sets, window sensitivity) needs the *same* per-gene top-chi2. Compute it once here and
write results/gene_stats/<label>.tsv, then have geneset_enrich consume the cache.

Output columns: gene chr start end length n_snp top_chi2 top_bp
The cache is (phenotype, flank)-specific: top_chi2/n_snp already reflect the ±flank window. The
window-sensitivity grid regenerates one cache per flank (still cheap vs. re-scanning per analysis).
"""
from __future__ import annotations

import argparse
import gzip
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

IL = {"CHR": 0, "BP": 1, "Z": 8}  # ILAE .tbl 0-based: CHR, BP, ..., Z-score at col 8
# FinnGen R12: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ... -> Z=beta/sebeta
FG = {"CHR": 0, "BP": 1, "BETA": 8, "SE": 9, "AF": 10}
FG_MAF_MIN = 0.01   # drop rare variants: tiny-SE rare imputed markers inflate beta/sebeta absurdly


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_genes(path: str):
    by_chr = defaultdict(list)  # chr -> [(start, end, gene)]
    meta = {}                   # gene -> (chr, start, end, length)
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            g, c, s, e = f[0], f[1], f[2], f[3]
            try:
                s, e = int(s), int(e)
            except ValueError:
                continue  # header
            by_chr[c].append((s, e, g))
            meta[g] = (c, s, e, e - s)
    for c in by_chr:
        by_chr[c].sort()
    return by_chr, meta


def index_snps(tbl: str, keep_chr, source="ilae"):
    pos_by_chr = defaultdict(list)
    chi_by_chr = defaultdict(list)
    with _open(tbl) as fh:
        next(fh)
        for ln in fh:
            # FinnGen is TAB-delimited with empty fields (rsID-less indels leave an empty rsids
            # column); whitespace split() would collapse the empty field and shift every later
            # column left, corrupting beta/sebeta/af. Must split on tab. ILAE .tbl has no empty
            # fields, so its whitespace split is safe and preserved (bit-identical reproduction).
            f = ln.rstrip("\n").split("\t") if source == "finngen" else ln.split()
            if len(f) < 10:
                continue
            c = f[IL["CHR"]].replace("chr", "")
            if c not in keep_chr:
                continue
            try:
                bp = int(f[FG["BP"] if source == "finngen" else IL["BP"]])
                if source == "finngen":
                    af = float(f[FG["AF"]])
                    if af < FG_MAF_MIN or af > 1 - FG_MAF_MIN:
                        continue
                    se = float(f[FG["SE"]])
                    if se == 0:
                        continue
                    z = float(f[FG["BETA"]]) / se
                else:
                    z = float(f[IL["Z"]])
            except (ValueError, ZeroDivisionError):
                continue
            pos_by_chr[c].append(bp)
            chi_by_chr[c].append(z * z)
    for c in pos_by_chr:
        order = sorted(range(len(pos_by_chr[c])), key=lambda i: pos_by_chr[c][i])
        pos_by_chr[c] = [pos_by_chr[c][i] for i in order]
        chi_by_chr[c] = [chi_by_chr[c][i] for i in order]
    return pos_by_chr, chi_by_chr


def compute(tbl: str, genes_path: str, flank_kb: int, source="ilae"):
    flank = flank_kb * 1000
    by_chr, meta = load_genes(genes_path)
    pos_by_chr, chi_by_chr = index_snps(tbl, set(by_chr), source)
    rows = []
    for c, windows in by_chr.items():
        pos, chi = pos_by_chr.get(c, []), chi_by_chr.get(c, [])
        if not pos:
            continue
        for s, e, g in windows:
            lo = bisect_left(pos, s - flank)
            hi = bisect_right(pos, e + flank)
            if hi <= lo:
                continue
            seg = chi[lo:hi]
            mx = max(seg)
            top_bp = pos[lo + seg.index(mx)]
            _, gs, ge, ln = meta[g]
            rows.append((g, c, gs, ge, ln, hi - lo, mx, top_bp))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", default="reference/genes_hg19.tsv")
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank-kb", type=int, default=50)
    ap.add_argument("--source", choices=["ilae", "finngen"], default="ilae")
    args = ap.parse_args(argv)

    print(f"[gene_stats] scanning {args.outcome} (flank {args.flank_kb}kb, {args.source}) ...", file=sys.stderr)
    rows = compute(args.outcome, args.genes, args.flank_kb, args.source)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("gene\tchr\tstart\tend\tlength\tn_snp\ttop_chi2\ttop_bp\n")
        for g, c, s, e, ln, n, chi, bp in rows:
            fh.write(f"{g}\t{c}\t{s}\t{e}\t{ln}\t{n}\t{chi:.6f}\t{bp}\n")
    print(f"[gene_stats] wrote {len(rows)} genes -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
