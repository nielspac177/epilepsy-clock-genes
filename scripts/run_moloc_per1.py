"""Three-trait colocalization at PER1: GGE (epilepsy) + chronotype + PER1 brain cis-eQTL.

All three inputs are on disk (no new download): GGE = ILAE region slice, chronotype = Jones-2019
region slice, eQTL = PsychENCODE PER1 cis-eQTL (slope + nominal p -> beta/se). Answers: does a
single variant at 17p13.1 drive epilepsy risk, morningness, AND PER1 expression (PP_abc), or are
the signals independent (PP_none_shared)? Also reports all three pairwise colocs.
"""
from __future__ import annotations
import math
from pathlib import Path

from epicirc.mr.coloc import _z2_from_p
from epicirc.mr.coloc_multi import (load_gwas_region, coloc_pairwise, moloc3,
                                    W_CC, W_QUANT)

OUT = Path("results/coloc_multi"); OUT.mkdir(parents=True, exist_ok=True)
PE_FILE = "data/raw/eqtl/psychencode_clock_cis.tsv"
PER1_ENSG = "ENSG00000179094"
# PsychENCODE columns (0-based): GENE0 ... varCHR8 varPOS9 ... nominal_p11 slope12
PE = {"GENE": 0, "CHR": 8, "POS": 9, "P": 11, "SLOPE": 12}


def load_eqtl_betase(path: str, ensg: str) -> dict:
    """(chr,pos) -> (beta=slope, se) for the gene's cis-eQTLs; se from |slope|/z, keep min-p/pos."""
    best_p: dict = {}
    out: dict = {}
    with open(path) as fh:
        for ln in fh:
            f = ln.split()
            if len(f) <= PE["SLOPE"] or not f[PE["GENE"]].startswith(ensg):
                continue
            try:
                key = (f[PE["CHR"]].replace("chr", ""), int(f[PE["POS"]]))
                p = float(f[PE["P"]]); slope = float(f[PE["SLOPE"]])
            except ValueError:
                continue
            if key in best_p and p >= best_p[key]:
                continue
            z2 = _z2_from_p(p)
            if z2 <= 0:
                continue
            se = abs(slope) / math.sqrt(z2)
            if se > 0:
                best_p[key] = p
                out[key] = (slope, se)
    return out


def main() -> int:
    gge = load_gwas_region("results/coloc_multi/gge_per1.tbl",
                           {"CHR": 0, "BP": 1, "BETA": 11, "SE": 12}, "17", 7_850_000, 8_250_000)
    chrono = load_gwas_region("results/coloc_multi/chronotype_per1.txt",
                              {"CHR": 1, "BP": 2, "BETA": 7, "SE": 8}, "17", 7_850_000, 8_250_000)
    eqtl = load_eqtl_betase(PE_FILE, PER1_ENSG)
    print(f"n SNPs: GGE={len(gge)} chronotype={len(chrono)} PER1-eQTL={len(eqtl)}")
    print(f"3-way overlap: {len(set(gge) & set(chrono) & set(eqtl))}")

    rows = []
    # pairwise
    for na, a, wa, nb, b, wb in [
        ("GGE", gge, W_CC, "chronotype", chrono, W_QUANT),
        ("GGE", gge, W_CC, "PER1_eQTL", eqtl, W_QUANT),
        ("chronotype", chrono, W_QUANT, "PER1_eQTL", eqtl, W_QUANT),
    ]:
        r = coloc_pairwise(a, b, wa, wb)
        rows.append((f"pair:{na}~{nb}", r))
        print(f"  {na}~{nb}: n={r['n_snps']} PP.H4={r.get('PP.H4')}")

    m = moloc3(gge, chrono, eqtl, (W_CC, W_QUANT, W_QUANT))
    rows.append(("moloc3:GGE,chronotype,PER1_eQTL", m))
    print(f"  moloc3: n={m['n_snps']} PP_abc={m.get('PP_abc')} "
          f"PP_any_pair={m.get('PP_any_pair')} PP_none_shared={m.get('PP_none_shared')}")

    with (OUT / "moloc_per1.tsv").open("w") as fh:
        fh.write("analysis\tkey\tvalue\n")
        for name, r in rows:
            for k, v in r.items():
                fh.write(f"{name}\t{k}\t{v:.6f}\n" if isinstance(v, float)
                         else f"{name}\t{k}\t{v}\n")
    print(f"wrote {OUT/'moloc_per1.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
