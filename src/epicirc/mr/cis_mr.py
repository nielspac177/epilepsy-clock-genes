"""cis-MR: does genetically-predicted expression of each clock gene affect epilepsy?

Given behavioural sleep-trait MR was null, this tests the molecular clock directly. For each core
clock gene we take its strongest blood cis-eQTL (eQTLGen) as a single instrument, match it to the
ILAE outcome by CHR:POS, harmonize to the assessed allele, and compute the Wald ratio.

APPROXIMATIONS (stated in any write-up): exposure effect on standardized expression is approximated
as beta_exp = Z / sqrt(N), se_exp = 1/sqrt(N) (eQTLGen ships Z, not beta+MAF); one instrument per
gene avoids LD but ignores allelic heterogeneity; blood eQTLs may not reflect brain (PER1/RORB have
no significant blood cis-eQTL here). A proper analysis uses brain eQTLs (GTEx/MetaBrain) +
colocalization to exclude LD-confounding — this is a screen, not proof.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from epicirc.data.harmonize import harmonize

# clock_cis_eqtls.tsv columns (0-indexed)
EQ = {"P": 0, "SNP": 1, "CHR": 2, "POS": 3, "EA": 4, "OA": 5, "Z": 6, "SYMBOL": 8, "N": 12}
IL = {"CHR": 0, "BP": 1, "A1": 3, "A2": 4, "BETA": 11, "SE": 12}


def top_cis_per_gene(eqtl_tsv: str) -> dict[str, dict]:
    """Strongest (smallest-P) cis-eQTL per gene."""
    best: dict[str, dict] = {}
    with open(eqtl_tsv) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= EQ["N"]:
                continue
            try:
                p = float(f[EQ["P"]]); z = float(f[EQ["Z"]]); n = float(f[EQ["N"]])
            except ValueError:
                continue
            sym = f[EQ["SYMBOL"]]
            if sym not in best or p < best[sym]["p"]:
                best[sym] = {"snp": f[EQ["SNP"]], "chr": f[EQ["CHR"]], "pos": int(f[EQ["POS"]]),
                             "ea": f[EQ["EA"]].upper(), "oa": f[EQ["OA"]].upper(),
                             "z": z, "n": n, "p": p}
    return best


def match_outcome(ilae_tbl: str, keys: set[tuple[str, int]]) -> dict[tuple[str, int], dict]:
    found = {}
    with open(ilae_tbl) as fh:
        next(fh)
        for ln in fh:
            f = ln.split()
            if len(f) < 13:
                continue
            try:
                key = (f[IL["CHR"]], int(f[IL["BP"]]))
            except ValueError:
                continue
            if key in keys:
                try:
                    found[key] = {"a1": f[IL["A1"]].upper(), "a2": f[IL["A2"]].upper(),
                                  "beta": float(f[IL["BETA"]]), "se": float(f[IL["SE"]])}
                except ValueError:
                    continue
    return found


def _p_from_z(z: float) -> float:
    return 2.0 * (0.5 * math.erfc(abs(z) / math.sqrt(2.0)))


def run(eqtl_tsv: str, ilae_tbl: str):
    cis = top_cis_per_gene(eqtl_tsv)
    keys = {(v["chr"], v["pos"]) for v in cis.values()}
    out = match_outcome(ilae_tbl, keys)
    rows = []
    for gene, e in cis.items():
        o = out.get((e["chr"], e["pos"]))
        if o is None:
            continue
        h = harmonize(o["beta"], o["a1"], o["a2"], e["ea"], e["oa"])
        if not h.usable:
            continue
        beta_exp = e["z"] / math.sqrt(e["n"])       # per-allele effect on standardized expression
        se_exp = 1.0 / math.sqrt(e["n"])
        wald = h.beta / beta_exp
        se_wald = o["se"] / abs(beta_exp)            # first-order delta method
        z = wald / se_wald
        rows.append({"gene": gene, "snp": e["snp"], "eqtl_p": e["p"], "beta_exp": beta_exp,
                     "se_exp": se_exp, "beta_out": h.beta, "se_out": o["se"],
                     "wald": wald, "se_wald": se_wald, "p": _p_from_z(z)})
    return sorted(rows, key=lambda r: r["p"])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eqtl", default="data/raw/eqtl/clock_cis_eqtls.tsv")
    ap.add_argument("--outcome", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    rows = run(args.eqtl, args.outcome)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["gene", "snp", "eqtl_p", "beta_exp", "beta_out", "se_out", "wald", "se_wald", "p"]
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                               for c in cols) + "\n")
    nsig = sum(1 for r in rows if r["p"] < 0.05)
    print(f"[cis-mr] {len(rows)} clock genes tested; {nsig} nominal p<0.05 -> {out}",
          file=sys.stderr)
    for r in rows[:6]:
        print(f"  {r['gene']:8} wald={r['wald']:+.3f} p={r['p']:.2e} (eQTL p={r['eqtl_p']:.1e})",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
