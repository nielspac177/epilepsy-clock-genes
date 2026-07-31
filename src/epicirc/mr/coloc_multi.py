"""Two-GWAS and three-trait colocalization at a locus (extends mr/coloc.py to arbitrary traits).

The existing `mr/coloc.py` is hardcoded to a clock-gene brain-eQTL (p-value ABF) vs the ILAE
epilepsy GWAS. This module generalizes the same Wakefield-ABF machinery so BOTH sides can be
full GWAS (beta/se), and adds a coloc-consistent three-trait test (moloc-style enumeration of the
15 sharing configurations for traits a,b,c).

Key simplification: for a GWAS with (beta, se) the Wakefield log-ABF is
    labf = 0.5*(log(1-r) + r*z^2),  r = W/(se^2 + W),  z^2 = (beta/se)^2
which needs NO MAF (unlike the sdY quantitative-eQTL approximation in coloc.py). So a two-GWAS
coloc only needs matched (chr:pos) beta/se on each side. Priors follow coloc/moloc defaults:
one causal variant in one trait -> P1; a variant shared by two traits -> P12; shared by three -> P123.
Independent variants multiply P1. This reduces EXACTLY to coloc.abf for two traits (H0..H4).

Sign-agnostic (uses z^2), so allele orientation between traits does not matter for the posterior.
"""
from __future__ import annotations

import argparse
import math
import sys
from itertools import product
from pathlib import Path

from epicirc.mr.coloc import _labf, _logsumexp  # reuse the audited ABF building blocks

# prior variances (same as coloc.py): case/control log-odds vs quantitative
W_CC = 0.2 ** 2      # binary trait (epilepsy) prior variance on log-OR
W_QUANT = 0.15 ** 2  # quantitative trait (chronotype, expression) prior variance
P1 = 1e-4            # a variant is causal for one trait
P12 = 1e-5           # a single variant is shared by two traits (coloc default)
P123 = 1e-6          # a single variant is shared by three traits


def load_gwas_region(path: str, cols: dict, chrom: str, start: int, end: int,
                     sep: str | None = None, gz: bool = False) -> dict[tuple[str, int], tuple[float, float]]:
    """(chr,pos) -> (beta, se) for rows in [start,end] on `chrom`.

    cols maps logical names CHR,BP,BETA,SE (0-based indices) to file columns. `chrom` is compared
    with any leading 'chr' stripped. OR columns can be given via cols['OR']=idx (converted to logOR).
    """
    import gzip
    opener = gzip.open if (gz or path.endswith(".gz")) else open
    ci, bi = cols["CHR"], cols["BP"]
    has_or = "OR" in cols
    beta_i = cols.get("BETA")
    si = cols["SE"]
    out: dict[tuple[str, int], tuple[float, float]] = {}
    with opener(path, "rt") as fh:
        next(fh)  # header
        for ln in fh:
            f = ln.split(sep) if sep else ln.split()
            if len(f) <= max(ci, bi, si, beta_i if beta_i is not None else 0,
                             cols["OR"] if has_or else 0):
                continue
            c = f[ci].replace("chr", "")
            if c != chrom:
                continue
            try:
                pos = int(f[bi])
                if pos < start or pos > end:
                    continue
                se = float(f[si])
                beta = math.log(float(f[cols["OR"]])) if has_or else float(f[beta_i])
            except (ValueError, KeyError):
                continue
            if se > 0:
                out[(c, pos)] = (beta, se)
    return out


def _labf_gwas(beta: float, se: float, w: float) -> float:
    return _labf((beta / se) ** 2, se ** 2, w)


def coloc_pairwise(a: dict, b: dict, w_a: float = W_CC, w_b: float = W_QUANT) -> dict:
    """coloc.abf for two GWAS (beta/se each). Returns PP.H0..H4 over matched SNPs."""
    shared = sorted(set(a) & set(b))
    la, lb = [], []
    for k in shared:
        ba, sa = a[k]; bb, sb = b[k]
        if sa <= 0 or sb <= 0:
            continue
        la.append(_labf_gwas(ba, sa, w_a))
        lb.append(_labf_gwas(bb, sb, w_b))
    n = len(la)
    if n < 5:
        return {"n_snps": n, "PP.H4": float("nan")}
    l1 = _logsumexp(la)
    l2 = _logsumexp(lb)
    l4 = _logsumexp([x + y for x, y in zip(la, lb)])
    l3 = math.log(max(math.exp(l1 + l2) - math.exp(l4), 1e-300))
    terms = [0.0, math.log(P1) + l1, math.log(P1) + l2,
             2 * math.log(P1) + l3, math.log(P12) + l4]
    denom = _logsumexp(terms)
    pp = [math.exp(t - denom) for t in terms]
    return {"n_snps": n, "PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2],
            "PP.H3": pp[3], "PP.H4": pp[4]}


# ---- three-trait moloc (coloc-consistent enumeration of the 15 configurations) ----
# A configuration assigns each of traits (a,b,c) a label in {0,1,2,3}: 0 = no causal variant,
# equal nonzero label = shared causal variant. logBF(config) = sum over groups (traits sharing a
# label) of logsumexp_i( sum_{t in group} labf[t][i] ); null traits contribute 0. Config prior =
# product over groups: singleton -> P1, pair -> P12, triple -> P123.
def _canonical_configs():
    """Return the 15 canonical (a,b,c) label tuples (dedup relabelings of group ids)."""
    seen, out = set(), []
    for labels in product(range(4), repeat=3):
        # canonicalize: relabel nonzero groups in order of first appearance to 1,2,3
        mapping, nxt, canon = {}, 1, []
        for x in labels:
            if x == 0:
                canon.append(0)
            else:
                if x not in mapping:
                    mapping[x] = nxt; nxt += 1
                canon.append(mapping[x])
        t = tuple(canon)
        if t not in seen:
            seen.add(t); out.append(t)
    return out


def _config_name(cfg):
    names = "abc"
    groups: dict[int, str] = {}
    for t, lab in enumerate(cfg):
        if lab:
            groups.setdefault(lab, "")
            groups[lab] += names[t]
    if not groups:
        return "H0"
    return "_".join("".join(sorted(g)) for g in (groups[k] for k in sorted(groups)))


def moloc3(a: dict, b: dict, c: dict, w=(W_CC, W_QUANT, W_QUANT)) -> dict:
    """Three-trait coloc over traits a,b,c (each a {(chr,pos):(beta,se)} dict).

    Returns per-configuration posteriors plus aggregates: PP_abc (all three share one variant),
    PP_any_pair (some pair shares), PP_none_shared.
    """
    keys = sorted(set(a) & set(b) & set(c))
    labf = [[], [], []]
    dicts = (a, b, c)
    for k in keys:
        ok = True
        vals = []
        for t in range(3):
            beta, se = dicts[t][k]
            if se <= 0:
                ok = False; break
            vals.append(_labf_gwas(beta, se, w[t]))
        if ok:
            for t in range(3):
                labf[t].append(vals[t])
    n = len(labf[0])
    if n < 5:
        return {"n_snps": n, "PP_abc": float("nan")}

    # per-trait marginal logsumexp (used for singleton groups)
    single = [_logsumexp(labf[t]) for t in range(3)]

    def group_logbf(group: tuple[int, ...]) -> float:
        if len(group) == 1:
            return single[group[0]]
        # shared variant: sum labf across the group at the same SNP, then logsumexp over SNPs
        return _logsumexp([sum(labf[t][i] for t in group) for i in range(n)])

    def group_prior(group):
        return {1: P1, 2: P12, 3: P123}[len(group)]

    terms, cfgs = [], _canonical_configs()
    for cfg in cfgs:
        groups: dict[int, list[int]] = {}
        for t, lab in enumerate(cfg):
            if lab:
                groups.setdefault(lab, []).append(t)
        logbf = sum(group_logbf(tuple(g)) for g in groups.values())
        logprior = sum(math.log(group_prior(tuple(g))) for g in groups.values())
        terms.append((cfg, logprior + logbf))

    denom = _logsumexp([t[1] for t in terms])
    post = {_config_name(cfg): math.exp(val - denom) for cfg, val in terms}
    # aggregates
    pp_abc = post.get("abc", 0.0)
    pp_any_pair = sum(v for name, v in post.items()
                      if any(len(g) >= 2 for g in name.split("_")) and name != "H0")
    return {"n_snps": n, "PP_abc": pp_abc, "PP_any_pair": pp_any_pair,
            "PP_none_shared": 1.0 - pp_any_pair, **{f"cfg.{k}": v for k, v in post.items()}}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["pair", "moloc3"], required=True)
    ap.add_argument("--out", required=True)
    # generic per-trait spec: path:CHR,BP,BETA,SE[,OR] indices ; use OR= for odds-ratio columns
    ap.add_argument("--trait", action="append", required=True,
                    help="name=path=chr,bp,beta,se[,or][,gz]  (repeat 2x for pair, 3x for moloc3)")
    ap.add_argument("--region", required=True, help="chr:start-end (hg19)")
    ap.add_argument("--wtype", action="append", default=None,
                    help="cc|quant per trait (order matches --trait); default cc,quant,quant")
    args = ap.parse_args(argv)

    chrom, rng = args.region.split(":")
    start, end = (int(x) for x in rng.split("-"))

    def parse_trait(spec):
        name, path, colspec = spec.split("=", 2)
        parts = colspec.split(",")
        cols = {"CHR": int(parts[0]), "BP": int(parts[1]), "SE": int(parts[3])}
        if "or" in parts:
            cols["OR"] = int(parts[2])
        else:
            cols["BETA"] = int(parts[2])
        gz = "gz" in parts
        return name, load_gwas_region(path, cols, chrom.replace("chr", ""), start, end, gz=gz)

    traits = [parse_trait(s) for s in args.trait]
    wmap = {"cc": W_CC, "quant": W_QUANT}
    wt = args.wtype or (["cc", "quant", "quant"][:len(traits)])
    ws = [wmap[x] for x in wt]

    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    if args.mode == "pair":
        (na, a), (nb, b) = traits
        r = coloc_pairwise(a, b, ws[0], ws[1])
        r["trait_a"], r["trait_b"] = na, nb
        print(f"[coloc pair] {na} x {nb}: PP.H4={r.get('PP.H4')} (n={r['n_snps']})", file=sys.stderr)
    else:
        (na, a), (nb, b), (nc, c) = traits
        r = moloc3(a, b, c, tuple(ws))
        r["traits"] = f"{na},{nb},{nc}"
        print(f"[moloc3] {na},{nb},{nc}: PP_abc={r.get('PP_abc')} "
              f"PP_any_pair={r.get('PP_any_pair')} (n={r['n_snps']})", file=sys.stderr)
    with out.open("w") as fh:
        for k, v in r.items():
            fh.write(f"{k}\t{v:.6f}\n" if isinstance(v, float) else f"{k}\t{v}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
