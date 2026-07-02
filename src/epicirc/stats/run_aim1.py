"""Aim 1 synthesis (Snakemake `aim1_heterogeneity` rule).

Combines the per-phenotype circadian gene-set enrichment (MAGMA arm) into:
  1. an enrichment table, power-gated by the pre-registered estimable_cells.tsv and BH-corrected
     across estimable cells (exploratory) with Bonferroni on the confirmatory family;
  2. the primary focal-vs-GGE heterogeneity test, using the shared-control correlation rho
     (estimated as the focal/GGE effect correlation) so the difference SE is not double-counted.
Non-significant differences are reported with the minimum detectable difference, never as
"no difference" (docs/adr/0003).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from epicirc.ldsc.run_rg import rg_from_files
from epicirc.magma.geneset import parse_gsa_out
from epicirc.stats.heterogeneity import heterogeneity_test
from epicirc.stats.multipletesting import benjamini_hochberg, bonferroni

CONFIRMATORY = {"focal_epilepsy", "generalized_gge", "all_epilepsy"}


def _load_estimable(path: str) -> dict[str, bool]:
    est = {}
    p = Path(path)
    if not p.exists():
        return est
    lines = p.read_text().splitlines()
    for ln in lines[1:]:
        f = ln.split("\t")
        est[f[0]] = f[-1].strip() == "yes"
    return est


def _core_beta(gsa_path: Path):
    for r in parse_gsa_out(gsa_path):
        if r.set_name == "circadian_core":
            return r
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/analysis_params.yaml")
    ap.add_argument("--magma-dir", default="results/magma")
    ap.add_argument("--processed-dir", default="data/processed")
    ap.add_argument("--estimable", default="results/estimable_cells.tsv")
    ap.add_argument("--out", default="results/aim1_heterogeneity.tsv")
    args = ap.parse_args(argv)

    estimable = _load_estimable(args.estimable)
    magma_dir = Path(args.magma_dir)

    # ---- enrichment per phenotype (circadian_core) ----
    rows = []
    for gsa in sorted(magma_dir.glob("circadian_core__*.gsa.out")):
        pheno = gsa.stem.replace("circadian_core__", "").replace(".gsa", "")
        res = _core_beta(gsa)
        if res is None:
            continue
        rows.append({"phenotype": pheno, "beta": res.beta, "se": res.beta_se, "p": res.p,
                     "estimable": estimable.get(pheno, True),
                     "family": "confirmatory" if pheno in CONFIRMATORY else "exploratory"})

    conf = [r for r in rows if r["family"] == "confirmatory"]
    expl = [r for r in rows if r["family"] == "exploratory" and r["estimable"]]
    for r, padj in zip(conf, bonferroni([r["p"] for r in conf])):
        r["p_adj"] = padj
    for r, padj in zip(expl, benjamini_hochberg([r["p"] for r in expl])):
        r["p_adj"] = padj
    for r in rows:
        r.setdefault("p_adj", float("nan"))  # not-estimable exploratory cells: no test

    enrich_out = Path(args.out).with_name("aim1_enrichment.tsv")
    with enrich_out.open("w") as fh:
        fh.write("phenotype\tfamily\testimable\tbeta\tse\tp\tp_adj\tstatus\n")
        for r in sorted(rows, key=lambda x: (x["family"], x["phenotype"])):
            status = ("not_estimable" if not r["estimable"]
                      else ("sig" if r["p_adj"] < 0.05 else "ns"))
            padj = "NA" if r["p_adj"] != r["p_adj"] else f"{r['p_adj']:.3e}"
            fh.write(f"{r['phenotype']}\t{r['family']}\t{int(r['estimable'])}\t{r['beta']:.4f}\t"
                     f"{r['se']:.4f}\t{r['p']:.3e}\t{padj}\t{status}\n")

    # ---- primary focal-vs-GGE heterogeneity ----
    focal = _core_beta(magma_dir / "circadian_core__focal_epilepsy.gsa.out")
    gge = _core_beta(magma_dir / "circadian_core__generalized_gge.gsa.out")
    het_out = Path(args.out)
    with het_out.open("w") as fh:
        fh.write("contrast\tbeta_focal\tse_focal\tbeta_gge\tse_gge\trho\tdelta\tse_delta\tz\tp\t"
                 "min_detectable_delta\tverdict\n")
        if focal and gge:
            fp = Path(args.processed_dir) / "focal_epilepsy.sumstats.gz"
            gp = Path(args.processed_dir) / "generalized_gge.sumstats.gz"
            rho = 0.0
            if fp.exists() and gp.exists():
                rho = max(-0.99, min(0.99, rg_from_files(str(fp), str(gp))["rg"]))
            res = heterogeneity_test(focal.beta, focal.beta_se, gge.beta, gge.beta_se, rho=rho)
            mdd = 1.96 * res.se_delta  # ~min detectable difference at 80% power (approx 2.8*se/1.4)
            verdict = ("detected" if res.p_two_sided < 0.05
                       else ("equivalent" if abs(res.delta) + 1.96 * res.se_delta < 0.2
                             else "inconclusive"))
            fh.write(f"focal_vs_gge\t{focal.beta:.4f}\t{focal.beta_se:.4f}\t{gge.beta:.4f}\t"
                     f"{gge.beta_se:.4f}\t{rho:.4f}\t{res.delta:.4f}\t{res.se_delta:.4f}\t"
                     f"{res.z:.3f}\t{res.p_two_sided:.3e}\t{mdd:.4f}\t{verdict}\n")
            print(f"[aim1] focal-vs-GGE delta={res.delta:.3f} p={res.p_two_sided:.2e} "
                  f"verdict={verdict} (rho={rho:.2f})", file=sys.stderr)
        else:
            fh.write("focal_vs_gge\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_inputs\n")
    print(f"[aim1] enrichment -> {enrich_out}; heterogeneity -> {het_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
