"""Aggregate results into report.html + dashboard.json (Snakemake `report` rule).

dashboard.json feeds dashboard/index.html; report.html is a standalone summary. Kept dependency
-free (stdlib) so it runs in the synthetic path.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def build_dashboard(results: Path) -> dict:
    het = _read_tsv(results / "aim1_heterogeneity.tsv")
    enrich = _read_tsv(results / "aim1_enrichment.tsv")
    estimable = _read_tsv(results / "estimable_cells.tsv")
    aim3 = _read_tsv(results / "aim3_severity.tsv")

    primary = het[0] if het else {}
    n_est = sum(1 for r in estimable if r.get("estimable") == "yes")
    # crude geneset percentile: rank of confirmatory GGE circadian beta vs all phenotypes
    betas = sorted(float(r["beta"]) for r in enrich if r.get("beta") not in (None, "NA"))
    gge = next((float(r["beta"]) for r in enrich if r["phenotype"] == "generalized_gge"), None)
    pct = (100.0 * sum(1 for b in betas if b <= gge) / len(betas)) if (betas and gge is not None) else None

    kpis = {
        "primary_verdict": primary.get("verdict", "—"),
        "geneset_percentile": f"{pct:.0f}th" if pct is not None else "—",
        "causal_status": "run MR (real data) — pending",
        "estimable_cells": f"{n_est}/{len(estimable)}" if estimable else "—",
    }
    status = [{"phenotype": r["phenotype"], "trait": "circadian_core",
               "neff": r.get("n_eff", "—"), "h2z": r.get("expected_h2_z", "—"),
               "status": "estimable" if r.get("estimable") == "yes"
                         else "not estimable — insufficient power"}
              for r in estimable]
    return {"kpis": kpis, "status": status,
            "primary_contrast": primary, "aim3": aim3[0] if aim3 else {}}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="results/report.html")
    args = ap.parse_args(argv)

    results = Path(args.results)
    dash = build_dashboard(results)

    # dashboard.json next to the static dashboard + in results
    (results / "dashboard.json").write_text(json.dumps(dash, indent=2) + "\n")
    dash_dir = Path("dashboard")
    if dash_dir.exists():
        (dash_dir / "dashboard.json").write_text(json.dumps(dash, indent=2) + "\n")

    p = dash["primary_contrast"]
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>epicirc results</title></head><body>
<h1>Circadian × epilepsy-type — results summary</h1>
<p><b>NOTE:</b> generated from the SYNTHETIC pipeline unless real data was supplied.</p>
<h2>Primary contrast (focal vs GGE, circadian core)</h2>
<ul>
<li>delta = {p.get('delta','—')} (SE {p.get('se_delta','—')}), rho = {p.get('rho','—')}</li>
<li>p = {p.get('p','—')} &rarr; <b>{p.get('verdict','—')}</b></li>
<li>min detectable difference = {p.get('min_detectable_delta','—')}</li>
</ul>
<p>KPIs: {dash['kpis']}</p>
<p>See dashboard/index.html for the interactive view; claims constrained by the adversarial review.</p>
</body></html>"""
    Path(args.out).write_text(html)
    print(f"[report] dashboard.json + {args.out} written; verdict={p.get('verdict','—')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
