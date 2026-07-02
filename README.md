# epicirc — Circadian genetics × epilepsy type

A reproducible, summary-statistics pipeline testing whether the common-variant genetic
architecture of circadian ("clock") genes differs between **focal** and **generalized (GGE)**
epilepsy, and whether circadian biology is **causally** linked to epilepsy type — plus a
pharmacoresistance (morbidity proxy) analysis.

> **Scope, honestly.** Confirmatory claims are limited to **focal vs GGE**; the 7 ILAE subtypes are
> underpowered and treated as exploratory. There is **no public GWAS of surgical outcome**, so that
> aim is a future-ready stub, not a result. See `docs/adr/` and `docs/adversarial_design_review.md`.

## Scientific aims
- **Aim 1** — Differential circadian gene-set involvement by type (MAGMA competitive + covariate-
  matched null + stratified-LDSC), with a shared-control-corrected focal-vs-GGE contrast and TOST
  equivalence testing.
- **Aim 2** — Causal triangulation: behavioural sleep-trait MR **and** cis-clock-gene MR + coloc,
  MVMR on metabolic + psychiatric confounders, overlap-aware estimators (MRlap/CAUSE),
  index-event-bias correction, network mediation (total/direct/indirect).
- **Aim 3** — Within-epilepsy pharmacoresistance (drug-resistant-epilepsy GWAS), type-conditioned.

## Layout
```
config/            pre-registered gene sets, data manifest, hardened analysis params
src/epicirc/       data · geneset · magma · ldsc · mr · stats · viz
workflow/Snakefile reproducible DAG
tests/             pytest (pure-Python core: stats, harmonize, MR, gene-set) — no external tools
docs/adr/          architecture decision records
docs/adversarial_design_review.md   45-agent design critique + mitigations
paper/             IMRaD manuscript
dashboard/         accessible results dashboard
```

## Quick start
```bash
uv venv --python 3.12 .venv          # or conda env create -f environment.yml
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m pytest           # 23 tests, no external tools needed
snakemake -n                         # dry-run the full analysis DAG
```

The statistical core (heterogeneity with shared-control correction, BH/Bonferroni, allele
harmonization, IVW/MR-Egger, gene-set hashing/matched-null) is pure Python and unit-tested. LDSC,
MAGMA, PLINK, GCTA, and R (TwoSampleMR/CAUSE/MRlap) are invoked as pinned external tools for the
full runs (see `environment.yml`).

## Reproducibility guarantees
- Gene sets frozen + hashed before outcome analysis (`results/geneset_lock.json`).
- Pre-registration tables (`estimable_cells`, `mr_power`, `feasibility`) written before data lock.
- External inputs referenced by checksum in `config/traits_manifest.yaml`; nothing large committed.
- Every result regenerates deterministically from pinned inputs; software versions logged.
