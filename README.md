# epicirc — Circadian genetics × epilepsy type

A reproducible, summary-statistics study of whether the common-variant architecture of circadian
("clock") genes differs between **focal** and **generalized (GGE)** epilepsy, and whether circadian
biology connects to epilepsy in *time* (sleep/genetic overlap) or in cortical *space* (the epileptic
network's molecular makeup). Everything here runs from public GWAS; nothing needs individual-level data.

> **Scope, honestly.** Confirmatory claims are limited to **focal vs GGE**; the seven ILAE subtypes
> are underpowered and stay exploratory. There is no public GWAS of surgical outcome, so that aim is
> a documented stub, not a result. See `docs/adr/` and `docs/adversarial_design_review.md`.

## Status: analyses complete

The obtainable analyses are done, independently re-verified by an adversarial agent team, and
robustness-checked. Two things remain, and both need external inputs rather than more work:
LDSC-SEG tissue/cell-type heritability (the precomputed annotations were not retrievable here) and
the MEG frequency-band maps (need Connectome Workbench). Neither changes the conclusions.

## Results summary

**In time, not in space.** Circadian genetics is tied to generalized epilepsy through when the brain
sleeps, not through where the epileptic network sits in the cortex.

- **Clock genes are enriched in GGE, not focal.** LD-aware MAGMA competitive test P = 1.7×10⁻⁴
  (β = 0.67), focal null (P = 0.73), housekeeping set null, robust to dropping the strongest gene
  (*PER1*, P = 1.4×10⁻⁴). Specific to the 23-gene core oscillator; the broad GO circadian annotation
  is null. It did not replicate in FinnGen R12 (1,690 cases, power- and phenotype-limited).
- **GGE shares genes with short sleep, not with chronotype.** LDSC genetic correlation:
  GGE ~ sleep duration rg = −0.12 (p = 0.002), specific to GGE (focal null, p = 0.70); GGE ~
  chronotype null (p = 0.47). The positive control rg(GGE, focal) = 0.61 (p = 2×10⁻¹⁶) validates the
  pipeline. GGE SNP-heritability 0.091 (intercept 1.058). Sleep deprivation, a classic generalized-
  seizure trigger, now has a genetic correlate.
- **No shared causal variant at PER1.** Two-GWAS colocalization puts GGE and chronotype at *distinct*
  variants (PP.H3 = 0.9997), and a three-trait analysis (GGE + chronotype + PER1 brain eQTL) gives a
  0.97 posterior of three independent signals. Even the anchor gene is not a shared regulatory switch.
- **Behavioural sleep MR is a robust null.** Chronotype does not causally drive either epilepsy type.
- **The LGS network has a molecular fingerprint, but not a circadian one.** Across 68 neuromaps
  annotations, the Lennox–Gastaut EEG-fMRI network overlaps glutamatergic mGluR5 density (r = 0.54)
  and glucose metabolism; both survive controlling for the dominant cortical gradient *and* a change
  of parcellation (Schaefer-400). A GABA-A/benzodiazepine signal appears but fails the parcellation
  check, so it is reported as suggestive only. Clock-gene expression has its own modest cortical
  topography yet maps onto neither the LGS network (an informative null: 80% power for |r| ≥ 0.31)
  nor the mGluR5 signature.

Full tables, figures, and the adversarial guardrails: [`docs/RESULTS_final.md`](docs/RESULTS_final.md).
A submission-ready manuscript (`.docx`, `.pdf`, and an Overleaf `.zip`) is in `submission/`.

## Study design

![Methods overview](docs/figures/methods_figure.png)

***Figure 1. Analysis overview.*** Public summary statistics feed the analyses; findings are
colour-coded by status (green positive, grey ruled out, amber open, blue conclusion). Every stage was
checked by an adversarial agent team, which caught real bugs along the way (a FinnGen tab-parsing
error; an unseeded spatial null) before they reached a result.

## What was actually run

- **Enrichment.** LD-aware MAGMA competitive gene-set analysis (1000G EUR, conditioning on gene
  size, density, sample size, MAC), a covariate-matched top-SNP null, and robustness checks
  (leave-one/k-out, in-body-only, median statistic, multi-seed, control gene sets); a shared-control
  focal-vs-GGE difference test.
- **Heritability and correlation.** LDSC SNP-h² and genetic correlation of GGE/focal with chronotype
  and sleep duration (belowlab LDSC, eur_w_ld_chr, w_hm3 merge-alleles), plus baselineLD partitioned
  h². Full UK Biobank sleep GWAS (Jones 2019, Dashti 2019) fetched from the GWAS Catalog.
- **Causality and mechanism.** Two-sample MR of behavioural sleep exposures, cis-eQTL MR,
  colocalization (coloc.abf plus a new two-GWAS and three-trait moloc extension) at PER1, and SuSiE
  fine-mapping.
- **Imaging.** Clock-gene AHBA expression and the LGS network against a 68-map neuromaps battery,
  with spin and variogram nulls, gradient-partial and Schaefer-400 robustness, and a power analysis.

## Layout
```
config/            pre-registered gene sets, data manifest, analysis params
src/epicirc/       data · geneset · magma · ldsc · mr (incl. coloc_multi) · stats · viz · render
scripts/           LDSC, coloc, neuromaps battery/dominance, figures, submission build
tests/             pytest (pure-Python core, no external tools); coloc_multi 8/8 is one
docs/              FINDINGS_*, RESULTS_final.md, ADRs, adversarial design review
paper/             IMRaD manuscript + figures/tables appendix
submission/        manuscript.docx · manuscript.pdf · Overleaf .zip
```

## Running it

The pure-Python statistical core (MR, coloc/moloc, harmonization, matched-null, FDR) is unit-tested
and needs nothing external:
```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e ".[dev]"
PYTHONPATH=src .venv/bin/python -m pytest
```
The LDSC and neuromaps steps run in two local venvs kept **off** any iCloud-synced folder (iCloud
eviction otherwise stalls imports for minutes):
```bash
uv venv ~/.epicirc-venv   && uv pip install --python ~/.epicirc-venv/bin/python "numpy<2" pandas scipy bitarray rich rich-argparse xopen modified-logger
uv venv ~/.epicirc-neuro  && uv pip install --python ~/.epicirc-neuro/bin/python "numpy<2" "scipy<1.13" nibabel nilearn matplotlib brainsmash neuromaps abagen
```
`numpy<2` matters: the belowlab LDSC fork calls `float(array)`, which numpy 2 removed. Scripts under
`scripts/` point at these venvs.

## Reproducibility
- Gene sets frozen and hashed before outcome analysis (`results/geneset_lock.json`).
- External inputs recorded in `config/traits_manifest.yaml`; nothing large is committed.
- Every headline number regenerates from pinned inputs and was re-derived independently by the
  verification team; figures rebuild from `scripts/make_figures.py`.
