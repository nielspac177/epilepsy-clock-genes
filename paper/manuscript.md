---
title: "Circadian gene common-variant architecture and its causal relationship to focal versus
  generalized epilepsy: a summary-statistics genetic study"
authors: [Rolston Lab]
target_journal: TBD (e.g. Brain, Epilepsia, or Neurology Genetics)
reporting_guideline: STREGA (genetic association) + MR-STROBE (for Aim 2)
status: PRE-DATA SKELETON — hypotheses and methods are pre-registered; Results/Discussion are
  placeholders filled only after the pipeline runs on locked data.
---

## Abstract
*(Structured; ≤250 words. Confirmatory scope = focal vs GGE. Subtype findings labelled
exploratory. No surgery/severity claims — Aim 3 is within-epilepsy pharmacoresistance.)*

## Introduction
- Circadian rhythm and seizures: chronoepileptology; GGE (esp. JME) seizures phase-locked to
  sleep/wake; generalized-seizure patients ~5× more likely to be late chronotypes.
- Prior circadian-gene work in epilepsy is animal expression or small candidate-gene studies
  (<100 patients, null); no systematic subtype-resolved genetic-architecture analysis.
- **Gap & aims** (one estimand per aim; "circadian" is not used as a shared construct):
  - Aim 1 — does circadian core-gene common-variant involvement differ between focal and GGE?
  - Aim 2 — is a circadian exposure *causally* related to epilepsy type, and does the signal
    localize to clock biology (cis) rather than sleep behaviour or metabolic/mood confounding?
  - Aim 3 — is circadian genetic burden associated with within-epilepsy pharmacoresistance?

## Methods
### Data
ILAE Consortium 2023 GWAS (focal, GGE, subtypes; epiGAD); chronotype/insomnia/sleep-duration/
napping GWAS; drug-resistant-epilepsy GWAS; MDD/BIP/ADHD; BMI/T2D; cis-eQTL/pQTL (GTEx/brain).
European ancestry primary. All summary-statistics; no individual genotypes.

### Pre-registration (before data lock)
Estimable-cell table (expected h2 Z from Neff), MR-power table (Brion/mRnd), feasibility table
(#instruments, F, h2 Z); frozen, hashed gene sets; single primary gene-set definition; SESOI for
equivalence tests; global multiplicity ledger. See `docs/adr/0003`.

### Aim 1 — gene-set involvement
MAGMA competitive analysis with base covariates **plus brain-expression and constraint (LOEUF/pLI)
gene covariates**; covariate-matched empirical null (1000+ sets); stratified-LDSC enrichment +
standardized τ*. Focal-vs-GGE difference via a **shared-control-corrected case-case GWAS** (or
control-resampling permutation), with GGE down-sampled to focal Neff. Robustness: leave-one-gene-
out, leave-known-epilepsy-out, leave-ion-channels-out, window sweep, LDSC-SEG.

### Aim 2 — causal triangulation
Behavioural-trait MR (IVW/Egger/weighted-median/PRESSO; **CAUSE** for correlated pleiotropy)
**and** cis-clock-gene MR with colocalization; **MVMR** on BMI, T2D, psychiatric liability;
overlap measured per pair (bivariate LDSC intercept) with **MRlap/CAUSE** where non-negligible;
winner's-curse shrinkage; **index-event-bias correction** (Slope-Hunter/CWLS) for type outcomes;
**network mediation** reporting total/direct/indirect with role-swap sensitivity; leave-clock-out
bridging. See `docs/adr/0004`.

### Aim 3 — pharmacoresistance
Type-conditioned DRE analysis (within-type / mtCOJO on focal-GGE liability), PGx-locus adjusted;
single pre-specified circadian-burden operationalization. See `docs/adr/0006`.

### Statistics
Bonferroni on the confirmatory family (focal, GGE, focal-vs-GGE); IHW on exploratory cells;
non-significance reported via TOST equivalence with minimum detectable difference.

## Results
*(placeholder — auto-populated from `results/` via `paper/tables/` and `paper/figures/`.)*

## Discussion
*(placeholder — claims constrained by docs/adversarial_design_review.md §"claims to soften".)*

## Limitations
- Confirmatory scope focal vs GGE only; subtypes exploratory/underpowered.
- Summary-statistics design; no individual-level PRS case-case classification.
- European ancestry; portability untested.
- No surgical-outcome data; Aim 3 is pharmacoresistance, with residual type-confounding.
- Behavioural sleep exposures index composite traits; causal "clock biology" claims require cis
  concordance.

## Data & code availability
All inputs are public (epiGAD + listed GWAS); pipeline at `github.com/<rolston-lab>/epicirc`,
reproducible via `snakemake` with pinned environment and checksummed inputs.
