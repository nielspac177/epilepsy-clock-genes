# Concept proposal — ENIGMA-Epilepsy imaging-genetics sub-study

**Working title.** Circadian clock-gene common variation and thalamocortical structure in genetic
generalized epilepsy.

**Proposers.** Rolston Lab (UCSF). *Contact: [PI], [email].*

## Background and rationale
The generalized epilepsies have a longstanding, quantified relationship to the circadian system:
generalized tonic-clonic seizures on awakening occur within two hours of waking in ~90% of cases,
patients with generalized epilepsy are disproportionately evening chronotypes (36% vs 11% in focal
epilepsy), and correcting circadian misalignment can markedly reduce seizures. Despite this, no
study has tested whether the *genetics* of the molecular clock relate to brain structure in
generalized epilepsy.

In preliminary work using the ILAE 2023 European GWAS, we found that common variation in 23 core
circadian ("clock") genes is enriched for association with **genetic generalized epilepsy (GGE)**
(competitive, LD-robust, covariate-matched enrichment ratio 1.53, 95% CI 1.05–2.04; empirical
p = 1.7×10⁻³) but **not focal epilepsy** (ratio 0.98). The signal is anchored by a genome-wide-
significant intronic variant in **PER1** (p = 5.8×10⁻¹⁰). Genetically-instrumented sleep *behaviour*
showed no robust causal effect, and brain-eQTL colocalization did not implicate steady-state
expression — indicating a genetic, gene-local effect whose brain correlate is unknown.

## Gap ENIGMA is uniquely positioned to fill
ENIGMA-Epilepsy is, to our knowledge, the only resource that pairs harmonized multi-site MRI with
merged genetic data at the scale needed (≈288 adults with IGE/GGE plus controls across 18 centres),
and it has established precedent for relating imaging phenotypes to epilepsy-risk-gene expression
axes and for multi-spectral diffusion mega-analyses in GGE.

## Hypothesis
Circadian clock-gene common-variant burden (and specifically the PER1 locus) is associated with
GGE-relevant **thalamocortical** structural phenotypes — thalamic volume, cortical thickness in
fronto-central regions, and structural-network integrity — with focal epilepsy serving as a
negative control (no association expected).

## Data requested
For existing ENIGMA-Epilepsy participants with both imaging and genotypes: harmonized subcortical
volumes, cortical thickness/surface area, and DTI/structural-network metrics; imputed genotypes at
the 23 clock-gene regions (or a pre-specified circadian gene-set score); diagnosis (IGE/GGE vs focal
vs control), age, sex, site. No raw imaging or individual genotypes need leave contributing sites —
the analysis fits the standard ENIGMA distributed/mega-analytic model.

## Analysis plan (brief)
1. Circadian gene-set burden score (and PER1 region score) per participant.
2. Association with thalamocortical ROIs and network metrics (linear models; site, age, sex,
   intracranial volume covariates), within GGE and, as a negative control, within focal epilepsy.
3. Pre-registered, multiple-testing controlled; sensitivity to gene-set definition; PER1 locus
   fine-mapping where genotype density allows.

## What we bring
A fully reproducible, unit-tested pipeline; a pre-registered gene set; and the summary-statistics
preliminary results above as the analytic and biological rationale.

## Ask
We propose to join the ENIGMA-Epilepsy working group and lead this analysis as a collaborative
sub-study, following ENIGMA authorship and data-use conventions. A two-page analysis plan and
pre-registration can be provided on request.
