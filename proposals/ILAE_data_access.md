# Data-access / collaboration brief — ILAE Consortium on Complex Epilepsies

**Working title.** Individual-level testing of circadian clock-gene contribution to genetic
generalized epilepsy: polygenic, subtype-resolved, and locus-level analyses.

**Requesters.** Rolston Lab (UCSF). *Contact: [PI], [email].*
**Directed to.** ILAE Consortium on Complex Epilepsies (enquiries: Prof. Sam Berkovic / Karen
Oliver); with reference to the ILAE Big Data Commission's data-sharing framework.

## Summary
Using the Consortium's publicly released 2023 European **summary statistics**, we found that common
variation in core circadian ("clock") genes is enriched for association with **genetic generalized
epilepsy (GGE)** but not focal epilepsy (competitive, LD-robust, covariate-matched enrichment ratio
1.53, 95% CI 1.05–2.04; empirical p = 1.7×10⁻³), anchored by a genome-wide-significant intronic
**PER1** variant (p = 5.8×10⁻¹⁰). Behavioural-sleep Mendelian randomization was null and brain-eQTL
colocalization did not implicate steady-state expression, pointing to a genetic, gene-local effect
whose mechanism is unresolved. Summary statistics cannot take this further; individual-level data
can.

## What summary statistics cannot do (and why we are writing)
1. **Case–case polygenic classification** — whether a circadian gene-set / PER1 polygenic score
   discriminates GGE from focal epilepsy at the individual level.
2. **Subtype-resolved analysis** — JME, CAE, JAE and GTCS-alone are underpowered in summary form
   (effective N ≈2,600–6,600); individual-level joint modelling and covariate adjustment would
   clarify the JME-specific hypothesis motivated by its morning-myoclonus phenotype.
3. **Locus-level dissection of PER1** — conditional/fine-mapping and, where available, splicing-QTL
   or coding annotation to test the non-eQTL mechanism our colocalization implies.
4. **Phenotype linkage** — if any contributing sites hold EEG (generalized spike-wave, photo-
   paroxysmal response) or imaging, testing circadian burden against these endophenotypes.

## Data requested
Individual-level imputed genotypes and harmonized phenotype labels (all-epilepsy, focal, GGE, and
the seven subtypes), age, sex, and ancestry principal components, for the European ILAE 2023
cohort — under the Consortium's standard data-use agreement and ethics framework, accessed via a
secure/managed environment as the Consortium prefers. Any available EEG/imaging endophenotypes
would be valuable but are not required for the core analyses.

## Analysis plan (brief)
Pre-registered polygenic and gene-set-burden models (circadian core set + PER1 region) for
case–case (GGE vs focal) and subtype-vs-control contrasts, with ancestry-PC and site adjustment,
appropriate multiple-testing control, and the same adversarial-verification workflow used in the
preliminary study. All code shared with the Consortium.

## What we bring
A complete, unit-tested, reproducible pipeline; a pre-registered gene set and analysis plan; a
documented adversarial design review; and the preliminary results above as rationale. We would
welcome full collaboration and would follow Consortium authorship and governance norms.

## Requested next step
A brief exploratory call to determine the appropriate access route (managed-access analysis,
collaborative sub-study, or Consortium membership) and the governance/ethics requirements.
