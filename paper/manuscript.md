---
title: "Circadian clock-gene common-variant architecture is enriched in genetic generalized
  epilepsy but not focal epilepsy: a summary-statistics genetic study"
authors: [Rolston Lab]
target_journal: "Epilepsia / Brain Communications / Neurology: Genetics (TBD)"
reporting: STREGA (association) + MR-STROBE (Mendelian randomization)
status: DRAFT — numbers are from the real analyses in results/ and docs/FINDINGS_realdata.md.
  Confirmatory pipeline items still pending are flagged [PENDING].
---

## Abstract

**Background.** Generalized epilepsies, particularly juvenile myoclonic epilepsy, show a strong
clinical relationship to the sleep–wake cycle, and generalized-seizure patients are more likely to
be evening chronotypes. Whether the circadian ("clock") gene system contributes to the *common-
variant* genetic architecture of epilepsy, and whether it does so differentially by epilepsy type,
has not been tested systematically.

**Methods.** Using summary statistics from the ILAE Consortium on Complex Epilepsies (European
analysis; focal epilepsy, genetic generalized epilepsy [GGE], and all-epilepsy), we tested a
pre-registered set of 23 core clock genes for enrichment of common-variant association using an
LD-robust (top-SNP-per-gene) competitive test against 10,000 gene-length- and SNP-count-matched
random gene sets. We then interrogated causality with two-sample Mendelian randomization (MR) of
four behavioural sleep exposures (five estimators including MR-Egger, weighted median, MR-PRESSO
and Steiger filtering), *cis*-eQTL MR of clock genes (blood and brain), and formal colocalization
(coloc.abf) using PsychENCODE prefrontal-cortex expression.

**Results.** Circadian core genes were enriched for GGE association (mean top-SNP χ² 10.70 vs
matched-null 7.00; ratio 1.53; empirical p = 1.7×10⁻³), driven by a genome-wide-significant
intronic *PER1* variant (rs2585398, p = 5.8×10⁻¹⁰) with additional intronic signals at *ARNTL* and
*NPAS2*. The enrichment was **absent in focal epilepsy** (ratio 0.98, p = 0.54) and attenuated in
all-epilepsy (ratio 1.12, p = 0.15). No behavioural sleep exposure showed a robust causal effect on
either epilepsy type: nominal inverse-variance-weighted signals did not survive weighted-median,
MR-PRESSO or Steiger analysis. A blood *cis*-eQTL signal for *ARNTL* did **not** colocalize with GGE
in brain (PP.H4 = 0.01), and *PER1* — despite a strong regional GGE signal (PP.H2 = 0.89) — showed
no expression colocalization, indicating a non-eQTL (likely splicing/regulatory) mechanism.

**Conclusions.** Circadian clock-gene common variation is enriched specifically in genetic
generalized epilepsy, anchored by *PER1*, but is not explained by genetically-instrumented sleep
behaviour and does not colocalize with steady-state gene expression. Circadian involvement in GGE
operates at the level of genetic architecture; the precise molecular mechanism remains open.

## Introduction

Seizures in the generalized epilepsies are tightly coupled to the sleep–wake transition: juvenile
myoclonic epilepsy (JME) is defined in part by myoclonus on awakening, sleep deprivation is a
canonical trigger, and patients with generalized seizures are several-fold more likely than
controls to report an evening chronotype. These observations motivate the "chronoepileptology"
hypothesis that circadian biology shapes seizure liability. Yet prior genetic work linking the
molecular clock to epilepsy has been limited to animal models of clock-gene expression and small
candidate-gene studies (typically <100 patients) that were underpowered and largely null; a recent
review explicitly identified well-powered, subtype-resolved analysis as an open question.

The ILAE Consortium's third genome-wide association study, comprising 29,944 cases and 52,538
controls, revealed markedly different common-variant architectures for focal versus generalized
epilepsy and now provides the statistical power to test the circadian hypothesis directly. We asked
three questions, each with a distinct estimand: (1) is the common-variant burden in core clock
genes differentially enriched between focal and generalized epilepsy; (2) is a circadian exposure
*causally* related to epilepsy type; and (3) if so, does the signal localize to the molecular clock
rather than to sleep behaviour or co-located genes. We deliberately treat these as triangulation on
separate estimands rather than a single "circadian" construct.

## Methods

**Data.** Outcome summary statistics were the ILAE 2023 European analyses for all-epilepsy, focal
epilepsy and GGE (effect allele = Allele1; per-marker effective N; GRCh37). Behavioural exposures
(GRCh37, GWAS Catalog) were morning chronotype (Jones 2019), sleep duration and short/long sleep
(Dashti 2019); insomnia (Hammerschlag 2017, UK Biobank) was underpowered and the larger Jansen 2019
GWAS is access-restricted. *cis*-eQTLs were eQTLGen blood (N = 31,684) and PsychENCODE prefrontal
cortex (N = 1,387, hg19). The clock-gene set (23 genes spanning the core transcription–translation
feedback loop) was pre-registered and hash-frozen before analysis.

**Gene-set enrichment (Aim 1).** For each protein-coding gene (GENCODE v19, hg19; ±50 kb window)
we took the single top-SNP χ² as an LD-robust gene statistic. The circadian set's mean was compared
to 10,000 random gene sets matched on gene length and SNP count (deciles), yielding an empirical
competitive p-value. This design controls the two dominant confounds of a naïve SNP-level screen —
within-gene LD and gene size/SNP density. [PENDING: LD-aware MAGMA and stratified-LDSC once a
macOS binary and LD reference panel are available; these were blocked at analysis time.]

**Mendelian randomization (Aim 2).** Instruments were selected at p < 5×10⁻⁸ and distance-clumped
(±1 Mb; an LD-panel-free approximation), matched to the outcome by chromosome:position, and
harmonized to the exposure effect allele (palindromic and incompatible variants dropped). We report
IVW, MR-Egger (with intercept test), weighted median, an MR-PRESSO-style outlier-corrected estimate,
and Steiger-filtered IVW. *cis*-MR used the strongest cis-eQTL per gene as a single instrument
(Wald ratio). Colocalization used coloc.abf (Wakefield approximate Bayes factors; priors
p₁=p₂=1×10⁻⁴, p₁₂=1×10⁻⁵), sign-agnostic and therefore robust to the absence of allele codes in the
brain-eQTL file; the GWAS arm used ILAE β/SE directly and minor-allele frequencies from the ILAE
data.

**Rigor.** All analyses are European-ancestry. Sample overlap between UK-Biobank-based sleep
exposures and the ILAE outcome (which contains no UK Biobank) is minimal and biases MR toward the
null. The full pre-registration, adversarial design review, and honest per-analysis verdicts are in
the project repository.

## Results

**Circadian genes are enriched in GGE, not focal epilepsy.** The 23 core clock genes showed a mean
top-SNP χ² of 10.70 in GGE versus a matched-null expectation of 7.00 (ratio 1.53; empirical p =
1.7×10⁻³). The signal was anchored by a genome-wide-significant intronic variant in *PER1*
(rs2585398, p = 5.8×10⁻¹⁰) with suggestive intronic signals at *ARNTL/BMAL1* (p ≈ 8×10⁻⁶) and
*NPAS2* (p ≈ 4×10⁻⁵). The identical gene set was unenriched in focal epilepsy (ratio 0.98, p =
0.54) and only weakly enriched in all-epilepsy (ratio 1.12, p = 0.15), consistent with dilution of
a GGE-specific effect. Because the focal analysis uses the same genes, the null focal result argues
against a generic gene-size or brain-expression artefact driving the GGE enrichment.

**No robust causal effect of sleep behaviour.** Across four behavioural exposures and both epilepsy
types, nominal IVW estimates did not survive sensitivity analysis. For example, sleep duration →
GGE was nominally protective by IVW (β = −0.18, p = 0.042) but reversed sign under weighted median
(+0.18) and attenuated to null under MR-PRESSO (−0.08, p = 0.38, four outliers removed) and Steiger
filtering (+0.01, p = 0.95), with high instrument heterogeneity (Cochran Q = 164 on 51 df). Morning
chronotype showed a directionally consistent but non-significant protective trend for GGE (IVW
β = −0.05, p = 0.10) that was larger than for focal epilepsy. We therefore find no reliable evidence
that genetically-instrumented sleep behaviour causally alters epilepsy risk of either type.

**The GGE signal does not colocalize with clock-gene expression.** A blood *cis*-eQTL Wald estimate
suggested higher *ARNTL* expression raises GGE risk (p = 0.003), but this did **not** replicate as a
shared causal variant in prefrontal cortex: formal colocalization gave PP.H4 = 0.01 for *ARNTL* and
< 0.03 for every clock gene tested. For *PER1*, colocalization placed almost all posterior mass on a
GGE association *without* an expression signal (PP.H2 = 0.89), and VEP annotation localized the lead
variant to a *PER1* intron. Together these indicate that the *PER1* GGE signal is genuine and
gene-local but is **not** mediated by steady-state expression, implicating a splicing or regulatory
mechanism. We note that two suggestive set members (the signals nearest *CSNK1E* and *PER3*)
annotate to neighbouring non-clock genes (*KCNJ4*, a potassium channel; *CAMTA1*), a limitation of
window-based assignment.

## Discussion

Three lines of evidence converge on a specific and defensible conclusion: the common-variant
architecture of the circadian clock is enriched in genetic generalized epilepsy, not focal
epilepsy, and this enrichment is a property of genetic association rather than of sleep behaviour or
steady-state gene regulation. The result places a long-standing clinical observation — the
sleep-linked phenotype of the generalized epilepsies — on a genetic footing, and it does so with the
appropriate specificity: the focal-null contrast, the matched-null design, and the failure of
behavioural MR each guard against an obvious artefact.

Equally important is what we do **not** claim. Genetically-instrumented sleep behaviour is not a
demonstrable cause of either epilepsy type in these data, and no clock gene's expression colocalizes
with GGE. The initially attractive *cis*-MR signal at *ARNTL/BMAL1* did not survive brain
colocalization and has been retired. The strongest and cleanest signal, an intronic
genome-wide-significant *PER1* variant, is not an expression quantitative-trait locus, pointing to a
mechanism — plausibly splicing — that steady-state eQTL catalogues cannot capture. This is a
concrete, falsifiable target for functional follow-up.

## Limitations

The enrichment test, while LD-robust and covariate-matched, uses a top-SNP statistic rather than a
full LD-aware gene model (MAGMA) or partitioned heritability (stratified-LDSC), both pending on
tool/reference availability. Window-based gene assignment attributes some suggestive signal to
neighbouring non-clock genes. MR used distance-based clumping in place of formal LD clumping.
Colocalization used a modestly-powered single brain tissue with p-value-based eQTL effect
approximations, so an absence of colocalization for *PER1/ARNTL* reflects, in part, limited eQTL
power (posterior mass on H0/H2, not H3) and is inconclusive rather than a refutation. Analyses are
European-ancestry; portability is untested. We had no individual-level genotypes, so no per-person
polygenic classification, and no epilepsy-surgery-outcome data.

## Conclusion

Circadian clock-gene common variation is enriched specifically in genetic generalized epilepsy,
anchored by an intronic *PER1* signal, and is independent of genetically-instrumented sleep
behaviour and of steady-state clock-gene expression. The circadian contribution to GGE is genetic
and gene-local; defining its molecular mechanism — beginning with *PER1* splicing — is the natural
next step.

## Data and code availability

All input GWAS are public (ILAE/epiGAD; GWAS Catalog; eQTLGen; PsychENCODE). The full reproducible
pipeline, pre-registration, adversarial design review, and per-analysis verdicts are in the project
repository; every result table regenerates deterministically from the pinned inputs.
