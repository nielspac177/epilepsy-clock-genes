# Methods

*Reported against STREGA (genetic association) and STROBE-MR (Mendelian randomization). All analyses
use publicly available summary statistics; no individual-level genotypes were accessed. Ancestry is
European throughout, and cross-ancestry portability is not claimed. The full pipeline, pre-registered
gene set, and per-analysis verdicts are in the repository.*

## Data sources
Outcome data were the ILAE Consortium on Complex Epilepsies 2023 European meta-analyses for
all-epilepsy, focal epilepsy, generalized (GGE), and seven subtypes (epiGAD; GRCh37; effect allele
= Allele1; per-marker effective sample size and Zhu-derived β/SE provided). Effective sample sizes,
read directly from the files, were approximately 65,100 (all-epilepsy), 43,800 (focal), and 23,400
(GGE); subtypes ranged from ≈6,600 (JME) to ≈2,600 (JAE). Behavioural exposures were morning
chronotype (Jones 2019; GCST007565), sleep duration and short/long sleep (Dashti 2019; GCST007561,
GCST007559/60), and insomnia (Hammerschlag 2017; GCST004695) — all European, GRCh37. Expression
instruments were eQTLGen blood cis-eQTLs (N = 31,684) and PsychENCODE prefrontal-cortex cis-eQTLs
(N = 1,387, hg19). The gene model was GENCODE v19 (hg19; 19,347 protein-coding genes). The ILAE
outcome contains no UK Biobank participants, so overlap with the UK-Biobank-based sleep exposures is
minimal and biases MR estimates toward the null.

## Pre-registered circadian gene set
A set of 23 core clock genes spanning the transcription–translation feedback loop (ARNTL, ARNTL2,
CLOCK, NPAS2, PER1–3, CRY1/2, NR1D1/2, RORA/B/C, CSNK1D/E, FBXL3, DBP, NFIL3, TIMELESS, BHLHE40/41,
CIART) was specified and hash-frozen before any outcome was analysed, with the rationale documented
per gene. Genome-wide coordinates (GRCh37) were obtained from Ensembl and GENCODE.

## Aim 1 — Gene-set enrichment (LD-robust, covariate-matched)
For each protein-coding gene, all SNPs within a ±50 kb window were assigned to that gene, and the
single largest χ² (top SNP) was taken as an LD-robust gene-level statistic — one approximately
independent value per gene rather than many linkage-correlated SNPs. The circadian set's mean
top-SNP χ² was compared with 10,000 random gene sets matched on gene length and SNP count (deciles),
giving a competitive empirical p-value that is calibrated against the pipeline's own false-positive
rate. A 95% confidence interval on the enrichment ratio was obtained by bootstrapping the 23-gene
set. This design controls the two dominant confounds of a naïve SNP-level screen — within-gene
linkage disequilibrium and gene size / SNP density. A full LD-aware gene model (MAGMA) and
stratified LD-score regression were not run because the required binaries and reference panels were
unavailable at analysis time; the matched-null competitive test is the substitute and is reported
as such.

## Aim 2 — Mendelian randomization
**Instrument selection and harmonization.** Instruments were exposure variants at p < 5×10⁻⁸,
distance-clumped at ±1 Mb (an approximation used because no linkage-disequilibrium reference panel
was available), matched to the outcome by chromosome:position (not rsID, to avoid identifier
mismatch), and harmonized to the exposure effect allele; palindromic and allele-incompatible
variants were dropped. Instrument counts and Cochran's Q are reported per analysis.

**Estimators.** The primary estimator was inverse-variance weighting. Sensitivity estimators were
MR-Egger (with an intercept test for directional pleiotropy), the weighted median, an MR-PRESSO-style
outlier-corrected estimate, and Steiger-filtered IVW, following STROBE-MR guidance that a single
estimate is insufficient. The three instrumental-variable assumptions (relevance, independence,
exclusion restriction) are addressed respectively by genome-wide-significant selection, the low
sample overlap noted above, and the pleiotropy/heterogeneity diagnostics.

**cis-instrument MR.** For each clock gene, the strongest blood cis-eQTL (eQTLGen) was used as a
single instrument and the Wald ratio computed against each epilepsy outcome, with the exposure
effect on standardized expression approximated as Z/√N (eQTLGen reports Z, not β with allele
frequency). This is a screen, reported alongside the colocalization analysis that supersedes it.

## Aim 2c — Colocalization
Colocalization used coloc.abf with Wakefield approximate Bayes factors (priors p₁ = p₂ = 1×10⁻⁴,
p₁₂ = 1×10⁻⁵). The method is sign-agnostic (it uses Z², so the absence of an allele column in the
PsychENCODE eQTL file is immaterial); the epilepsy arm used ILAE β and SE directly, and the
expression arm used the eQTL p-value with minor-allele frequency (from ILAE) and N. Posterior
probabilities PP.H0–H4 are reported per gene, with PP.H4 the probability of a shared causal variant
and PP.H3 the probability of distinct (linked) variants.

## Variant annotation
The lead GGE variant in each clock-gene region was annotated with Ensembl VEP (GRCh37) for its most
severe consequence and mapped gene, to distinguish genuine clock-gene signals from co-located
neighbours.

## Statistical reporting
Enrichment significance is the matched-null empirical p; per-gene association is reported as nominal
χ²/p without genome-wide correction beyond the gene-set test. MR effects are reported as β with 95%
confidence intervals (β ± 1.96·SE). Colocalization is reported as posterior probabilities. No
individual-level polygenic scoring or case–case classification was performed, as only summary
statistics were available.

## Reproducibility
Inputs are referenced by accession and checksum; every result table regenerates deterministically
from pinned inputs. The analysis code is unit-tested and the study design was subjected to a
structured adversarial review before analysis (documented in the repository).
