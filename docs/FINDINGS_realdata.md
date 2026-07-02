# Real-data preliminary findings

Two independent analyses on the real ILAE3 European data both point to **GGE-specific circadian
involvement**. These are FIRST LOOKS with tool-free/approximate methods; the confirmatory pipeline
(MAGMA + stratified-LDSC + robust/cis MR) is still required (see caveats).

## 1. Circadian gene-set enrichment (Aim 1 first look)
`results/real_circadian_enrichment.tsv` — mean-χ² of SNPs in ±50 kb of 23 core clock genes vs genome.

| Phenotype | Enrichment ratio | z | Reading |
|---|---|---|---|
| **GGE** | **1.49** | 25.3 | strongly enriched |
| Focal | 0.91 | −5.1 | not enriched |
| All-epilepsy | 0.99 | −0.7 | null |

Per-gene audit (`results/real_circadian_pergene_generalized_gge.tsv`): the GGE signal is **broad
across ≥6 independent chromosomes** — **PER1 genome-wide significant (P=5.8×10⁻¹⁰)**, CSNK1E
(9×10⁻⁷), PER3 (4×10⁻⁶), ARNTL/BMAL1 (8×10⁻⁶), RORA (2×10⁻⁵), NPAS2 (4×10⁻⁵). The **same genes are
null in focal**, controlling for gene-size/brain-expression confounds on the difference.

## 2. Two-sample MR of sleep exposures (Aim 2) — behavioural traits
`results/mr_real/*.tsv` — 5 estimators, distance-clumped (±1 Mb), CHR:POS-matched, harmonized.
Exposures fetched GRCh37 from GWAS Catalog: chronotype (Jones 2019), sleep duration / short / long
sleep (Dashti 2019). Insomnia (Hammerschlag 2017) had too few instruments (underpowered, 113k);
the larger Jansen insomnia is access-restricted (23andMe).

**Robustness is decisive — the nominal IVW hits do NOT survive sensitivity analysis** (β [p]):

| Exposure → GGE | IVW | weighted median | MR-PRESSO | Steiger | verdict |
|---|---|---|---|---|---|
| sleep duration | −0.18 [0.04] | **+0.18** [0.23] | −0.08 [0.38] | +0.01 [0.95] | **not robust** (sign flips; 4 outliers; Q=164) |
| short sleep | +0.92 [0.03] | −0.11 [0.86] | +0.23 [0.59] | −0.39 [0.57] | **not robust** (unstable; borderline Egger pleiotropy) |
| chronotype | −0.05 [0.10] | −0.04 [0.48] | −0.04 [0.24] | −0.02 [0.53] | directionally stable (protective) but **ns** |
| long sleep | +1.96 [0.10] | — | — | — | uninterpretable (3 instruments) |

**Conclusion:** no behavioural sleep/circadian trait has a **robust** causal effect on GGE or focal.
Nominal IVW signals were driven by heterogeneous/outlier instruments (high Cochran Q). Chronotype
shows a consistent-but-non-significant protective trend for GGE > focal. **Interpretation:** the
strong GGE circadian gene-set enrichment (Result 1) is **not explained by causal sleep behaviour**,
which points toward the molecular clock / shared genetic architecture — motivating cis-clock-gene MR.

## Caveats (why these are not yet claims)
- Enrichment all-SNP z is **LD-inflated**; needs MAGMA (gene-based, LD-aware) + covariate-matched
  null (brain-expression/constraint) and stratified-LDSC (ADR-0005).
- MR uses **distance-clumping** (LD approximation, no panel) and only IVW/Egger; the nominal
  sleep-duration→GGE hit has **high Cochran Q** (instrument heterogeneity) → needs weighted-median,
  MR-PRESSO, and **CAUSE** before trust (ADR-0004).
- **Multiple testing**: 4 MR tests + enrichment; nominal p-values are not corrected.
- Behavioural exposures ≠ clock biology: **cis-eQTL/pQTL MR on core clock genes + colocalization**
  is required to attribute causality to the molecular clock, plus **MVMR** on BMI/T2D/psychiatric.
- **Insomnia** (Jansen 2019) full stats are access-restricted (23andMe); GWAS Catalog hosts only a
  no-full-p version. Needs a DUA or the UKB-only release.

## Next steps
1. Robust MR (weighted median, MR-PRESSO, CAUSE) + reverse MR + Steiger on the two exposures.
2. cis-clock-gene MR + colocalization; MVMR conditioning on BMI/T2D/psychiatric.
3. Install MAGMA (macOS) + LD panel → real gene-set analysis with matched null; stratified-LDSC rg.
4. Fetch daytime-napping + short/long-sleep exposures (GCST007559/007560) to complete the panel.
