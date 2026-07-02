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

## 2. Two-sample MR of sleep exposures (Aim 2 first look)
`results/mr_real/*.tsv` — IVW, distance-clumped (±1 Mb), CHR:POS-matched, harmonized.

| Exposure → outcome | IVW β | p | Egger intercept | Cochran Q |
|---|---|---|---|---|
| chronotype → GGE | −0.050 | 0.10 | ~0 (ns) | 233/98 (high) |
| chronotype → focal | −0.030 | 0.19 | ~0 (ns) | 112/97 |
| **sleep duration → GGE** | **−0.181** | **0.042** | ~0 (ns) | 164/51 (high) |
| sleep duration → focal | +0.091 | 0.17 | ~0 (ns) | 56/50 |

**Pattern:** longer sleep duration and morningness both trend **protective for GGE**, and the sleep
-duration effect **reverses sign in focal** — a GGE-specific, direction-consistent echo of the
enrichment result. Sleep duration → GGE is nominally significant.

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
