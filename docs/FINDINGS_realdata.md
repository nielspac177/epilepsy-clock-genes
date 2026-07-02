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

## 3. cis-clock-gene MR (Aim 2, molecular) — the informative one
`results/mr_real/cis_clock__*.tsv` — top blood cis-eQTL (eQTLGen) per clock gene → ILAE, Wald ratio.

| Gene → GGE | Wald | p | direction |
|---|---|---|---|
| **ARNTL (BMAL1)** | **+0.100** | **0.003** | higher expression → higher GGE risk |
| NR1D1 (REV-ERBα) | −0.159 | 0.072 | higher expression → lower GGE (BMAL1's repressor — coherent) |
| others (NPAS2, CRY1/2, PER3, RORA, ARNTL2) | — | ns | — |

Focal: nothing nominal (top CRY1 p=0.11). **ARNTL survives Bonferroni across the 8 testable genes
and is GGE-specific**; it was also a hit in the enrichment (Result 1). Only 8 genes had a cis-eQTL
present in ILAE; **PER1 (the enrichment lead) has no blood cis-eQTL** — likely brain-specific.

**Convergence:** the molecular clock (ARNTL/BMAL1, NR1D1) shows a GGE-specific causal signal
*exactly where behavioural sleep MR was null* — evidence that GGE circadian involvement operates at
the level of clock-gene regulation, not sleep behaviour.

## Caveats (why these are not yet claims)
- Enrichment all-SNP z is **LD-inflated**; needs MAGMA (gene-based, LD-aware) + covariate-matched
  null (brain-expression/constraint) and stratified-LDSC (ADR-0005).
- MR uses **distance-clumping** (LD approximation, no panel). Robust methods now run; behavioural
  hits did **not** survive (see Result 2).
- **cis-MR is a screen, not proof**: single blood cis-eQTL per gene, Z/√N effect approximation, **no
  colocalization** (ARNTL Wald could be LD-confounded by a neighbouring gene), blood≠brain.
- **Multiple testing**: nominal p-values across the MR grid are not globally corrected (ARNTL does
  survive Bonferroni within the 8-gene cis panel).
- Insomnia (Jansen 2019) full stats are access-restricted (23andMe); Hammerschlag UKB version is
  underpowered (0 usable instruments). Daytime sleepiness (Wang) not cleanly on the FTP.

## Next steps (priority reordered after these results)
1. **Colocalization** (coloc/SMR-HEIDI) on ARNTL & NR1D1 to rule out LD-confounding — the single
   most important check before the cis result means anything.
2. **Brain cis-eQTLs** (GTEx cortex / MetaBrain) for the clock genes — esp. PER1 (no blood eQTL);
   multi-instrument cis-MR where allelic heterogeneity allows.
3. Install MAGMA (macOS) + LD panel → LD-aware gene-based test with covariate-matched null;
   stratified-LDSC rg (replaces the LD-inflated enrichment screen).
4. MVMR conditioning on BMI/T2D/psychiatric; CAUSE for the behavioural exposures.
