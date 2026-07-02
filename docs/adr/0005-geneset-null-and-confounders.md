# ADR-0005: Gene-set enrichment — null model and confounders

**Status:** accepted · **Date:** 2026-07-02

## Context
Review failure mode #4: MAGMA's competitive test conditions on gene size, density, and LD, but not
on **brain expression** or **constraint** (LOEUF/pLI). Core clock genes are long, broadly
brain-expressed, and constrained, so a "circadian enrichment" can be generic brain-gene signal.
Random (or random-brain) negative controls cannot reject that null.

## Decision
1. **Covariate-matched empirical null:** 1000+ gene sets matched to the circadian core on set size,
   gene length, SNP count, LD score, GTEx brain-expression decile, and LOEUF decile; report the
   circadian set's empirical percentile against this null.
2. **MAGMA conditional model** adding brain-expression and constraint gene covariates
   (`--gene-covar`).
3. **Definitional-validity reruns:** leave-known-epilepsy-genes-out, leave-ion-channels-out; GO set
   restricted to **experimental evidence codes (exclude IEA)**; gene-Z vs PubMed annotation-count
   check. "Circadian-specific" language only if signal survives all.
4. **Stability:** leave-one-gene-out jackknife for every significant set; window sweep (0/10/35 kb);
   orthogonal LDSC-SEG check.
5. **Type difference:** never a z-test of independent per-trait betas. Build a **shared-control-
   corrected focal-vs-GGE case-case GWAS** and run one gene-set analysis, OR a control-resampling
   permutation null; **down-sample GGE to focal Neff** to rule out power-driven differences.

## Consequences
`src/epicirc/stats/heterogeneity.py` (with the `rho` shared-control correction) applies to the
MAGMA case-case arm; the LDSC rg arm uses the within-replicate jackknife difference, which already
absorbs the covariance.
