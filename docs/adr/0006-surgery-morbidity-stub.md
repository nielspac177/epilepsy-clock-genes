# ADR-0006: Morbidity/severity via pharmacoresistance proxy; surgery as a stub

**Status:** accepted · **Date:** 2026-07-02

## Context
The original goal included **response to epilepsy surgery**. No public GWAS of surgical outcome
exists; that analysis needs individual-level data the lab does not currently hold. The review
(threat #20) further showed the drug-resistant-epilepsy (DRE) GWAS is **type-confounded** and is
not a clean "severity/surgery" proxy.

## Decision
- Aim 3 estimand is restated as **within-epilepsy pharmacoresistance**, not severity or surgery.
- DRE analyses are **type-conditioned** (within-type, or mtCOJO conditioning DRE on focal/GGE
  liability) and adjusted for known pharmacogenomic loci (HLA, drug transporters); crude and
  conditioned results are both reported. rg omitted if DRE h2 Z < 4.
- "Circadian burden" for Aim 3 is a single pre-specified summary-stat operationalization
  (MAGMA gene-level Z aggregated over the core set); variants are labelled sensitivity.
- **Surgical-outcome analysis is a documented pipeline stub** (`src/epicirc/mr/` + a data schema)
  that accepts individual-level surgical/PRS data when a cohort (Epi25 / UKB / in-house) is
  obtained. It is not executed and is not claimed in the manuscript.

## Consequences
The manuscript drops "surgery / morbidity / severity" language for Aim 3 and frames it as
pharmacoresistance, with residual type-confounding flagged as a limitation.
