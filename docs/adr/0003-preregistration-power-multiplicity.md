# ADR-0003: Pre-registration, power gating, and multiplicity

**Status:** accepted · **Date:** 2026-07-02

## Context
Review failure modes #5 and #6: a data-dependent FDR family and "null = no effect" readings turn
an underpowered grid into false negatives dressed as findings; a live gene-set-definition degree
of freedom inflates false positives.

## Decision
Written **before any outcome is touched** (checksummed):
1. `results/estimable_cells.tsv` — expected h2 Z from Neff per phenotype×trait; the FDR family is
   fixed to *estimable* cells only.
2. `results/mr_power.tsv` — Brion/mRnd power + minimum-detectable-OR per exposure→outcome; only
   powered pairs are confirmatory.
3. `results/feasibility.tsv` — per-subtype #GWS instruments, mean F, h2 Z (gates reverse-MR/mediation).
4. Gene sets are hashed and frozen (`results/geneset_lock.json`); one **primary definition**
   (circadian core), a fixed secondary hierarchy, and a promotion rule.

**Error control:** small Bonferroni family for confirmatory cells (focal, GGE, focal-vs-GGE);
IHW (covariate = Neff / h2 Z) for exploratory cells; the gene-set-definition axis enters the
correction; an effective-number-of-tests (eigendecomposition of the estimate correlation matrix)
is reported so the ~28-cell grid is not over-read. A **global multiplicity ledger** spans Aims 1–3.

**Non-significance:** reported via **TOST equivalence** against a pre-registered SESOI
(|Δrg| ≥ 0.2) with the minimum detectable difference; verdicts are detected / equivalent /
inconclusive. "Did not differ" is forbidden unless the TOST CI excludes the SESOI.

## Consequences
BH/IHW remain valid under the positive dependence induced by shared controls (verified). Power is
not diluted by structurally unestimable cells.
