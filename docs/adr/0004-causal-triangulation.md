# ADR-0004: Causal triangulation for Aim 2 (the "why circadian")

**Status:** accepted · **Date:** 2026-07-02

## Context
Review failure modes #1, #2, #3. Behavioural sleep-trait MR licenses only "sleep behaviour is
causal", not "the molecular clock is causal" (lead instruments sit in metabolic/mood loci, e.g.
FTO, MTNR1B). mtCOJO mediation assumes psychiatric liability is a confounder, not a mediator —
that assumption alone can create or destroy a "direct effect". MR onto a case-only "type" outcome
opens an index-event (collider) path via seizure-timing → classification.

## Decision
Claim "circadian **biology** is causally linked" only if **all** hold:
1. **Behavioural-trait MR** (chronotype, insomnia, sleep duration, napping) — IVW primary; Egger,
   weighted median, MR-PRESSO; **CAUSE primary for correlated pleiotropy**.
2. **cis-eQTL/pQTL MR on core clock genes** (ARNTL, CLOCK, PER1-3, CRY1/2, NR1D1/2, DBP, NFIL3)
   with **colocalization** (coloc / SMR-HEIDI). Behavioural and cis-clock evidence must concur.
3. **MVMR** conditioning jointly on **BMI, T2D, and psychiatric liability** (MDD/BIP/ADHD).
4. **Overlap handling:** bivariate LDSC intercept measured per pair; if non-negligible, **MRlap /
   CAUSE** become primary. Winner's-curse shrinkage (conditional-likelihood / FIQT) reported.
5. **Index-event-bias correction** (Slope-Hunter / CWLS) for any MR onto "type"; causal claims
   anchored on **subtype-vs-control** outcomes, not case-vs-case.
6. **Mediation** reported as **total / direct (MVMR) / indirect (two-step network MR)** with a
   role-swap sensitivity (confounder vs mediator vs collider). mtCOJO is a labelled sensitivity
   scenario only, never "the true psychiatric-independent effect".
7. **Bridging:** leave-clock-out MR tests whether the signal is carried by polygenic behavioural
   background rather than the clock.

## Consequences
Estimator agreement across IVW/Egger/median/PRESSO is explicitly **not** counted as pleiotropy or
overlap control. Confirmatory MR runs only for adequately powered pairs (focal-all, GGE-all).
