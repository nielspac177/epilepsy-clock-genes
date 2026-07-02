# Adversarial Design-Hardening Report

*Produced by a 45-agent adversarial review (5 attack dimensions × verify × synthesize).
27 threats survived adversarial verification. This document is the authoritative list of design
changes; `config/analysis_params.yaml` and the ADRs implement them.*

## Five structural failure modes (executive summary)

1. **The causal engine (Aim 2) is not identified as designed.** mtCOJO/genomic-SEM "removes"
   psychiatric liability on the untestable assumption it is a *confounder*, not a *mediator* — the
   choice determines the sign/existence of any "direct circadian effect." MR-side sample overlap
   and winner's curse are uncorrected (the LDSC-intercept guard protects only rg). Behavioral sleep
   exposures do not localize causation to the molecular clock; metabolic confounders (BMI/T2D)
   are never conditioned out.
2. **The "differs by type" contrast is collider-contaminated.** Focal-vs-GGE is case-only;
   circadian genotype plausibly acts on seizure timing → EEG capture → classification, manufacturing
   a type difference with zero etiologic effect — invisible to LDSC intercept, PRESSO, Steiger.
3. **The design is powered for ~2–3 phenotypes, not 7 subtypes.** Effective N (shared ~52k
   controls): JAE ~2.4k, CAE ~5.5k, JME ~12.8k, GGE ~26k, focal ~50k. The h2 z>4 gate structurally
   deletes most subtype rg cells; MR into small subtypes runs near type-I power; reverse-MR/mediation
   instruments do not exist for JME/CAE/JAE (0–2 GWS loci).
4. **Gene-set enrichment localizes to a gene LIST, not to circadian biology.** MAGMA conditions on
   gene size/density but not brain expression or constraint (LOEUF/pLI); clock genes are long,
   brain-expressed, constrained, so "circadian" enrichment can be generic brain-gene signal. Random
   negative controls cannot reject this.
5. **Multiplicity accounting is mis-specified.** The FDR denominator counts structurally unestimable
   cells (diluting real power); gene-set-definition choice is an unmodeled degree of freedom.

## What the review VALIDATED (do not over-engineer)

- **BH-FDR stays valid under positive dependence**; the shared-control case-case *difference* test
  is **conservative, not anti-conservative** (confirms the fix in `stats/heterogeneity.py`).
- **LDSC jackknife block-difference already captures shared-control covariance** when the difference
  is formed *within* each leave-one-block-out replicate — no extra machinery needed for the rg arm.
  (The explicit `rho` correction is needed for the **MAGMA** case-case arm only.)
- MAGMA already conditions on gene length/SNP density (drop that sub-claim; keep brain-expr/constraint).
- No core clock gene sits in a long-range-LD/inversion or MHC region (NR1D1 at 17q21.1 is ~5 Mb from
  the 17q21.31 inversion) — the inversion/MHC critique was refuted.

## Confirmed threats → mitigations (ranked)

| # | Threat | Sev | Mitigation |
|---|--------|-----|-----------|
| 1 | Mediation identification failure (mtCOJO assumes confounder role) | Critical | Report **total / direct / indirect** effects (network MR) with role-swap sensitivity; never report mtCOJO as "the true psychiatric-independent effect" |
| 2 | Focal-vs-GGE collider / index-event bias | Critical | Publish DAG; **index-event correction (Slope-Hunter/CWLS)**; anchor causal claims on **subtype-vs-control**; ascertainment-matched positive controls (age-of-onset, EEG-access) |
| 3 | "Why circadian" label leakage (behavioral sleep ≠ clock biology) | Critical | Add **cis-eQTL/pQTL MR on core clock genes + colocalization (coloc/SMR-HEIDI)**; require behavioral-MR AND cis-clock concordance before "circadian biology" language |
| 4 | Gene-set = gene-LIST not biology (brain-expr/constraint confound) | Critical→Major | **Covariate-matched null** (length, SNP count, LD, GTEx brain-expr decile, LOEUF/pLI; 1000+ sets); MAGMA **conditional model with gene covariates** |
| 5 | h2 gate deletes rg grid; FDR denominator fiction | Critical | Pre-compute expected h2 z from Neff; **pre-register estimable-cell family**; gated cells = "not estimable", never "null" |
| 6 | Null difference read as "no differential architecture" at ~0 power | Critical | **TOST equivalence** with pre-registered SESOI (|Δrg|≥0.2); report minimum detectable difference; 3-way verdict detected/equivalent/inconclusive |
| 7 | MR into small subtypes underpowered; suite worsens it | Critical | Pre-register **per-pair MR power** (Brion/mRnd); confirmatory MR only for powered pairs (focal-all, GGE-all) |
| 8 | Untreated MR sample overlap | Critical→cond. | Measure **bivariate LDSC intercept per pair**; if non-negligible use **MRlap/CAUSE** as primary; report mean/conditional F |
| 9 | Winner's curse in instrument selection | Major | Conditional-likelihood/FIQT shrinkage; prefer MRlap (joint overlap+WC) |
| 10 | Directional/correlated pleiotropy via excitability channels; metabolic confounders | Major | **CAUSE** primary; **MVMR on BMI+T2D+psychiatric**; estimator agreement ≠ pleiotropy control |
| 11 | MAGMA between-type diff test has no cross-phenotype covariance | Critical→Major | Build **shared-control-corrected focal-vs-GGE GWAS** and run one gene-set analysis, OR control-resampling permutation null |
| 12 | Beta difference confounded with differential h2/N | Major | Set-by-type interaction vs matched background; **down-sample GGE to focal Neff** |
| 13 | GO membership confounded by pleiotropy + annotation bias | Major | Core set primary; GO **experimental evidence codes only (exclude IEA)**; leave-known-epilepsy-out, leave-ion-channels-out |
| 17 | Live DoF over which gene-set definition is primary | Major | Name **single primary definition** + fixed hierarchy + promotion rule; version-pin GO release |
| 18 | Small core set → gene-level instability | Major | **Leave-one-gene-out jackknife**; window sweep 0/10/35 kb; LDSC-SEG orthogonal check |
| 19 | Aim 1 and Aim 2 estimate disjoint quantities | Major | One estimand per aim; **leave-clock-out MR** bridging analysis; present as triangulation |
| 20 | Aim 3 DRE not a clean severity/surgery proxy; type-confounded | Critical→Major | **Within-type / type-conditioned** DRE analysis; restate as **pharmacoresistance**, drop surgery/severity language |
| 24 | Aim 3 "circadian burden" undefined multiverse | Minor | Single pre-specified summary-stats operationalization (MAGMA gene-Z aggregation over core set) |

*(Threats 14–16, 21–23, 25 fold into the above; see git history for the full verbatim table.)*

## Claims that must be softened in the manuscript

- "Circadian biology is causally linked to epilepsy type" → only if behavioral-trait MR **and**
  cis-clock-gene MR concur and survive index-event correction + MVMR. Else: "genetically-instrumented
  sleep behavior is associated with…".
- "Differs across the 7 subtypes" → restrict to **focal vs GGE**; subtype statements = exploratory.
- "Independent of psychiatric comorbidity" → forbidden as a point claim; report total/direct/indirect.
- Any "did not differ" / "no differential involvement" → only if a TOST CI excludes the SESOI;
  else "inconclusive/underpowered" with the minimum detectable difference.
- "No genetic correlation" for gated subtypes → "not estimable — insufficient power."
- "Directionality confirmed" (Steiger) → "consistent with, not established by" where instruments are weak.
- "Predicts refractory/surgical disease" → "within-epilepsy pharmacoresistance"; drop surgery framing.
- Estimator agreement across IVW/Egger/median/PRESSO → NOT evidence against pleiotropy/overlap.
