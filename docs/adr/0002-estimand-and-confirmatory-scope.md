# ADR-0002: Estimands and confirmatory-vs-exploratory scope

**Status:** accepted · **Date:** 2026-07-02 · **Supersedes** the plan's "7-subtype" framing.

## Context
The adversarial review (docs/adversarial_design_review.md, failure mode #3) showed the data are
powered for ~2–3 phenotypes, not 7 subtypes. Effective N with shared ~52k controls: JAE ~2.4k,
CAE ~5.5k, JME ~12.8k, GGE ~26k, focal ~50k. A 7-subtype confirmatory claim is not identifiable.

## Decision
- **One confirmatory contrast:** focal vs GGE, exposure = circadian core set.
- **Confirmatory phenotypes:** focal, GGE (all-epilepsy for context).
- **All per-subtype analyses are exploratory** and labelled as such in title and abstract.
- Each analysis aim estimates **one named estimand**; the bare word "circadian" is never used as a
  shared label across Aim 1 (gene-set involvement) and Aim 2 (causal effect of a circadian
  exposure). They are presented as **triangulation on distinct estimands**, not one construct.
- The focal-vs-GGE contrast is a **case-only** contrast; its estimand is stated as "difference
  conditional on diagnosed type" and carries the index-event-bias caveat (see ADR-0004).

## Consequences
Honest scope; the confirmatory family is small (3 cells → Bonferroni). Underpowered cells are
reported as "not estimable", never "null" (see ADR-0003).
