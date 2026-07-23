# Tier 1–3 hardening: circadian core-oscillator enrichment in GGE

This document records the Tier 1→3 programme run to stress-test, replicate, and mechanistically
probe the headline finding (circadian gene-set enrichment specific to genetic generalized epilepsy).
Every result was produced on real ILAE3 / FinnGen R12 / GTEx data. **No result is protected**; the
negative replication and the pleiotropic-locus caveat are reported as prominently as the positives.

## Headline (updated)

The common-variant burden of the **23-gene core circadian oscillator** is enriched in **GGE**
(not focal epilepsy), and the enrichment survives the gold-standard LD-aware test:

| Test | GGE | Focal | Note |
|---|---|---|---|
| **MAGMA competitive (LD-aware, gene-density-conditioned)** | **P = 1.7×10⁻⁴** (β 0.67) | **P = 0.73** (β −0.11) | primary; housekeeping null 0.86 (GGE) |
| Matched-null top-SNP (±50 kb) | ratio 1.53, emp_p 1.9×10⁻³ (median of 50 seeds) | 0.98, p=0.54 | screen |
| In-body only (flank 0) | ratio 1.74, emp_p 1.7×10⁻³ | 0.97, p=0.57 | co-location-free |
| Rank-based (scale-free) | percentile 0.62, p=0.033 | 0.56, p=0.27 | cross-cohort-comparable |
| Difference test (GGE−focal) | Δ 0.545, 2-sided p 3.1×10⁻³ (50-seed median) | — | shared-control-robust |

**PER1** is the strongest single gene (MAGMA gene p = 7.9×10⁻⁷, genome-wide significant, rank
24/16,640; the windowed lead SNP rs2585398 is intronic to PER1 but is shared identically with the
co-located HES7/VAMP2/TMEM107 — "genome-wide significant" is a locus-level fact, and which gene is
causal is unresolved, see §3). **The set signal does not rest on PER1 alone, and is robust under the
primary (MAGMA) test:** dropping PER1 leaves P = 1.4×10⁻⁴, and dropping PER1+PER3 leaves P = 2.3×10⁻³.
(The matched-null top-SNP *screen* looked more fragile — leaving out the top 2 fell to p≈0.065 — but
that statistic over-weights a few co-located genes; MAGMA's gene-body assignment is the reliable
test and shows the signal is genuinely distributed, with PER3, ARNTL, NR1D2, NPAS2, NFIL3, RORB all
in the top ~15% of the genome.) Note the MAGMA/in-body tests use **22** genes: BHLHE40 has no SNP in
the 1000G EUR reference panel and is dropped.

## Tier 1 — internal hardening (no new data)

**1.1 Leave-one-out.** No single gene is necessary: all 23 one-gene drops stay emp_p<0.05 (worst =
PER1 drop, 0.021). **Leave-k-out in the matched-null screen:** dropping the top 2 (PER1+CSNK1E,
p=0.065) or top 3 (+PER3, p=0.157) is non-significant. But this fragility is *specific to the
top-SNP screen*, whose window statistic over-weights the co-located genes (CSNK1E's top SNP is in
KCNJ4; see 1b). **Under the primary MAGMA test the enrichment is robust to the same removals**
(drop PER1: P=1.4×10⁻⁴; drop PER1+PER3: P=2.3×10⁻³), because MAGMA assigns SNPs to gene bodies and
does not credit CSNK1E with the KCNJ4 signal. The honest reading: the *screen* is carried by ~2–3
genes; the *LD-aware test* is not.

**1.2 Difference test.** GGE enrichment significantly exceeds focal (Δ=0.545, p=3.1×10⁻³ multi-seed
median), robust across ±10–100 kb windows (Δ p from 1.3×10⁻³ to 9.5×10⁻³). Construction audited by
an independent statistician: valid competitive null, type-I ≤ nominal (not liberal).

**1.3 Controls.** Positive control (epilepsy genes) enriched in the matched-null screen (GGE p=0.040,
focal p=0.017); **housekeeping null in both** (0.14, 0.82) — FPR calibrated. Under strict MAGMA the
positive control is ns (0.12; expected — Mendelian genes carry weak common-variant GGE signal),
housekeeping cleanly null (0.86).

**1.4 Gene-set definition.** The signal is **specific to the 23 core TTFL genes**: the broad GO:0007623
set (185 genes) is only borderline (1.08, p=0.058) and **not GGE-specific** (Δ p=0.72); the 165
peripheral (non-core) circadian genes are **null in GGE** (1.02, p=0.32). So this is the core
molecular oscillator, not circadian biology broadly.

**Window sensitivity.** Robust across ±10/20/50/100 kb; *stronger at tighter windows* (flank10 ratio
1.72 vs flank100 1.41) — the signature of a gene-local effect, not a co-location artifact.

## Tier 1b — the co-location attack and its resolution

The adversarial audit's most serious concern: 15/23 clock top-SNPs lie in the flank, and the two
biggest window contributors are neighbour-driven (CSNK1E top-SNP inside **KCNJ4**; PER3 top-SNP inside
**CAMTA1**). Three independent lines resolve it:

1. **In-body-only enrichment** (flank 0, excludes all neighbour SNPs): GGE ratio **1.74**, p=1.7×10⁻³,
   focal null, Δ p=1.4×10⁻³ — *stronger* without the flank. An artifact would weaken, not strengthen.
2. **MAGMA gene-level** (assigns SNPs to gene bodies): **CSNK1E becomes null (p=0.40)** — MAGMA does
   not credit it with the KCNJ4 signal — yet the set stays enriched (P=1.7×10⁻⁴), carried by PER1,
   PER3, ARNTL, NR1D2, NPAS2 (all genuine gene-body signal).
3. **Median (robust) statistic**: GGE p=4.0×10⁻³ — not an artifact of a few outliers.

## Tier 2 — external validation

**2.1 FinnGen R12 replication — NEGATIVE (power-limited).** Independent Finnish biobank
(GE: 1,690 cases / 484,703 controls). The clock set does **not** replicate, and two independent
statistics agree on clean data: mean-χ² ratio 0.85 (emp_p 0.99) and scale-free rank percentile 0.40
(emp_p 0.99); GE_STRICT and FE are also null. **This is reported as a genuine non-replication.**
Context: FinnGen GE is ~2× smaller (effective) than ILAE GGE, registry-phenotyped (ICD-code
generalized epilepsy, broader/noisier than ILAE's expert-classified GGE), and a population isolate;
a positive-control check shows it detects the strongest established GGE loci only weakly (VRK2
χ² 79→14.8, SCN1A →13.2, STX1B →12.6; PER1 →5.7, gene-level 1-df p≈0.017, and the ILAE lead SNP
rs2585398 itself is non-significant). **Interpretation:** consistent with insufficient power +
phenotype/ancestry difference rather than clear refutation, but the finding is not independently
replicated and must be stated as such.

> **Reproducibility note / correction.** An earlier version of this analysis reported a spurious
> "inflation tail" (1,266 FinnGen genes χ²>50) and corrupted enrichment ratios (GE 0.30). These were
> a parsing bug in `gene_stats.py` — a whitespace split on the TAB-delimited FinnGen file collapsed
> empty rsID fields (indels), shifting columns and turning `sebeta` into `beta` for those rows. Fixed
> (tab split); caches regenerated. Correctly parsed, **FinnGen GE is a clean, uninflated GWAS** (max
> gene χ²≈31, zero genes >50). The non-replication conclusion is unchanged; the numbers above are the
> corrected values, and the mean-χ² and rank tests now agree.

**2.2 MAGMA (LD-aware) — CONFIRMED.** Real MAGMA v1.10 with the 1000G EUR panel, conditioning on
gene size, gene density, sample size, inverse MAC and their logs. **CIRCADIAN_CORE P = 1.7×10⁻⁴**
(β 0.67, 22 genes — BHLHE40 absent from the panel) in GGE, **null in focal (P = 0.73, β −0.11)**;
housekeeping null in both. Robust to dropping PER1 (P=1.4×10⁻⁴) and PER1+PER3 (P=2.3×10⁻³). PER1
is GGE-specific under MAGMA (gene p 7.9×10⁻⁷ vs focal 0.41). This is the confirmation the README
listed as pending, directly rebuts the gene-size / gene-density / within-gene-LD confounds, and
establishes GGE-specificity under the gold-standard method. Now the primary result. *(P = 1.7×10⁻⁴
is MAGMA's competitive p; a simpler covariate-adjusted competitive test reproduces the effect size
β≈0.68 at p≈2–3×10⁻³ — robust in direction and significance, exact exponent is model-specific.)*

## Tier 3 — mechanism at the PER1 anchor

**3.1 Locus content.** rs2585398 is genuinely **intronic to PER1** (chr17:8054860 ∈ PER1 body). But
17p13.1 is gene-dense: the association spans ~150 kb (HES7, VAMP2, AURKB, CTC1, PFAS all elevated by
LD from a single peak centred on PER1). The apparent shoulder near **ALOX12B** (windowed χ²=29,
−79 kb) is **not** an independent second peak — ALOX12B has essentially no in-body signal (χ²=0.69),
and χ² decays monotonically from the PER1/HES7/VAMP2 peak. PER1 is the best positional candidate
(lead intronic to it) but fine-mapping is required to resolve it against its co-located neighbours
(the genuine competitor is VAMP2, not ALOX12B — see 3.2).

**3.2 QTL mechanism — PER1 sQTL (GTEx v8), plus co-located brain eQTLs.**
- rs2585398 is a strong **PER1 splice-QTL** in peripheral tissues (Artery_Tibial p=6×10⁻²⁴; ~10
  tissues; GTEx v8) — matching the coloc result that it is *not* a steady-state eQTL. Splicing
  mechanism supported, but observed peripherally.
- The same SNP is an **eQTL for VAMP2** (synaptobrevin-2, a synaptic-vesicle SNARE) in **10 brain
  regions** (Putamen p=3×10⁻¹⁴, Hypothalamus 9×10⁻¹³, Hippocampus 1×10⁻¹⁰; GTEx v10), and for CTC1
  and others. **Caveat, stated carefully:** VAMP2 is a *biologically plausible* competing candidate
  (synaptic gene with a brain eQTL), but this rests on plausibility, not statistics — there is **no
  demonstrated colocalization** with the GGE signal, rs2585398 is a QTL for ~11 genes at this locus
  (the CTC1 whole-blood eQTL, p=7×10⁻²⁷, is far stronger than any VAMP2 effect), and VAMP2's brain
  eQTL effects are modest (|NES| 0.11–0.40). No PER1 **brain** sQTL was detected, but that is
  uninformative about true absence (GTEx v8 sQTL, low brain power, only cerebellum covered at this
  locus). The splice-vs-expression contrast is also drawn across two GTEx releases (v8 sQTL vs v10
  eQTL). **Net:** the causal gene at the PER1 locus is genuinely unresolved (PER1 / VAMP2 / CTC1 /
  HES7 all in play); do not assert the brain action is VAMP2 rather than PER1. The *set-level*
  circadian enrichment is robust to dropping PER1 (MAGMA P=1.4×10⁻⁴), so this locus-level ambiguity
  qualifies the mechanism story without touching the gene-set result.

**3.3 JME / GTCS-only — exploratory, null (underpowered).** Neither subtype shows enrichment
(JME ratio 1.02, p=0.42; GTCS 1.03, p=0.37); PER1 χ² is much lower (JME 6.1, GTCS 12.6 vs GGE 38.4).
Subtypes are severely underpowered (JME N≈6,600) and the GGE signal aggregates across subtypes, so
the null is consistent with power. The circadian signal is a **pan-GGE** property, not JME-specific.

## What changed vs the pre-Tier-1 write-up

- **Primary result is now MAGMA** (P=1.7×10⁻⁴), not the matched-null screen.
- The claim is narrowed to the **core oscillator** (not broad circadian annotation).
- Report **multi-seed** p-values (medians), not single lucky seeds; relabel the bootstrap CI.
- Add the **co-location resolution** (in-body + MAGMA) and the honest **leave-k-out** fragility.
- Add the **FinnGen non-replication** prominently.
- Add the **VAMP2 pleiotropy** caveat at the PER1 locus.
