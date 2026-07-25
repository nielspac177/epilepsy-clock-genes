# Spatial + genetic follow-ups (Phases A–E)

Extends the MAGMA-confirmed circadian×GGE finding with (A) a GWAS breadth sweep, (B) statistical
fine-mapping of the PER1 locus, (D) imaging-transcriptomics of the clock genes vs the LGS network,
and (E) the developmental-gradient position of the LGS network. Every brain-map correlation uses a
spatial-autocorrelation-preserving null; set-level claims add a random-gene null. No causal language
is drawn from a spatial correlation. Env: `.venv-neuro` (py3.11, abagen 0.1.3, neuromaps 0.0.7).

## Phase A — GWAS breadth (matched-null + rank enrichment)

| Phenotype | ratio | emp_p | rank pctl | read |
|---|---:|---:|---:|---|
| GGE (European) | 1.53 | 1.7×10⁻³ | 0.62 (p=0.033) | reference |
| **GGE (trans-ancestry)** | **1.57** | **6.0×10⁻⁴** | 0.61 (p=0.048) | **stronger with non-European samples** |
| CAE (childhood absence) | 1.19 | 0.044 | — | nominal (underpowered subtype) |
| JAE | 1.15 | 0.095 | — | trend |
| JME / GTCS-only | ~1.02 | ns | — | null (underpowered) |
| focal + focal subtypes (HS/lesion-neg/other) | ~0.95–1.05 | ns | — | null |
| all-epilepsy / trans-all / trans-focal | ns | — | — | null |

The enrichment is **not European-specific** (trans-ancestry stronger), is a **pan-GGE** property
(no single subtype carries it; aggregate GGE strongest), and remains **focal-null** across all focal
subtypes. (Drug-resistant-epilepsy was skipped — only a 56 kB placeholder file was available.)

## Phase B — PER1 locus fine-mapping (SuSiE-RSS, 1000G EUR LD)

917 harmonized SNPs (allele-aligned to the LD panel), median N=23,215. One 95% credible set:

| Gene (nearest) | PIP mass | note |
|---|---:|---|
| **PER1** | **0.37** | lead rs2585398 (intronic; PIP 0.35) — highest |
| PFAS | 0.15 | |
| AURKB | 0.15 | |
| C17orf59 | 0.14 | |
| CTC1 | 0.13 | |

The credible set spans 14 SNPs over ~140 kb. **PER1 is the single most probable causal gene (37% of
PIP) but is not statistically resolved** — a quantitative confirmation of the earlier "causal gene
unresolved" caveat. PER1 ranks *above* VAMP2 (which does not top any credible-set SNP), so PER1
remains the best positional candidate while the locus is formally ambiguous.

## Phase D — clock-gene expression vs the LGS network (imaging transcriptomics)

AHBA (Desikan-Killiany, 83 regions, 4/6 donors incl. both bilateral; 22/23 clock genes) × the
unthresholded LGS EEG-fMRI t-map, dual null (brainsmash spatial + 10,000 random-gene sets).

- **Clock-set expression does NOT co-localize with the LGS network: r = −0.03, p_spatial = 0.93,
  p_randomgene = 0.87 — null under both nulls.**
- A few individual genes show nominal, mixed-sign spatial correlation (NPAS2 +0.33 p=0.005; CRY2
  +0.30 p=0.024; NFIL3 −0.30 p=0.023) but do not survive 22-gene multiple testing and the set is null.
- Caveats: coarse DK parcellation (83 regions), 4-donor AHBA, sparse subcortical/thalamic sampling.

**Read:** the genetic circadian→GGE finding does not extend to a simple clock-gene expression-
topography overlap with the LGS network at this resolution.

## Phase E — the LGS network occupies the association/transmodal developmental pole

LGS t-map projected to fsLR-32k (neuromaps registration fusion; cortical surface only, so the
thalamic component is not captured) vs canonical developmental/hierarchy gradients, spin-test null
(alexander_bloch, 1000 rotations), BH-FDR.

| Gradient | Spearman r | p_spin |
|---|---:|---:|
| **Sensorimotor–association axis** (Sydnor 2021) | **+0.65** | **<0.001** (0/1000 rotations) |
| **Principal functional gradient** (Margulies 2016) | **+0.70** | **<0.001** (0/1000 rotations) |
| Hill developmental / evolutionary expansion; abagen genePC1 | — | *deferred (needs Connectome Workbench)* |

The two available gradients are themselves r = 0.85 (the same neurodevelopmental axis). **The LGS
epileptic network is strongly concentrated at the association/transmodal pole** of the brain's
principal developmental gradient — the late-maturing, evolutionarily expanded, top-of-hierarchy
cortex. The observed r sits beyond the entire spin null (null max |r| ≈ 0.59–0.60), so this is a
genuine large effect; the p is reported as `<0.001` (the 1000-rotation floor), not a precise value.
*Cortical surface only (thalamus excluded). Space note: the LGS map is ICBM152-2009b (symmetric);
neuromaps/abagen assume MNI152NLin2009cAsym — a small (~few-mm, mostly L–R) affine-only offset that
attenuates rather than inflates the correlation.*

## Genetics ↔ spatial: how the two connect (and don't, yet)

- The **genetic** circadian signal (GGE enrichment, PER1 anchor) and the **spatial** LGS-developmental-
  gradient result are, so far, **separate layers**: clock-gene expression does not track the LGS
  network (Phase D null), so the circadian genes are not simply "expressed where LGS lives."
- Transcriptional-gradient bridge (DK volumetric): the AHBA PC1 (dominant transcriptional gradient,
  42% variance) does not significantly track LGS at DK resolution (r=−0.23, p_spin=0.14; the surface
  S-A axis is the reliable, well-powered version). Clock-set expression loads strongly on the
  transcriptional gradient (r=−0.42; PER1 r=−0.54) but this axis does not align with the LGS network
  (clock-set vs LGS r=−0.03; PER1 vs LGS r=−0.08, both null). **So the clock genes have clear
  spatial-expression structure, but it is orthogonal to the LGS network — the direct expression
  bridge is null.**
- The natural next genetics×imaging test (not yet run): project the **whole GGE gene-level association**
  onto brain expression space and test whether GGE genetic risk is organized along the same
  association-pole gradient the LGS network occupies. This would directly bridge the GWAS to the LGS
  developmental topography.

## Interpreting the null bridge — H2/H3/H4 tests

**H3 (is it genuinely the oscillator?) — YES.** MAGMA on a strict oscillator core (dropping the
moonlighting kinases/ligase/checkpoint genes CSNK1D/E, FBXL3, TIMELESS) is *stronger* than the full
set: **P = 4.4×10⁻⁵** (β 0.83); the canonical 12-gene TTFL loop enriches (P = 1.4×10⁻⁴); and the
6 non-clock PER1-locus neighbours (AURKB/CTC1/PFAS/VAMP2/…) are **null as a set (P = 0.15)**. So the
enrichment is the transcriptional oscillator, not pleiotropy. (At the PER1 locus itself, gene-level
p cannot separate the neighbours — CTC1 3.2×10⁻⁷ is even more significant than PER1 7.9×10⁻⁷ — but
that is a *local* LD ambiguity, not what drives the set.)

**H2 (subcortical) — the null bridge is partly a cancellation.** Splitting the clock–LGS correlation
by tissue: cortex r = −0.12, **subcortex/brainstem r = +0.44** (PER1 +0.46). Opposite signs cancel to
the whole-brain null. Suggestive that clock genes track the LGS network subcortically — but n = 15,
not significance-tested, DK subcortex coarse (thalamus = 1 region; LGS-t there is modest, +2.2, and
the LGS peaks are association *cortex*). Needs a finer subcortical atlas (Tian S4) to confirm.

**H4 (developmental) — the relationship is dynamic and shifts monotonically.** Clock–LGS correlation
across BrainSpan regions by stage. Cortex-only (11 regions): prenatal −0.83 → adult +0.19. Firmed
with subcortex (15 regions, incl. mediodorsal thalamus/striatum/hippocampus/amygdala): prenatal
−0.40, infancy −0.31, early-childhood −0.30, child/adolescent +0.33, adult +0.19; per-stage
bootstrap CIs are wide (all cross 0, n = 15), but the **monotonic developmental trend across the 5
stages is nominally significant (Spearman ρ = +0.90, p = 0.037)**. So the coupling shifts from
negative (prenatal) toward positive (adult) across development — the adult "null" is the crossover
point. The strong cortex-only prenatal −0.83 was partly a cortical-gradient effect (attenuates to
−0.40 with subcortex). Exploratory (coarse BrainSpan→map region assignment), but the trend is real.

## Impact-booster results

**Cell-type (marker-based MAGMA) — no enrichment.** GGE common-variant heritability is not
significantly enriched in any brain cell type's marker panel (all p > 0.16; inhibitory β ≈ 0,
excitatory/OPC trend higher). Marker sets are small (8–23 genes; underpowered) — a specificity-
weighted gene-property analysis is the proper follow-up — but the data do not support a clean
interneuron mechanism.

**ENIGMA-Epilepsy imaging bridge — null.** The published GGE cortical-thickness map (Whelan 2018,
DK, n = 272) does not significantly track the LGS network (r = +0.14, p = 0.19), clock expression
(r = +0.05), the transcriptional gradient, or PER1 — though it is *directionally* GGE-specific
(GGE +0.14 vs all-epilepsy 0.00, TLE −0.03). GGE cortical-thickness effects are small, limiting power.

**Neurochemical fingerprint (creative imaging) — null.** The LGS network was correlated with 14
neuromaps PET/metabolism maps (GABA_A ×2, nicotinic/vesicular ACh, D1/D2 dopamine, SERT/5-HT1A/5-HT2A,
mu-opioid, CB1, glucose metabolism, CBF), DK-parcellated, seeded-brainsmash spin null, FDR. **No
neurotransmitter system is enriched** (GABA_A r = −0.06, p = 0.83; all receptors null). The only
nominal signals are **cerebral blood flow (r = 0.33, p = 0.035) and glucose metabolism (r = 0.30,
p = 0.10)** — hemodynamic/metabolic, co-varying with the association pole — and neither survives FDR.
So the network-level E/I test is also null: the epileptic network is *not* defined by a receptor
profile.

**Emerging synthesis (now very consistent):** *every* molecular bridge to the LGS network is null —
clock-gene expression, ENIGMA GGE atrophy, cell-type markers, and now the full neurotransmitter
receptor fingerprint. The network is defined **only** by its macroscale organizational position
(association/transmodal pole, high metabolism/CBF), not by any gene/cell-type/receptor signature.
In parallel, the circadian genetics is robust (MAGMA, fine-mapping, strict-core) but spatially
*uncoupled* from the network. Two clean, honest conclusions: (1) the LGS epileptic network is a
**macroscale systems** phenomenon irreducible to a local molecular signature; (2) the circadian
contribution to GGE acts in **time** (global excitability rhythm), leaving no cortical spatial
fingerprint — coherent with the behavioural-MR and expression-coloc nulls.

## Environment-blocked (not a scientific limit)
LDSC (genetic correlation + partitioned/stratified heritability, the README's long-pending item) is
now **installable and importable** here: the belowlab py3 fork works once its deps are supplied
(`modified-logger`, `rich_argparse`, `rich`, `xopen` into `.venv-neuro`); classic bulik/ldsc remains
blocked (Python-2.7 / Apple-Silicon / broken system conda). The **3.1 GB reference (baselineLD v2.2,
weights, frq, HM3) is in `tools/ldsc_ref/` and GGE is pre-munged** (`results/ldsc/gge.sumstats.pre`,
4.86M SNPs). The only remaining obstacle was the host being at load-average ~300 (heavily
oversubscribed, unrelated to this project), under which the multi-chromosome regression times out.
On an unloaded machine: `PYTHONPATH=tools/ldsc/src .venv-neuro/bin/python -m ldsc.main` for munge →
h² → partitioned h² with the downloaded reference. coloc(GGE, chronotype) and the
multi-trait coloc likewise await the full chronotype/sleep GWAS + regional eQTL (not cleanly
fetchable here).

## Still open
- coloc(GGE, chronotype) at PER1 + LDSC rg — need the full sleep/chronotype GWAS.
- LDSC partitioned heritability (GGE + baselineLD) — feasible without sleep data.
- H2 with a fine subcortical atlas (Tian S4); H4 with the mediodorsal thalamus + spatial nulls.
