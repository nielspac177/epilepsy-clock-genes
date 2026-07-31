# Circadian genetics × epilepsy — consolidated results

*Every number below was independently recomputed by an adversarial verification team (5-agent
workflow); values reproduced to <10⁻¹⁴ relative error unless noted. Local analysis venvs were
rebuilt off iCloud to run LDSC and neuromaps reproducibly.*

## Headline

The circadian–epilepsy link is **genetic and temporal, not spatial-molecular**. Generalized
epilepsy (GGE) shows circadian gene-set enrichment (MAGMA) and a **specific genetic correlation
with short sleep duration**, but **no shared causal variant** with chronotype at PER1 and **no
circadian signature** in the cortical map of the Lennox–Gastaut (LGS) epileptic network. That LGS
network is not molecularly featureless, though: it carries a **receptor-specific excitation/
inhibition fingerprint** (glutamatergic mGluR5 + GABA-A/benzodiazepine) that survives control for
the dominant cortical gradient.

---

## 1. Genetics

### 1.1 SNP heritability (LDSC, univariate eur_w_ld_chr)

| Trait | h² (obs) | SE | Intercept | λ_GC | Note |
|:--|--:|--:|--:|--:|:--|
| GGE | 0.091 | 0.004 | 1.058 | 1.29 | minimal confounding |
| Focal | 0.091 | 0.004 | 1.054 | — | |

*The baselineLD 97-category partitioned run gives a model-dependent total (0.576) that is **not**
the reportable h²; only its per-category enrichments are used.*

### 1.2 Genetic correlation (LDSC rg) — the new positive finding

| Pair | rg | 95% CI | p | Verdict |
|:--|--:|:--:|--:|:--|
| **GGE ~ focal** *(positive control)* | **0.61** | 0.46, 0.75 | **2×10⁻¹⁶** | control **passes** ✓ |
| **GGE ~ sleep duration** | **−0.12** | −0.20, −0.04 | **0.0022** | **significant, GGE-specific** |
| GGE ~ chronotype | 0.02 | −0.04, 0.08 | 0.47 | null |
| Focal ~ sleep duration | −0.03 | −0.16, 0.11 | 0.70 | null |
| Focal ~ chronotype | 0.05 | −0.07, 0.17 | 0.43 | null |

`GGE ~ sleep-duration = −0.12` (shorter sleep ↔ higher GGE liability) survives Bonferroni across
the four primary sleep-trait tests (α = 0.0125) and is **specific to GGE** (focal null, p = 0.70).
The *behavioural* circadian axis (chronotype) shows no genetic overlap; the *homeostatic* axis
(sleep duration) does. Sleep deprivation is a classic GGE seizure trigger — this gives that a
genetic correlate.

> **Guardrails (adversarial):** rg is not causal. The `GGE ~ focal` control shares ILAE controls
> (cross-trait intercept z = 32) — it validates the pipeline but is **not** an independent-sample
> correlation. All sleep-trait pairs are overlap-free (UKB sleep GWAS vs ILAE; gcov_int ≈ 0).

![rg forest](../figures/fig_rg_forest.png)

### 1.3 Colocalization at PER1 (17p13.1)

Two-GWAS coloc and a coloc-consistent 3-trait moloc (`src/epicirc/mr/coloc_multi.py`, 8/8 unit
tests; re-implemented from scratch by the audit team → identical PPs):

| Test | Result | Meaning |
|:--|:--|:--|
| coloc(GGE, chronotype) | **PP.H3 = 0.9997**, PP.H4 = 4.6×10⁻⁶ | both associated, **distinct** variants |
| coloc(GGE, PER1-eQTL) | PP.H4 = 0.011 | no shared variant |
| coloc(chronotype, PER1-eQTL) | PP.H4 = 0.015 | no shared variant |
| **moloc(GGE, chronotype, PER1-eQTL)** | **PP(all-shared) ≈ 0**, PP(none-shared) = **0.97** | three independent signals |

Even at the PER1 anchor gene, epilepsy risk, morningness and PER1 expression are three independent
association signals — the circadian–epilepsy link is not a single shared regulatory variant.

---

## 2. Imaging — the LGS network fingerprint (neuromaps, 68 maps)

The LGS EEG-fMRI network was correlated (Spearman, spin-test nulls, BH-FDR within family) against
68 neuromaps annotations projected into each map's native space by registration fusion.

### 2.1 A-priori excitation/inhibition receptor family (pre-registered)

| Receptor map | r | p_spin | FDR | **partial r \| gradient** | p_spin(\|grad) |
|:--|--:|:--:|--:|--:|:--:|
| **mGluR5 (ABP688)** ×3 tracers | 0.54 | <0.001 | 0.0015 | **0.49** | **<0.001** ✓ |
| GABA-A/BZ (Ro15-4513) | 0.42 | <0.001 | 0.0015 | 0.28 | 0.016 ✓ |
| GABA-A/BZ (flumazenil, dukart) | 0.37 | 0.004 | 0.0048 | 0.31 | 0.019 ✓ |
| GABA-A/BZ (flumazenil, norgaard) | 0.21 | 0.067 | 0.067 | — | *n.s. (disclosed)* |
| glucose metabolism (CMRglc) | 0.48 | <0.001 | 0.004 | 0.45 | <0.001 ✓ |
| myelin (T1w/T2w) *— gradient control* | −0.27 | <0.001 | — | 0.14 | 0.25 ✗ |

**The key rigor test:** LGS's single strongest correlate is the principal FC gradient (r = 0.70),
and 36/68 maps are FDR-significant — the fingerprint of a broadly gradient-aligned map. So each
receptor was re-tested with the **principal gradient partialled out**. mGluR5, both GABA-A/BZ
tracers, and metabolism **remain spin-significant** after gradient control; **myelin does not**
(it was gradient-driven). The LGS network therefore carries a *receptor-specific* E/I signature
(glutamatergic + GABAergic), not merely a gradient echo.

![LGS network](../figures/fig_lgs_surface.png)
![GABA-A/BZ (flumazenil)](../figures/fig_gaba_surface.png)

### 2.2 Circadian is absent from this fingerprint

Clock-gene AHBA expression and clock-gene developmental trajectories do **not** align with the LGS
map (dual spatial nulls; prior `FINDINGS_spatial.md`). So the LGS network has a rich molecular
identity (mGluR5, GABA-A, metabolism, transmodal gradient) — circadian genetics simply is not part
of it.

> **Guardrails (adversarial):** p_spin = 0.001 is the 1000-permutation floor → reported as
> "<0.001". The norgaard flumazenil tracer is non-significant; GABA-A/BZ evidence is 2-of-3
> tracers. Receptor overlaps are correlational (spatial), now shown gradient-independent.

---

## 3. Multiple-testing ledger (headline tests)

| # | Test | Statistic | p / FDR | Verdict |
|--:|:--|:--|:--|:--|
| 1 | MAGMA circadian oscillator in GGE | P = 1.7×10⁻⁴ | — | positive (prior) |
| 2 | rg(GGE, sleep duration) | −0.12 | p = 0.0022 (Bonf ✓) | **positive** |
| 3 | rg(GGE, chronotype) | 0.02 | p = 0.47 | null |
| 4 | rg(GGE, focal) control | 0.61 | p = 2×10⁻¹⁶ | control ✓ |
| 5 | coloc(GGE, chronotype) @PER1 | H4 = 5×10⁻⁶ | — | null (H3) |
| 6 | moloc @PER1 | none-shared 0.97 | — | null |
| 7 | LGS × mGluR5 (gradient-adj) | 0.49 | p < 0.001 | **positive** |
| 8 | LGS × GABA-A/BZ (gradient-adj) | 0.28–0.31 | p = 0.016–0.019 | positive |
| 9 | LGS × clock-gene expression | — | dual-null | null |
| 10 | Behavioural MR (chronotype→GGE) | β = −0.05 | p = 0.10 | null (prior) |

---

## 4. Triangulation

- **Enrichment (+):** circadian oscillator over-represented in GGE association.
- **Genetic correlation:** GGE ↔ short sleep duration **(+)**; GGE ↔ chronotype (null).
- **Locus mechanism:** PER1 coloc/moloc **null** — no shared causal variant.
- **Behavioural MR:** chronotype → GGE **null**.
- **Cortical space:** clock ↔ LGS **null**; but LGS ↔ mGluR5/GABA-A **(+, gradient-independent)**.

→ **Circadian genetics acts on epilepsy in *time* (sleep-homeostatic genetic overlap, oscillator
enrichment), not in cortical *space*; the LGS network's molecular identity is an
excitation/inhibition receptor signature, into which circadian expression does not map.**

---

## Limitations
- rg reflects genetic overlap, not causation; sleep-duration effect is modest (rg −0.12).
- Chronotype is UKB-only (Jones 2019, 23andMe portion access-restricted) — ample for rg/coloc.
- MEG delta/theta band maps require Connectome Workbench (unavailable) — the one blocked a-priori
  comparison; documented, not hidden.
- Receptor overlaps are spatial correlations (now gradient-adjusted), not molecular measurements in
  epileptic tissue.
