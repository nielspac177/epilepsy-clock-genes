
\newpage

## Tables

**Table 1. LDSC genetic correlation of epilepsy with sleep traits.**

| Pair | rg | 95% CI | p |
|:--|--:|:--:|--:|
| GGE ~ focal (positive control) | 0.61 | 0.46, 0.75 | 2×10⁻¹⁶ |
| GGE ~ sleep duration | −0.12 | −0.20, −0.04 | 0.0022 |
| GGE ~ chronotype | 0.02 | −0.04, 0.08 | 0.47 |
| Focal ~ sleep duration | −0.03 | −0.16, 0.11 | 0.70 |
| Focal ~ chronotype | 0.05 | −0.07, 0.17 | 0.43 |

GGE, genetic generalized epilepsy. Unconstrained LDSC intercept; UK Biobank sleep GWAS have no
subject overlap with the ILAE epilepsy GWAS (cross-trait intercept ≈ 0). The GGE–focal control shares
ILAE controls and is not an independent-sample correlation.

**Table 2. Colocalization at the PER1 locus (17p13.1).**

| Test | Result | Interpretation |
|:--|:--|:--|
| coloc(GGE, chronotype) | PP.H3 = 0.9997; PP.H4 = 5×10⁻⁶ | both associated, distinct causal variants |
| coloc(GGE, PER1 eQTL) | PP.H4 = 0.011 | no shared variant |
| coloc(chronotype, PER1 eQTL) | PP.H4 = 0.015 | no shared variant |
| moloc(GGE, chronotype, PER1 eQTL) | P(all-shared) ≈ 0; P(none-shared) = 0.97 | three independent signals |

**Table 3. LGS network receptor overlap and its robustness.**

| Map | fsLR-32k (spin) | gradient-adjusted | Schaefer-400 (variogram) | Robust |
|:--|--:|--:|--:|:--:|
| mGluR5 (ABP688) | 0.54, <0.001 | 0.49, <0.001 | 0.36, 0.028 | yes |
| glucose metabolism (CMRglc) | 0.48, <0.001 | 0.45, <0.001 | 0.25, 0.005 | yes |
| GABA-A/BZ (flumazenil) | 0.37, 0.004 | 0.31, 0.019 | 0.06, 0.72 | no |
| GABA-A/BZ (Ro15-4513) | 0.42, <0.001 | 0.28, 0.016 | 0.17, 0.49 | no |
| cortical myelin (control) | −0.46, <0.001 | 0.14, 0.25 | — | no |

Values are Spearman r and p. p_spin reported at the 1000-permutation floor as "<0.001". Only mGluR5
and metabolism survive both gradient control and a change of parcellation and null family.

\newpage

## Figures

![**Figure 1.** The Lennox–Gastaut EEG-fMRI network (unthresholded t-map) projected onto the inflated
cortical surface.](figures/fig_lgs_surface.png)

![**Figure 2.** GABA-A/benzodiazepine receptor density (flumazenil PET) on the cortical surface, for
visual comparison with the LGS network.](figures/fig_gaba_surface.png)

![**Figure 3.** Forest plot of LDSC genetic correlations of GGE and focal epilepsy with sleep
duration and chronotype (rg with 95% CI; dashed line = no correlation).](figures/fig_rg_forest.png)

![**Figure 4.** The circadian × epileptic-network coupling flips sign across development. Rows are
BrainSpan developmental stages; left column is mean clock-gene expression at that stage, right column
the fixed adult LGS network, with the per-stage Spearman r. The relationship runs from prenatal
anticorrelation (clock high where the future network is low) toward adult alignment; the trajectory
panel shows the monotonic trend (cortex ρ = 0.80, firming to ρ = 0.90, p = 0.037 with subcortex).
Cortical surface at BrainSpan-region resolution; exploratory.](figures/fig_developmental_signflip.png)
