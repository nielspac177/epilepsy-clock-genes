# Reporting standards → which tables & figures this paper needs

*Literature-grounded map of the tables/figures expected for a genetic-association + Mendelian-
randomization + colocalization study, and how ours satisfy them. Sources: STREGA (STROBE extension
for genetic association), STROBE-MR (JAMA 2021), coloc reporting conventions, MAGMA/gene-set
practice, and NEJM/JAMA in-table effect-size style.*

## The three governing guidelines
- **STREGA** — STROBE extension, 22 items; genetics-specific additions on **population
  stratification/ancestry, rationale for choice of genes/variants, statistical methods, multiple
  testing, and descriptive/outcome data** ([PLOS Med 2009](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000022)).
- **STROBE-MR** — 20 items/30 subitems; requires reporting the **three IV assumptions** (relevance,
  independence, exclusion-restriction), **instrument strength (F/conditional-F)**, **sample
  overlap**, and a **panel of sensitivity estimators** (Egger, weighted median, MR-PRESSO, Steiger)
  with estimates + CIs ([JAMA 2021](https://jamanetwork.com/journals/jama/fullarticle/2785494)).
- **coloc** — report **PP.H0–H4** per locus; H4 = posterior probability of a shared causal variant;
  high H3 = distinct variants (linkage). ([coloc convention](https://www.bioconductor.org/packages/release/bioc/vignettes/INTACT/inst/doc/INTACT.html))

## NEJM/JAMA table style (what the user asked for)
Effect estimates are shown **in-cell as `point (95% CI)`**, one row per comparison, with an
optional thin **inline graphical CI** column — i.e. a "table that contains the forest plot" rather
than a separate forest-plot figure. We render the inline CI with a fixed-scale unicode bar so it
reads on GitHub without image assets.

## Table set for THIS paper (maps guideline → our table)
| Table | Purpose | Guideline driver |
|---|---|---|
| **T1 Data sources** | every GWAS/eQTL input: N cases/controls, ancestry, build, accession, role | STREGA (data sources, ancestry, volume-of-data); STROBE-MR (data sources) |
| **T2 Gene-set enrichment** | circadian set per phenotype: obs, matched-null [95%], ratio (95% CI), empirical p | STREGA (statistical methods, multiple testing); MAGMA/gene-set practice |
| **T3 Lead clock-gene signals** | per-gene lead SNP, consequence, mapped gene, χ², P | STREGA (descriptive/outcome data, rationale for variants) |
| **T4 Mendelian randomization** | exposure→outcome × 5 estimators, β (95% CI) + inline CI, p, Q, Egger-intercept, F/instruments | STROBE-MR (sensitivity panel, instrument strength, pleiotropy) |
| **T5 Colocalization** | PP.H0–H4 per gene + verdict | coloc convention |
| **S1–S3 Supplement** | full per-gene enrichment; all MR pairs×methods; pre-registration power/estimability | STROBE-MR (all analyses), STREGA (volume-of-data) |

## Figures (Mermaid, per user)
- **F1 Analysis flow** — data sources → Aim 1/2/3 → adversarial verification → verdicts (a
  STROBE-style participant/analysis-flow diagram, adapted to datasets).
- **F2 Triangulation schematic** — enrichment (+), behavioural MR (null), cis-MR/coloc (null) →
  conclusion. Communicates the convergent-evidence logic STROBE-MR asks authors to make explicit.

## Reporting requirements we must satisfy in prose (checklist hooks)
- **Ancestry** stated (European) and portability limited — STREGA population-stratification item.
- **Rationale for the circadian gene set** pre-registered + hashed — STREGA gene/variant-choice item.
- **MR assumptions + F-stat/instrument count + sample overlap (ILAE has no UKB)** — STROBE-MR.
- **Multiple-testing** approach (matched-null empirical p; per-gene nominal) stated — STREGA.
- **coloc caveats** (single tissue, p-value ABF approximation, power) — coloc convention.
