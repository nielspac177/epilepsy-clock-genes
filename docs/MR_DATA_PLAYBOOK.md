# MR data-fetch playbook (adapted for circadian → epilepsy)

Reusable recipe for fetching public GWAS and running two-sample MR with **Python + curl only**
(no PLINK / genomics libraries). Adapted from a prior intracranial-aneurysm project. Implemented
in `epicirc.mr.clump` + `epicirc.mr.real_mr`.

## The setup here
- **Outcome = ILAE3 epilepsy** (already local: `data/raw/final_sumstats/ILAE3_*_final.tbl`,
  GRCh37, effect allele = Allele1). Use the **Caucasian_focal** and **Caucasian_GGE** files.
- **Exposure = circadian/sleep traits** — fetch from public repositories (below).
- **Sample overlap is low & conservative**: ILAE has **no UK Biobank**; UKB-based sleep exposures
  overlap the ILAE controls only marginally, biasing toward the null. Still report the bivariate
  LDSC intercept once LDSC is installed (ADR-0004).

## ⚠️ Build safety (the #1 gotcha)
ILAE is **GRCh37/hg19**. The exposure MUST be GRCh37, matched on **CHR:POS, not rsID** (many sleep
files use `chr:pos` IDs). A GWAS-Catalog **harmonised `*.h.tsv.gz` is often GRCh38** → it will match
~0 SNPs. Prefer the author's native GRCh37 file, or a `build37` release, or liftOver.

## Fetching exposures (login-free)
- **GWAS Catalog EBI FTP by accession** — resolve the study, then pull its summary-stats file:
  ```bash
  curl -s "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId=<PMID>"
  # summary-stats live under:
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST<range>/GCST<id>/
  ```
  Sleep traits: chronotype/morningness (Jones 2019), insomnia (Jansen 2019), sleep duration
  (Dashti 2019), daytime napping (Dashti 2021). Resolve each accession from its PMID via the API
  (do not hard-code — accessions drift).
- **Figshare API** (some datasets): `curl -s https://api.figshare.com/v2/articles/<id>` → parse
  `files[].download_url`.
- **Direct-download page scrape**: `curl -sL <page> | grep -ioE 'href="[^"]*(gz|zip)"'`.

## Pre-filter big files before loading
Exposure files are 10M+ rows. Filter to instruments first:
```bash
zcat exposure.tsv.gz | awk -F'\t' 'NR==1 || $<pcol> < 5e-8' > exposure.gwsig.tsv
```

## Run
```bash
PYTHONPATH=src python -m epicirc.mr.real_mr \
  --exposure exposure.gwsig.tsv \
  --outcome  data/raw/final_sumstats/ILAE3_Caucasian_GGE_final.tbl \
  --exp-chr CHR --exp-pos POS --exp-ea Effect_Allele --exp-oa Other_Allele \
  --exp-beta Effect --exp-p Pval \
  --clump-kb 1000 --out results/mr_real/chronotype__gge.tsv
```
Distance-clumping (±1 Mb) approximates LD pruning (no panel) — state this as a limitation; swap in
PLINK clumping once the LD panel is installed.

## Positive control (do this BEFORE trusting any null)
Run a known causal exposure→outcome pair through the identical driver and confirm it recovers the
expected effect. If it doesn't, the build/pipeline is broken, not the biology. (For MR method
validation generally; e.g. BMI→T2D.)

## Escalation to defensible causal claims (adversarial review, ADR-0004)
Behavioural-trait MR alone licenses "sleep behaviour is causal", NOT "clock biology". To claim
circadian biology, add cis-eQTL/pQTL MR on core clock genes + colocalization, MVMR on BMI/T2D +
psychiatric, overlap-aware MR (MRlap/CAUSE), and index-event-bias correction for the type outcome.
