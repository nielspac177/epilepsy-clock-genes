# Data access — how to obtain the real GWAS inputs

The pipeline runs in **synthetic mode** out of the box. To run on real data, download each dataset
and drop it in `data/raw/` under the exact filename the manifest expects, then re-run Snakemake.
Nothing here needs a credit card; a few need a free registration or a data-use agreement (DUA).

---

## 1. epiGAD — ILAE Consortium 2023 epilepsy GWAS (the core outcome data)

**What:** European-ancestry summary statistics for all-epilepsy, focal, generalized (GGE), and the
7 subtypes (JME, CAE, JAE, GTCS-alone, focal-HS, focal-lesion-negative, focal-other-lesion).
Paper: *Nature Genetics* 2023, DOI 10.1038/s41588-023-01485-w.

**Steps:**
1. Go to **https://www.epigad.org/** (EPIGAD v5.0 — the ILAE Genetics Commission portal).
2. Find the **2023 ILAE Consortium on Complex Epilepsies** GWAS release. Choose the
   **European-only** meta-analysis files (our pipeline is EUR-primary).
3. epiGAD typically serves these as a direct download or after a short **registration form**
   (name + institutional email). If prompted, use your UCSF/Rolston-lab email. There is usually no
   long DUA for the ILAE common-variant summary stats, but read the terms shown.
4. Download the per-phenotype files. You'll likely get one gzipped table per phenotype (focal,
   GGE, subtypes) plus the all-epilepsy file.

> If you'd rather not register, tell me and I'll check whether the release is mirrored on the
> GWAS Catalog (https://www.ebi.ac.uk/gwas/) under the ILAE study accession, which is login-free.

**Where to put them:** rename to match `config/traits_manifest.yaml` keys, e.g.
`data/raw/focal_epilepsy.sumstats.gz`, `data/raw/generalized_gge.sumstats.gz`, etc. (I'll add the
column-mapping step so whatever header epiGAD uses is normalized to our schema — send me one file's
header and I'll wire it exactly.)

---

## 2. Circadian / sleep exposures (for Aim 2 MR + Aim 1 rg) — login-free

- **Morning chronotype** (Jones 2019) and **insomnia** (Jansen 2019), **sleep duration** (Dashti
  2019), **daytime napping** (Dashti 2021): all on the **Sleep Disorder Knowledge Portal**
  (https://sleep.hugeamp.org/) and the **GWAS Catalog**. Download the EUR summary stats.
- Save as `data/raw/chronotype_morning.sumstats.gz`, `insomnia.sumstats.gz`, etc.

## 3. Aim 3 + confounders

- **Drug-resistant epilepsy** (eBioMedicine 2025): summary stats via the paper's data-availability
  section (GWAS Catalog accession). → `data/raw/drug_resistant_epilepsy.sumstats.gz`.
- **Psychiatric (MDD/BIP/ADHD)** and **metabolic (BMI/T2D)**: PGC (https://pgc.unc.edu/) and GIANT/
  DIAGRAM. Note: **PGC MDD full release contains 23andMe and needs a DUA**; the 23andMe-excluded
  version is open. BIP/ADHD are largely open.

## 4. cis instruments (to claim "clock biology", ADR-0004) — login-free

- **GTEx v8** (brain tissues) and **eQTLGen** cis-eQTLs; **deCODE/UKB-PPP** pQTLs where available.

## 5. Reference panels + tools (I attempted these; status below)

| Asset | Status | Note |
|---|---|---|
| PLINK 1.9 (macOS) | ✅ in `tools/` | runs (x86_64 via Rosetta) |
| GCTA (macOS) | ✅ in `tools/` | x86_64 |
| MAGMA (macOS) | ❌ | fetched a Linux build by mistake; get the macOS binary from https://cncr.nl/research/magma/ |
| MAGMA aux (NCBI37, g1000_eur) | ❌ | URLs moved; on the MAGMA site's SURFdrive |
| LDSC eur LD scores + w_hm3.snplist | ❌ | Zenodo path changed; get from the LDSC repo's current release |
| LDSC itself | ⚠️ blocker | Python-2 tool; needs its own `conda create -n ldsc python=2.7` env (system conda's solver is currently broken — `libarchive` missing) |

---

## After you drop files in `data/raw/`
```bash
cd "<repo>"
.venv/bin/snakemake -s workflow/Snakefile --cores 4   # will use real files where present
```
Send me one epiGAD file header and I'll add the exact column-normalization so the real run is
turnkey.
