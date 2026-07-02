# ADR-0001: Snakemake as the reproducibility engine

**Status:** accepted · **Date:** 2026-07-02

## Context
Every result must regenerate deterministically from pinned public summary statistics, with a
provenance trail acceptable to a genetics journal and a reviewer re-running the pipeline.

## Decision
Use **Snakemake** for the analysis DAG (download → munge/harmonize → Aim 1/2/3 → figures →
report). External binaries (LDSC, MAGMA, PLINK, GCTA, Rscript) are pinned via `environment.yml`
and a container; large inputs are referenced by checksum in `config/traits_manifest.yaml`, never
committed.

## Alternatives
- **Nextflow** — better cloud/container scaling but heavier; overkill at summary-stat scale.
- **Plain Make/Python** — weaker provenance and re-run guarantees.

## Consequences
`snakemake -n` reproduces the DAG from a clean checkout before any large download; CI can dry-run
the graph nightly. Python-native rules keep the wrappers (`src/epicirc/*`) directly importable and
unit-tested.
