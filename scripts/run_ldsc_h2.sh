#!/bin/zsh
# Real-data LDSC univariate + partitioned heritability for GGE.
# Self-contained: needs only results/ldsc/gge.sumstats.pre + tools/ldsc_ref/*.
# Runs as a background job (no 2-min foreground cap). Writes a run log to results/ldsc/RUN.log.
set -u
cd "/Users/nielspacheco/Desktop/Research/Rolston lab/Epilepsy_clock_Genes"
export PYTHONPATH=tools/ldsc/src
PY=/Users/nielspacheco/.epicirc-venv/bin/python   # local (non-iCloud) venv; instant, stable
LOG=results/ldsc/RUN.log
: > "$LOG"
echo "[$(date)] import smoke-test..." >>"$LOG"
$PY -c "import ldsc.ldscore.sumstats; print('import OK')" >>"$LOG" 2>&1 || { echo "IMPORT FAILED" >>"$LOG"; exit 3; }

echo "[$(date)] === STEP 1: munge GGE ===" >>"$LOG"
$PY -m ldsc.main munge_sumstats \
  --sumstats results/ldsc/gge.sumstats.pre \
  --signed-sumstats Z,0 --snp SNP --a1 A1 --a2 A2 --p P --N-col N \
  --out results/ldsc/gge_real >>"$LOG" 2>&1
echo "[$(date)] munge exit=$?  outputs:" >>"$LOG"
ls -la results/ldsc/gge_real* >>"$LOG" 2>&1

MUNGED=results/ldsc/gge_real.sumstats.gz
[ -f "$MUNGED" ] || MUNGED=$(ls results/ldsc/gge_real*.sumstats.gz 2>/dev/null | head -1)
echo "[$(date)] using munged=$MUNGED" >>"$LOG"

echo "[$(date)] === STEP 2: partitioned h2 (baselineLD v2.2) ===" >>"$LOG"
$PY -m ldsc.main ldsc h2 \
  --h2 "$MUNGED" \
  --ref-ld-chr tools/ldsc_ref/baselineLD. \
  --w-ld-chr tools/ldsc_ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --frqfile-chr tools/ldsc_ref/1000G_Phase3_frq/1000G.EUR.QC. \
  --overlap-annot --print-coefficients \
  --out results/ldsc/gge_partitioned >>"$LOG" 2>&1
echo "[$(date)] h2 exit=$?  outputs:" >>"$LOG"
ls -la results/ldsc/gge_partitioned* >>"$LOG" 2>&1
echo "[$(date)] DONE" >>"$LOG"
