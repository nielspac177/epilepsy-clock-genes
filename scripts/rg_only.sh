#!/bin/zsh
set -u
cd "/Users/nielspacheco/Desktop/Research/Rolston lab/Epilepsy_clock_Genes"
export PYTHONPATH=tools/ldsc/src
PY=/Users/nielspacheco/.epicirc-venv/bin/python
E="tools/ldsc_ref/eur_w_ld_chr/"
LOG=results/ldsc/RUN_rgonly.log; : > "$LOG"
rg(){ local out="results/ldsc/$1"; shift; local args=();
  for f in "$@"; do args+=(--rg "results/ldsc/${f}_m.sumstats.gz"); done
  echo "[$(date)] rg $out : $*" >>"$LOG"
  $PY -m ldsc.main ldsc rg "${args[@]}" --ref-ld-chr "$E" --w-ld-chr "$E" --out "$out" >>"$LOG" 2>&1; }
rg rg_gge2   gge_real chronotype sleepduration focal
rg rg_focal2 focal    chronotype sleepduration
echo "=== SUMMARIES ===" >>"$LOG"
for f in results/ldsc/rg_gge2*.log results/ldsc/rg_focal2*.log; do
  echo "--- $f ---" >>"$LOG"; grep -A6 "Summary of Genetic Correlation" "$f" >>"$LOG" 2>&1; done
echo "RGONLY_DONE" >>"$LOG"
