#!/bin/zsh
# Re-munge with --merge-alleles (w_hm3.snplist) to fix the "incompatible alleles" rg error, then
# run the rg matrix. Waits for the first LDSC pipeline (RUN_rg.log) to finish to avoid file races.
set -u
cd "/Users/nielspacheco/Desktop/Research/Rolston lab/Epilepsy_clock_Genes"
export PYTHONPATH=tools/ldsc/src
PY=/Users/nielspacheco/.epicirc-venv/bin/python
REF=tools/ldsc_ref
HM3=$REF/eur_w_ld_chr/w_hm3.snplist
E="$REF/eur_w_ld_chr/"
LOG=results/ldsc/RUN_rg2.log; : > "$LOG"
say(){ echo "[$(date)] $*" >>"$LOG"; }

# wait (up to ~30 min) for the first pipeline to release the munged files
for i in $(seq 1 300); do
  grep -q RG_PIPELINE_DONE results/ldsc/RUN_rg.log 2>/dev/null && break
  sleep 6
done
say "first pipeline done; re-munging with --merge-alleles $HM3"

mungem(){ local name="$1"; shift; say "munge(m) $name";
  $PY -m ldsc.main munge_sumstats "$@" --merge-alleles "$HM3" --out "results/ldsc/${name}_m" >>"$LOG" 2>&1;
  say "  -> $(ls -la results/ldsc/${name}_m.sumstats.gz 2>/dev/null | awk '{print $5}') bytes"; }

mungem gge_real     --sumstats results/ldsc/gge.sumstats.pre   --signed-sumstats Z,0 --snp SNP --a1 A1 --a2 A2 --p P --N-col N
mungem focal        --sumstats results/ldsc/focal.sumstats.pre --signed-sumstats Z,0 --snp SNP --a1 A1 --a2 A2 --p P --N-col N
mungem chronotype   --sumstats data/raw/sleep_full/chronotype_full.txt.gz \
                    --snp SNP --a1 ALLELE1 --a2 ALLELE0 --signed-sumstats LOGOR,0 --p P_BOLT_LMM --N 449734
mungem sleepduration --sumstats data/raw/sleep_full/sleepdurationsumstats.txt \
                    --snp SNP --a1 ALLELE1 --a2 ALLELE0 --signed-sumstats BETA_SLEEPDURATION,0 --p P_SLEEPDURATION --N 446118

rg(){ local out="results/ldsc/$1"; shift; local args=();
  for f in "$@"; do args+=(--rg "results/ldsc/${f}_m.sumstats.gz"); done
  say "rg $out : $*"
  $PY -m ldsc.main ldsc rg "${args[@]}" --ref-ld-chr "$E" --w-ld-chr "$E" --out "$out" >>"$LOG" 2>&1; }

rg rg_gge2   gge_real chronotype sleepduration focal
rg rg_focal2 focal    chronotype sleepduration

say "=== RG SUMMARIES ==="
for f in results/ldsc/rg_gge2*.log results/ldsc/rg_focal2*.log; do
  echo "--- $f ---" >>"$LOG"; grep -A6 "Summary of Genetic Correlation" "$f" >>"$LOG" 2>&1
done
say "RG2_DONE"
