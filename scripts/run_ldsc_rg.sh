#!/bin/zsh
# LDSC genetic correlation (rg) of epilepsy (GGE, focal) x sleep traits (chronotype, sleep duration),
# plus partitioned h2 for focal & all-epilepsy (contrast to GGE). Local venv, belowlab LDSC.
set -u
cd "/Users/nielspacheco/Desktop/Research/Rolston lab/Epilepsy_clock_Genes"
export PYTHONPATH=tools/ldsc/src
PY=/Users/nielspacheco/.epicirc-venv/bin/python
REF=tools/ldsc_ref
LOG=results/ldsc/RUN_rg.log; : > "$LOG"
say(){ echo "[$(date)] $*" >>"$LOG"; }

# ---- 0. eur_w_ld_chr reference for rg ----
if [ ! -d "$REF/eur_w_ld_chr" ]; then
  say "downloading eur_w_ld_chr"
  curl -sSL --max-time 600 -o "$REF/eur_w_ld_chr.tar.gz" \
    "https://zenodo.org/api/records/8182036/files/eur_w_ld_chr.tar.gz/content" >>"$LOG" 2>&1
  tar -xzf "$REF/eur_w_ld_chr.tar.gz" -C "$REF" >>"$LOG" 2>&1
fi
say "eur_w_ld_chr files: $(ls $REF/eur_w_ld_chr/ 2>/dev/null | wc -l)"

munge(){  # munge <name> <extra args...>
  local name="$1"; shift
  say "munge $name"
  $PY -m ldsc.main munge_sumstats "$@" --out "results/ldsc/$name" >>"$LOG" 2>&1
  say "  -> $(ls -la results/ldsc/$name.sumstats.gz 2>/dev/null | awk '{print $5}') bytes"
}

# ---- 1. build pre-munge for focal + all epilepsy (SNP A1 A2 Z N P) from ILAE .tbl ----
for pair in "focal:ILAE3_Caucasian_focal_epilepsy_final.tbl" "allepi:ILAE3_Caucasian_all_epilepsy_final.tbl"; do
  nm="${pair%%:*}"; tbl="data/raw/final_sumstats/${pair##*:}"
  if [ ! -f "results/ldsc/$nm.sumstats.pre" ]; then
    say "pre-munge $nm from $tbl"
    awk 'NR==1{print "SNP\tA1\tA2\tZ\tN\tP"; next}{print $3"\t"$4"\t"$5"\t"$9"\t"$8"\t"$10}' \
      "$tbl" > "results/ldsc/$nm.sumstats.pre" 2>>"$LOG"
  fi
done

# ---- 2. munge everything ----
munge gge_real       --sumstats results/ldsc/gge.sumstats.pre   --signed-sumstats Z,0 --snp SNP --a1 A1 --a2 A2 --p P --N-col N
munge focal          --sumstats results/ldsc/focal.sumstats.pre --signed-sumstats Z,0 --snp SNP --a1 A1 --a2 A2 --p P --N-col N
munge allepi         --sumstats results/ldsc/allepi.sumstats.pre --signed-sumstats Z,0 --snp SNP --a1 A1 --a2 A2 --p P --N-col N
munge chronotype     --sumstats data/raw/sleep_full/chronotype_full.txt.gz \
                     --snp SNP --a1 ALLELE1 --a2 ALLELE0 --signed-sumstats LOGOR,0 --p P_BOLT_LMM --N 449734

if [ ! -f data/raw/sleep_full/sleepdurationsumstats.txt ]; then
  say "unzip sleep duration"; unzip -o -d data/raw/sleep_full data/raw/sleep_full/sleepduration_full.txt.zip >>"$LOG" 2>&1
fi
munge sleepduration  --sumstats data/raw/sleep_full/sleepdurationsumstats.txt \
                     --snp SNP --a1 ALLELE1 --a2 ALLELE0 --signed-sumstats BETA_SLEEPDURATION,0 --p P_SLEEPDURATION --N 446118

# ---- 3. rg matrix (unconstrained intercept; UKB not in ILAE -> low overlap) ----
E="$REF/eur_w_ld_chr/"
rg(){  # rg <outstem> <p1> <p2...>
  local out="results/ldsc/$1"; shift
  local args=(); for f in "$@"; do args+=(--rg "results/ldsc/$f.sumstats.gz"); done
  say "rg $out : $*"
  $PY -m ldsc.main ldsc rg "${args[@]}" --ref-ld-chr "$E" --w-ld-chr "$E" --out "$out" >>"$LOG" 2>&1
}
rg rg_gge   gge_real chronotype sleepduration focal   # rg(GGE, chrono/sleepdur/focal[ctrl])
rg rg_focal focal    chronotype sleepduration          # rg(focal, chrono/sleepdur)

# ---- 4. partitioned h2 for focal + all (contrast to GGE) ----
h2(){  # h2 <name>
  say "partitioned h2 $1"
  $PY -m ldsc.main ldsc h2 --h2 "results/ldsc/$1.sumstats.gz" \
    --ref-ld-chr $REF/baselineLD. \
    --w-ld-chr $REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --frqfile-chr $REF/1000G_Phase3_frq/1000G.EUR.QC. \
    --overlap-annot --print-coefficients --out "results/ldsc/${1}_partitioned" >>"$LOG" 2>&1
}
h2 focal
h2 allepi

say "collecting rg logs"
grep -A3 "Summary of Genetic Correlation" results/ldsc/rg_*.log >> "$LOG" 2>&1
say "RG_PIPELINE_DONE"
