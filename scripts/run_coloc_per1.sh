#!/bin/zsh
# coloc(GGE, chronotype) at the PER1 locus (17p13.1). Pure-python (base .venv), no neuro venv.
# PER1 hg19 chr17:8,043,790-8,059,824; window ±~200kb.
set -u
cd "/Users/nielspacheco/Desktop/Research/Rolston lab/Epilepsy_clock_Genes"
OUT=results/coloc_multi; mkdir -p "$OUT"
LOG=$OUT/RUN_per1.log; : > "$LOG"
CHR=17; LO=7850000; HI=8250000

echo "[$(date)] slicing GGE chr$CHR:$LO-$HI" >>"$LOG"
awk -v c=$CHR -v lo=$LO -v hi=$HI 'NR==1{print;next} ($1==c && $2>=lo && $2<=hi){print}' \
  data/raw/final_sumstats/ILAE3_Caucasian_GGE_final.tbl > "$OUT/gge_per1.tbl" 2>>"$LOG"
echo "[$(date)] GGE region lines: $(wc -l < $OUT/gge_per1.tbl)" >>"$LOG"

echo "[$(date)] slicing chronotype chr$CHR:$LO-$HI" >>"$LOG"
gzcat data/raw/sleep_full/chronotype_full.txt.gz | \
  awk -v c=$CHR -v lo=$LO -v hi=$HI 'NR==1{print;next} ($2==c && $3>=lo && $3<=hi){print}' \
  > "$OUT/chronotype_per1.txt" 2>>"$LOG"
echo "[$(date)] chronotype region lines: $(wc -l < $OUT/chronotype_per1.txt)" >>"$LOG"

echo "[$(date)] running coloc_pairwise (GGE cc x chronotype quant)" >>"$LOG"
PYTHONPATH=src .venv/bin/python -m epicirc.mr.coloc_multi \
  --mode pair \
  --trait "GGE=$OUT/gge_per1.tbl=0,1,11,12" \
  --trait "chronotype=$OUT/chronotype_per1.txt=1,2,7,8" \
  --region "17:$LO-$HI" --wtype cc --wtype quant \
  --out "$OUT/coloc_gge_chronotype_per1.tsv" >>"$LOG" 2>&1
echo "[$(date)] result:" >>"$LOG"; cat "$OUT/coloc_gge_chronotype_per1.tsv" >>"$LOG" 2>&1
echo "COLOC_PER1_DONE" >>"$LOG"
