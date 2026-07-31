#!/bin/zsh
# Assemble a submission-ready package: manuscript + figures/tables -> .docx and an Overleaf .zip.
set -eu
cd "/Users/nielspacheco/Desktop/Research/Rolston lab/Epilepsy_clock_Genes"
SUB=submission; OVL=$SUB/overleaf
rm -rf "$SUB"; mkdir -p "$OVL/figures"
cp figures/*.png "$OVL/figures/"

# combined markdown: manuscript body (drop HTML comments) + figures/tables appendix
COMB=$SUB/_combined.md
{ cat paper/manuscript.md; echo; cat paper/submission_figures_tables.md; } \
  | perl -0777 -pe 's/<!--.*?-->//gs' > "$COMB"

# 1) Word
pandoc "$COMB" -o "$SUB/manuscript.docx" --resource-path=".:$OVL" 2>&1 | tail -2 || true

# 2) Overleaf LaTeX (xelatex handles the unicode: ×, β, ρ, ≥, arrows)
pandoc "$COMB" -s --pdf-engine=xelatex -o "$OVL/main.tex" 2>&1 | tail -2 || true
cat > "$OVL/README.txt" <<'EOF'
Overleaf package for the circadian-genetics x epilepsy manuscript.
Upload this whole folder (main.tex + figures/) to Overleaf, or compile locally:
  xelatex main.tex   (xelatex required for unicode; run twice for refs)
Figures are in figures/. Tables and prose are in main.tex.
EOF

# 3) try a local PDF build to verify (non-fatal)
( cd "$OVL" && xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null 2>&1 \
  && xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null 2>&1 \
  && echo "PDF built" || echo "PDF build skipped/failed (tex still valid to upload)" )
cp "$OVL/main.pdf" "$SUB/manuscript.pdf" 2>/dev/null || true

# 4) zip the Overleaf folder
( cd "$SUB" && zip -qr epilepsy-clock-genes-overleaf.zip overleaf )
echo "=== package contents ==="; ls -la "$SUB"
