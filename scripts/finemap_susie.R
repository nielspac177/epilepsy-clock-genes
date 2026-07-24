#!/usr/bin/env Rscript
# SuSiE-RSS fine-mapping of the GGE 17p13.1 (PER1) locus, using harmonized z + 1000G EUR LD.
# Outputs per-SNP PIP and 95% credible sets, annotated with position so the causal gene can be read.
suppressMessages(library(susieR))

z   <- scan("results/finemap/z_aligned.txt", quiet = TRUE)
R   <- as.matrix(read.table("results/finemap/ld_aligned.txt"))
snp <- read.table("results/finemap/snps_aligned.tsv", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)
n   <- round(median(snp$neff))
cat(sprintf("SNPs=%d  n=%d  lead|z|=%.2f\n", length(z), n, max(abs(z))))

fit <- susie_rss(z = z, R = R, n = n, L = 10,
                 estimate_residual_variance = FALSE, max_iter = 500)

snp$pip <- fit$pip
cs <- fit$sets$cs
snp$cs <- NA_integer_
if (!is.null(cs)) for (k in seq_along(cs)) snp$cs[cs[[k]]] <- as.integer(k)

# genes overlapping the credible-set window (for annotation, using hg19 coords)
write.table(snp[order(-snp$pip), c("snp","bp","z","pip","cs")],
            "results/finemap/susie_pip.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== credible sets ===\n")
if (is.null(cs) || length(cs) == 0) {
  cat("no credible set (no confidently fine-mapped signal under L=10)\n")
} else {
  for (k in seq_along(cs)) {
    idx <- cs[[k]]
    ss  <- snp[idx, ]
    ss  <- ss[order(-ss$pip), ]
    cat(sprintf("CS%d: %d SNPs, span %d-%d, purity=%.2f, top: %s (bp %d, PIP %.3f)\n",
                k, length(idx), min(ss$bp), max(ss$bp), fit$sets$purity$min.abs.corr[k],
                ss$snp[1], ss$bp[1], ss$pip[1]))
  }
}
cat("\ntop 10 SNPs by PIP:\n")
top <- snp[order(-snp$pip), ][1:10, c("snp","bp","z","pip","cs")]
print(top, row.names = FALSE)
