#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

message("Starting volcano plot...")

res_path <- "results/tables/result.txt"
out_png  <- "results/figures/volcano.png"

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load DESeq2 results
res <- read.table(res_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Basic checks
required <- c("log2FoldChange", "pvalue", "padj", "gene_id")
missing <- setdiff(required, colnames(res))
if (length(missing) > 0) {
  stop("Missing columns in result.txt: ", paste(missing, collapse = ", "))
}

# Prepare data
res$padj_num <- suppressWarnings(as.numeric(res$padj))
res$pvalue_num <- suppressWarnings(as.numeric(res$pvalue))
res$log2FC <- suppressWarnings(as.numeric(res$log2FoldChange))

# Use padj if available, otherwise fallback to pvalue for plotting
res$neglog10 <- -log10(ifelse(!is.na(res$padj_num), res$padj_num, res$pvalue_num))

# Mark significance
padj_thr <- 0.05
lfc_thr  <- 1

res$status <- "not_significant"
res$status[!is.na(res$padj_num) & res$padj_num < padj_thr & res$log2FC >=  lfc_thr] <- "up"
res$status[!is.na(res$padj_num) & res$padj_num < padj_thr & res$log2FC <= -lfc_thr] <- "down"

# Select genes to label: top by padj among significant; fallback to pvalue if needed
lab_n <- 10
sig <- res[!is.na(res$padj_num) & res$padj_num < padj_thr, ]
if (nrow(sig) >= 1) {
  sig <- sig[order(sig$padj_num), ]
} else {
  sig <- res[!is.na(res$pvalue_num), ]
  sig <- sig[order(sig$pvalue_num), ]
}
label_df <- head(sig, lab_n)
# Plot
p <- ggplot(res, aes(x = log2FC, y = neglog10)) +
  geom_point(aes(alpha = status), size = 1.2) +
  scale_alpha_manual(values = c(
    up = 0.9,
    down = 0.9,
    not_significant = 0.35
  ), guide = "none") +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
  geom_text_repel(
    data = label_df,
    aes(label = gene_id),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano plot (fermentation vs control)",
    x = "log2FoldChange",
    y = expression(-log[10]("padj"))
  ) +
  theme_minimal()

ggsave(out_png, plot = p, width = 10, height = 7, dpi = 200)

message("Done.")
message("Volcano: ", out_png)
