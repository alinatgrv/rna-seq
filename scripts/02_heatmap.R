#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pheatmap)
})

message("Starting heatmap pipeline...")

res_path  <- "results/tables/result.txt"
norm_path <- "results/tables/norm-matrix-deseq2.txt"
out_pdf   <- "results/figures/heatmap_top50.pdf"

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Read DE results
res <- read.table(res_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Keep only rows with padj
res <- res[!is.na(res$padj), ]
res <- res[order(res$padj, res$pvalue), ]

top <- head(res$gene_id, 50)
if (length(top) < 2) stop("Too few genes for heatmap. Check result.txt.")

# Read normalized matrix
norm <- read.table(norm_path, header = TRUE, sep = "\t", check.names = FALSE)

# If first column is gene id, turn it into rownames.
# (In our output from DESeq2 counts(), rownames are stored as first column when writing table.)
if (!any(rownames(norm) %in% top)) {
  # Try treating the first column as gene IDs
  rownames(norm) <- norm[, 1]
  norm <- norm[, -1, drop = FALSE]
}

# Subset top genes
mat <- norm[top, , drop = FALSE]

# Z-score per gene (row-wise): center/scale across samples
zmat <- t(scale(t(as.matrix(mat))))

# Replace any NaN (happens if a gene has zero variance across samples)
zmat[is.nan(zmat)] <- 0

# Plot
pdf(out_pdf, width = 7, height = 10)
pheatmap(
  zmat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "Top 50 DE genes (z-score, DESeq2 normalized counts)"
)
dev.off()

message("Done.")
message("Heatmap: ", out_pdf)
