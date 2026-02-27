#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

# Paths
counts_path  <- "data/external/simple_counts.txt"
samples_path <- "metadata/samples.tsv"
out_res      <- "results/tables/result.txt"
out_norm     <- "results/tables/norm-matrix-deseq2.txt"

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ---- Read sample sheet ----
samples <- read.csv(samples_path, stringsAsFactors = FALSE)
stopifnot(all(c("sample","condition") %in% colnames(samples)))
samples$condition <- factor(samples$condition)

# ---- Read counts (featureCounts-like) ----
# skip commented lines starting with '#'
raw <- read.table(counts_path,
                  header = TRUE,
                  sep = "",
                  quote = "",
                  comment.char = "#",
                  check.names = FALSE,
                  stringsAsFactors = FALSE)

# Expect first column = Geneid (or similar)
if (!("Geneid" %in% colnames(raw))) {
  stop("Expected a column named 'Geneid' in count file.")
}

# Keep only Geneid + sample columns
# In your file headers are like: SRR941816_out.bam, ...
wanted_cols <- c("Geneid", paste0(samples$sample, "_out.bam"))
missing <- setdiff(wanted_cols, colnames(raw))
if (length(missing) > 0) {
  message("Available columns:\n", paste(colnames(raw), collapse = ", "))
  stop("Missing expected columns: ", paste(missing, collapse = ", "))
}

counts_df <- raw[, wanted_cols, drop = FALSE]
rownames(counts_df) <- counts_df$Geneid
counts_df$Geneid <- NULL

# Rename columns to clean sample names
colnames(counts_df) <- samples$sample

# Ensure integer matrix
count_mat <- as.matrix(counts_df)
mode(count_mat) <- "numeric"
if (any(is.na(count_mat))) stop("NA values found in count matrix after parsing.")
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

# ---- DESeq2 ----
coldata <- data.frame(row.names = samples$sample,
                      condition = samples$condition)

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)

# Results: fermentation vs control
res <- results(dds, contrast = c("condition", "fermentation", "control"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Sort by adjusted p-value (FDR)
res_df <- res_df[order(res_df$padj, res_df$pvalue), ]

# Write outputs
write.table(res_df, file = out_res, sep = "\t", quote = FALSE, row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file = out_norm, sep = "\t", quote = FALSE)

message("Done.")
message("DE results: ", out_res)
message("Normalized counts: ", out_norm)