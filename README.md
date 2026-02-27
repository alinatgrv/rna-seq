# RNA-seq Differential Expression Analysis — Yeast Fermentation

Reproducible RNA-seq differential expression workflow comparing *Saccharomyces cerevisiae*
before and during fermentation.

## What is inside
- DESeq2 differential expression
- Heatmap of top 50 DE genes
- Optional volcano plot
- GO Slim interpretation (SGD)

## Structure
- `data/external/` — input count matrix (`simple_counts.txt`)
- `metadata/` — sample sheet (`samples.tsv`)
- `scripts/` — analysis scripts
- `results/` — output tables and figures
