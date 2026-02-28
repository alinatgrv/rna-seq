# Differential RNA expression analysis — Yeast fermentation (RNA-seq)

**Student:** Alina Tagirova

## 1. Goal

The goal of this project was to identify genes that change their expression levels when yeast cells transition from normal conditions to fermentation. Yeast are facultative anaerobes and can switch metabolism depending on oxygen availability: under aerobic conditions they rely on mitochondrial respiration, while under low oxygen they shift toward fermentation with ethanol and CO₂ production.

We used RNA-seq read count data from two conditions (control vs fermentation) with two replicates each, applied DESeq2 to find differentially expressed genes, and visualized the strongest expression changes using a heatmap.

## 2. Project structure and reproducibility

Structure (separating input data, scripts, results, and report):

- data/external/ — input count matrix (simple_counts.txt)
- metadata/ — sample sheet (samples.tsv)
- scripts/ — reproducible analysis scripts
- results/tables/ — output tables (DE results and normalized counts)
- results/figures/ — plots
- environment/ — reproducible environment specification
- report/ — written report (report.md)

### 2.1 Conda/Mamba environment

A dedicated environment was created to ensure that DESeq2 and plotting dependencies are reproducible.

Environment activation:

```bash
mamba activate rna_seq
```

Environment export:
```bash
mamba env export -n rna_seq --no-builds > environment/rna_seq.yml
```


---

## 3. Input data

### 3.1 Count matrix

We used a gene-by-sample count matrix simple_counts.txt containing raw RNA-seq read counts for four samples:

- SRR941816 — control replicate 1
- SRR941817 — control replicate 2
- SRR941818 — fermentation replicate 1
- SRR941819 — fermentation replicate 2

Each row corresponds to a gene, and each column corresponds to the number of reads mapped to that gene in the corresponding sample.

The file was placed into:

- data/external/simple_counts.txt

### 3.2 Sample sheet (experimental design)

A sample sheet mapping sample IDs to conditions was created.

Path:

- metadata/samples.tsv

Content:

```text
sample,condition
SRR941816,control
SRR941817,control
SRR941818,fermentation
SRR941819,fermentation
```

This encoding defines the first two samples as controls and the last two samples as fermentation samples. The condition factor is later used by DESeq2 to estimate expression differences between fermentation and control.


---

## 4. Differential expression analysis (DESeq2)

Differential expression analysis was performed in R using DESeq2, which models raw counts using a negative binomial distribution and performs statistical testing with multiple testing correction.

### 4.1 Script and command

The analysis was implemented as a reproducible script:

- scripts/01_deseq2.R

Run command:

```bash
Rscript --vanilla scripts/01_deseq2.R
```

### 4.2 What DESeq2 does in this step

DESeq2 performs:

- library size normalization (size factors)
- dispersion estimation
- model fitting for each gene
- statistical testing for differential expression
- multiple testing correction (adjusted p-values, padj)

The comparison is defined as:

- fermentation vs control

Interpretation of log2FoldChang:

- log2FoldChange > 0 → higher expression in fermentation
- log2FoldChange < 0 → higher expression in control

### 4.3 Outputs produced

The script writes two key outputs into results/tables/.

1) Differential expression results  
Path: results/tables/result.txt

This table includes (for each gene):

- baseMean
- log2FoldChange
- standard error (lfcSE)
- test statistic (stat)
- p-value (pvalue)
- adjusted p-value (padj)
- gene_id

2) Normalized counts matrix  
Path: results/tables/norm-matrix-deseq2.txt

This file contains DESeq2-normalized counts used for visualization.

During execution, DESeq2 reported the expected internal steps such as:

- estimating size factors

- estimating dispersions

- fitting model and testing

### 4.4 Question

**Question 1**: Based on the colData, which samples are controls and which are fermentation samples? How can you tell?

Answer: Controls are SRR941816 and SRR941817 because they are assigned condition = control in metadata/samples.tsv. Fermentation samples are SRR941818 and SRR941819 because they are assigned condition = fermentation. DESeq2 uses this design table to group replicates and estimate the fermentation vs control effect.

**Question 2**: For the top 6 genes (smallest padj), which are upregulated in fermentation and which are downregulated?

How it was checked: the top of result.txt is already sorted by padj, and the sign of log2FoldChange indicates direction.

Example check:

```bash
head -n 7 results/tables/result.txt
```

Interpretation rule:

- positive log2FoldChang → upregulated in fermentation

- negative log2FoldChange → downregulated in fermentation (i.e., higher in control)

At this stage, the most significant genes at the top showed strongly positive log2FoldChange, indicating upregulation in fermentation.


---

## 5. Heatmap visualization (top 50 DE genes)

A heatmap was generated to visualize expression patterns of the most significant DE genes and to check that replicates cluster by condition.
### 5.1 Script and command

Script:

- scripts/02_heatmap.R

Run:

```bash
Rscript --vanilla scripts/02_heatmap.R
```


---

### 5.2 Method

Input:

- normalized counts (results/tables/norm-matrix-deseq2.txt)
- top 50 genes selected from the DESeq2 results (results/tables/result.txt)

For visualization, counts were transformed into gene-wise Z-scores (row-wise standardization) so that each gene’s relative up/down pattern across samples is visible regardless of the absolute level of expression.

Output plot:

- results/figures/heatmap_top50.pdf
- results/figures/heatmap_top50.png

### Heatmap of top 50 DE genes

![Heatmap of top 50 DE genes](../results/figures/heatmap_top50.png)

### 5.3 Heatmap interpretation

The two control replicates (SRR941816, SRR941817) cluster together, and the two fermentation replicates (SRR941818, SRR941819) cluster together, indicating good reproducibility between replicates.

The heatmap shows clear blocks of genes with low expression in controls and high expression in fermentation samples (and a smaller block with the opposite pattern), confirming strong transcriptional differences between conditions.

## 6. Gene Ontology (GO) analysis (SGD GO Slim Mapper)

To interpret biological functions of the most significant DE genes, Gene Ontology Slim analysis was performed using the Saccharomyces Genome Database (SGD) GO Slim Mapper (Biological Process).

### 6.1 Preparing gene lists

The DESeq2 results file (results/tables/result.txt) is sorted by adjusted p-values (padj). We extracted the top 50 most significant genes and prepared gene lists for GO analysis.

Top 50 gene list:

```bash
scripts/03_top50_genes.sh
# output: results/tables/top50_genes.txt
```

Because the top 50 genes were mostly upregulated, we also split top 50 genes by direction (sign of log2FoldChange):

```bash
scripts/04_split_top50_up_down.sh
# outputs:
# results/tables/top50_up_genes.txt
# results/tables/top50_down_genes.txt
```

To obtain a robust list of downregulated genes, we additionally selected the top 50 downregulated genes across the whole results table (lfc < 0, sorted by padj):

```bash
scripts/05_top50_down_by_padj.sh
# output: results/tables/top50_down_by_padj.txt
```


---

### 6.2 GO Slim Mapper settings

- Tool: Gene Ontology Slim Term Mapper
- GO set: Yeast GO-Slim: process
- Terms: SELECT ALL TERMS

Inputs:

- results/tables/top50_up_clean.txt (upregulated genes)
- results/tables/top50_down_clean.txt (downregulated genes)

For upregulated genes (48 annotated genes):

- rRNA processing: 11/48 genes (22.92%) vs 5.18% genome background
- ribosomal large subunit biogenesis: 16.67% vs 2.01%
- ribosome assembly: 12.5% vs 1.25%

These categories are strongly overrepresented compared to genome background, indicating activation of ribosome biogenesis and RNA processing during fermentation.

For downregulated genes (50 genes):

- carbohydrate metabolic process: 18% vs 2.35%
- lipid metabolic process: 16% vs 4.61%
- cellular respiration: 10% vs 1.73%

These terms are also enriched relative to the genome, suggesting suppression of respiratory and oxidative metabolism.

Full GO Slim output tables are available in:

- results/go/2d1746615348f030f856e71f23c6a33d.html
- results/go/2fc3548b8b37f3ccc01300f278560809.html

**Biological interpretation:**

The transcriptomic profile indicates a clear metabolic reprogramming during the transition to fermentation.

Upregulation of ribosome biogenesis, rRNA processing and RNA polymerase I activity suggests activation of translational capacity and ribosomal production. This may reflect the need for rapid protein synthesis during adaptation to a new metabolic state.

Downregulation of carbohydrate metabolism, lipid metabolism and cellular respiration pathways is consistent with reduced oxidative metabolism and a shift away from mitochondrial respiration toward fermentative ATP production.

Together, these results support the expected biological model: yeast cells suppress aerobic respiration and reprogram gene expression to support fermentation.