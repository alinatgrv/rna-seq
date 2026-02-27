#!/usr/bin/env bash
set -euo pipefail

RES="results/tables/result.txt"
OUT_UP="results/tables/top50_up_genes.txt"
OUT_DN="results/tables/top50_down_genes.txt"

mkdir -p results/tables

# Find column numbers for gene_id and log2FoldChange from header
gene_col=$(head -n 1 "$RES" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="gene_id") print i}')
lfc_col=$(head -n 1 "$RES" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="log2FoldChange") print i}')

if [[ -z "${gene_col}" || -z "${lfc_col}" ]]; then
  echo "Could not detect gene_id or log2FoldChange columns in $RES" >&2
  exit 1
fi

# Take top 50 rows after header (already sorted by padj), split by sign of log2FC
tail -n +2 "$RES" | head -n 50 | \
awk -F'\t' -v g="$gene_col" -v l="$lfc_col" '
  $l > 0 {print $g > "'"$OUT_UP"'"}
  $l < 0 {print $g > "'"$OUT_DN"'"}
'

# Report counts
echo "Wrote: $OUT_UP ($(wc -l < "$OUT_UP" 2>/dev/null || echo 0) genes)"
echo "Wrote: $OUT_DN ($(wc -l < "$OUT_DN" 2>/dev/null || echo 0) genes)"
