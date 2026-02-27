#!/usr/bin/env bash
set -euo pipefail

RES="results/tables/result.txt"
OUT="results/tables/top50_down_by_padj.txt"

mkdir -p results/tables

gene_col=$(head -n 1 "$RES" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="gene_id") print i}')
lfc_col=$(head -n 1 "$RES" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="log2FoldChange") print i}')
padj_col=$(head -n 1 "$RES" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="padj") print i}')

if [[ -z "${gene_col}" || -z "${lfc_col}" || -z "${padj_col}" ]]; then
  echo "Could not detect required columns in $RES" >&2
  exit 1
fi

# Filter lfc < 0, sort by padj (numeric), take top 50, output gene_id
tail -n +2 "$RES" | \
awk -F'\t' -v g="$gene_col" -v l="$lfc_col" -v p="$padj_col" '
  $l < 0 && $p != "NA" {print $p "\t" $g}
' | sort -g -k1,1 | head -n 50 | cut -f2 > "$OUT"

echo "Wrote: $OUT"
wc -l "$OUT"
