#!/usr/bin/env bash
set -euo pipefail

RES="results/tables/result.txt"
OUT="results/tables/top50_genes.txt"

mkdir -p results/tables

# result.txt is tab-separated with header, gene_id in last column
# take top 50 rows after header, print gene_id column
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="gene_id") col=i; next}
            NR<=51 {print $col}' "$RES" > "$OUT"

echo "Wrote: $OUT"
wc -l "$OUT"
