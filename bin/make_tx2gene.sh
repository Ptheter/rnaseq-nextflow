#!/usr/bin/env bash
# Builds a transcript -> gene ID -> gene symbol table from the GENCODE transcript FASTA headers,
# so it always matches the Salmon index built from the same file.
set -euo pipefail
fa="${1:-reference/transcripts.fa.gz}"
out="${2:-reference/tx2gene.tsv}"
zcat "$fa" | grep '^>' | tr -d '>' \
  | awk -F'|' 'BEGIN{OFS="\t"} {print $1,$2,$6}' | sort -u > "$out"
wc -l "$out"
