#!/usr/bin/env bash
# Builds the Salmon reference: GENCODE transcripts restricted to those present in the chromosome-only
# annotation (drops alt/patch-sequence copies), plus a matching transcript -> gene -> symbol table.
# usage: prepare_reference.sh gencode.vXX.transcripts.fa.gz gencode.vXX.annotation.gtf.gz [outdir]
set -euo pipefail
fa="${1:?transcripts fasta (.gz) required}"; gtf="${2:?annotation gtf (.gz) required}"; out="${3:-reference}"
mkdir -p "$out"
tmp=$(mktemp)
zcat "$gtf" | awk -F'\t' '$3=="transcript"' | sed -E 's/.*transcript_id "([^"]+)".*/\1/' | sort -u > "$tmp"
zcat "$fa" | awk -v list="$tmp" 'BEGIN{while((getline l<list)>0) keep[l]=1} /^>/{split(substr($0,2),a,"|"); p=(a[1] in keep)} p' | gzip > "$out/transcripts_chr.fa.gz"
zcat "$out/transcripts_chr.fa.gz" | grep '^>' | tr -d '>' | awk -F'|' 'BEGIN{OFS="\t"}{print $1,$2,$6}' | sort -u > "$out/tx2gene_chr.tsv"
echo "transcripts in annotation:     $(wc -l < "$tmp")"
echo "transcripts in filtered FASTA: $(zcat "$out/transcripts_chr.fa.gz" | grep -c '^>')"
rm -f "$tmp"
