#!/usr/bin/env bash
set -euo pipefail
srr="$1"; tsv="$2"

line=$(awk -F'\t' -v id="$srr" '
  NR==1 { for (i=1; i<=NF; i++) col[$i]=i; next }
  $col["run_accession"]==id { print $col["fastq_ftp"] "\t" $col["fastq_md5"] }
' "$tsv")

if [ -z "$line" ]; then echo "ERROR: $srr not found in $tsv" >&2; exit 1; fi

IFS=';' read -r url1 url2 <<< "$(echo "$line" | cut -f1)"
IFS=';' read -r md51 md52 <<< "$(echo "$line" | cut -f2)"

curl -fsSL --retry 5 --retry-delay 10 -o "${srr}_1.fastq.gz" "https://${url1}"
curl -fsSL --retry 5 --retry-delay 10 -o "${srr}_2.fastq.gz" "https://${url2}"

printf '%s  %s\n' "$md51" "${srr}_1.fastq.gz" "$md52" "${srr}_2.fastq.gz" > md5.txt
md5sum -c md5.txt
