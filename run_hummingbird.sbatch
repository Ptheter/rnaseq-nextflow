#!/bin/bash
#SBATCH --job-name=nf-rnaseq
#SBATCH --partition=128x24
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=nextflow_%j.log

# usage: NF_PROFILE=modules sbatch run_hummingbird.sbatch assets/samplesheet_mcf7.csv results_mcf7
SAMPLESHEET="${1:?samplesheet required}"
OUTDIR="${2:-results}"
PROFILE="${NF_PROFILE:-slurm,singularity}"

module load miniconda3 nextflow singularity-ce
export NXF_SINGULARITY_CACHEDIR="$PWD/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" reports
STAMP=$(date +%Y%m%d_%H%M%S)

nextflow run . -profile "$PROFILE" \
    --input "$SAMPLESHEET" \
    --transcriptome reference/transcripts.fa.gz \
    --outdir "$OUTDIR" \
    -with-report   "reports/report_${STAMP}.html" \
    -with-trace    "reports/trace_${STAMP}.txt" \
    -with-timeline "reports/timeline_${STAMP}.html" \
    -resume
