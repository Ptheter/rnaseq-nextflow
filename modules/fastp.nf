process FASTP {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed_R*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),                emit: json
    path "*.html",                                  emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread ${task.cpus}
    """
}
