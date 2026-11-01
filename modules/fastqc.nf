process FASTQC {
    tag "$meta.id"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.{zip,html}"), emit: reports

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
