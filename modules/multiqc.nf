process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path multiqc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc . --force
    """
}
