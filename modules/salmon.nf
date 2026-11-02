process SALMON_INDEX {
    tag "$transcripts"
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'

    input:
    path transcripts

    output:
    path "salmon_index", emit: index

    script:
    """
    gunzip -c ${transcripts} > transcripts.fa
    salmon index \\
        --transcripts transcripts.fa \\
        --index salmon_index \\
        --gencode \\
        --threads ${task.cpus}
    """
}

process SALMON_QUANT {
    tag "$meta.id"
    container 'quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: quant

    script:
    """
    salmon quant \\
        --index ${index} \\
        --libType A \\
        --mates1 ${reads[0]} \\
        --mates2 ${reads[1]} \\
        --validateMappings \\
        --threads ${task.cpus} \\
        --output ${meta.id}
    """
}
