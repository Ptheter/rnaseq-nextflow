process FETCH_ENA {
    tag "$meta.id"

    input:
    tuple val(meta), val(srr)
    path ena_tsv

    output:
    tuple val(meta), path("${srr}_{1,2}.fastq.gz"), emit: reads

    script:
    """
    fetch_ena.sh ${srr} ${ena_tsv}
    """
}
