nextflow.enable.dsl = 2

include { FASTQC } from './modules/fastqc.nf'

workflow {
    reads_ch = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple([id: row.sample],
                  [file(row.fastq_1, checkIfExists: true),
                   file(row.fastq_2, checkIfExists: true)])
        }

    FASTQC(reads_ch)
}
