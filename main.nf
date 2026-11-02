nextflow.enable.dsl = 2

include { FASTQC                     } from './modules/fastqc.nf'
include { FASTP                      } from './modules/fastp.nf'
include { SALMON_INDEX; SALMON_QUANT } from './modules/salmon.nf'

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
    FASTP(reads_ch)

    SALMON_INDEX(file(params.transcriptome, checkIfExists: true))
    SALMON_QUANT(FASTP.out.reads, SALMON_INDEX.out.index)
}
