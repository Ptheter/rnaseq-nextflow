nextflow.enable.dsl = 2

include { FETCH_ENA                  } from './modules/fetch_ena.nf'
include { FASTQC                     } from './modules/fastqc.nf'
include { FASTP                      } from './modules/fastp.nf'
include { SALMON_INDEX; SALMON_QUANT } from './modules/salmon.nf'
include { MULTIQC                    } from './modules/multiqc.nf'

workflow {
    rows = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .branch { row ->
            ena:   row.srr
            local: true
        }

    FETCH_ENA(
        rows.ena.map { row ->
            tuple([id: row.sample, cell_line: row.cell_line, condition: row.condition], row.srr)
        },
        file(params.ena_metadata)
    )

    local_reads = rows.local.map { row ->
        tuple([id: row.sample],
              [file(row.fastq_1, checkIfExists: true),
               file(row.fastq_2, checkIfExists: true)])
    }

    reads_ch = local_reads.mix(FETCH_ENA.out.reads)

    FASTQC(reads_ch)
    FASTP(reads_ch)

    SALMON_INDEX(file(params.transcriptome, checkIfExists: true))
    SALMON_QUANT(FASTP.out.reads, SALMON_INDEX.out.index)

    multiqc_ch = FASTQC.out.reports.map { meta, files -> files }
        .mix(FASTP.out.json.map { meta, json -> json })
        .mix(SALMON_QUANT.out.quant.map { meta, dir -> dir })
        .flatten()
        .collect()

    MULTIQC(multiqc_ch)
}
