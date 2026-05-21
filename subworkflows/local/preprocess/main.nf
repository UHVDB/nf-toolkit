include { DEACON_INDEXFETCH } from '../../../modules/local/deacon/indexfetch'
// include { READ_DOWNLOAD     } from '../../../modules/local/read/download'
// include { READ_PREPROCESS   } from '../../../modules/local/read/preprocess'

workflow PREPROCESS {
    take:
    deacon_index_name   // string, name of deacon index
    input_fastqs        // channel: [ [ meta ], [ read1.fastq.gz, read1.fastq.gz? ] ]
    input_sras          // channel: [ [ meta ], sra ]

    main:
    def ch_preprocessed_spring = channel.empty()

    //
    // MODULE: Download deacon index
    //
    DEACON_INDEXFETCH(
        deacon_index_name
    )

    //
    // MODULE: Download, QC, and remove human reads, then compress with spring
    //
    READ_DOWNLOAD(
        input_sras,
        DEACON_INDEXFETCH.out.index.collect()
    )
    ch_preprocessed_spring = READ_DOWNLOAD.out.spring
        .combine(READ_DOWNLOAD.out.pe_count, by:0)
        .map { meta, spring, pe_count ->
            meta.single_end = (pe_count == 1)
            return [ meta, spring ]
        }
        .mix(ch_preprocessed_spring)

    //
    // MODULE: QC reads, remove human reads, and compress with spring
    //
    READ_PREPROCESS(
        input_fastqs,
        DEACON_INDEXFETCH.out.index.collect()
    )
    ch_preprocessed_spring = ch_preprocessed_spring
        .mix(READ_PREPROCESS.out.spring)

    emit:
    preprocessed_spring = ch_preprocessed_spring
}
