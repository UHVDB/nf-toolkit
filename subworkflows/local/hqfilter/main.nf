include { TRTRIMMER                             } from '../../../modules/local/trtrimmer/main'
include { FIND_CONCATENATE                      } from '../../../modules/nf-core/find/concatenate/main'
include { FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER  } from '../../../subworkflows/nf-core/fasta_vclust_prefilter_align_cluster/main'
include { CSVTK_UNIQ                            } from '../../../modules/local/csvtk/uniq/main'

workflow HQFILTER {
    take:
    virus_fna_gz    // channel: [ val(meta), fna ]
    complete_fna_gz // channel: [ val(meta), fna ]
    class_tsv_gz    // channel: [ val(meta), tsv ]
    checkv_db       // channel: [ val(meta), checkv_db ]

    main:
    //
    // MODULE: Trim DTRs from complete sequences
    //
    TRTRIMMER(
        complete_fna_gz
    )

    //
    // MODULE: Combine complete sequences
    //
    FIND_CONCATENATE(
        TRTRIMMER.out.fna_gz.map { _meta, fasta -> [ fasta] }.collect().map { fastas -> [ [ id: 'complete_viruses' ], fastas ] }
    )

    // Define checkv database in case subworkflow is skipped
    ch_checkv_db = checkv_db

    // Define vclust input if therer are enough complete sequences for an update
    ch_vclust_input = FIND_CONCATENATE.out.file_out
        .map { meta, fasta -> [ meta, fasta, fasta.countFasta() ] }
        .filter { _meta, _fasta, count -> count >= params.min_checkv_update || workflow.stubRun }
        .map { meta, fasta, _count -> [ meta, fasta ] }

    //
    // SUBWORKFLOW: Cluster complete viruses for update
    //
    FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER(
        ch_vclust_input,
        false,
        'ani',
        [],
        [],
        0.95
    )

    //
    // MODULE: Extract representative sequence IDs from vclust output
    //
    CSVTK_UNIQ(
        FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER.out.clusters
    )
}
