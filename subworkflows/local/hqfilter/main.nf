include { TRTRIMMER } from '../../../modules/local/trtrimmer/main'
include { countFastAs } from '../../../subworkflows/local/functions/main'

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

    // If not enough complete sequences, use current CheckV database
    if ( countFastAs(TRTRIMMER.out.fna_gz) < params.min_checkv_update ) {
        ch_checkv_db = checkv_db
    } else {
        VCLUST_
    }

}
