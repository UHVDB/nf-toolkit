// include { GENOMAD_DOWNLOADHALLMARKS } from '../../modules/local/genomad/downloadhallmarks/main'
include { CSVTK_FILTER2             } from '../../../modules/local/csvtk/filter2/main'
include { SEQKIT_GREP               } from '../../../modules/nf-core/seqkit/grep/main'
include { HMMER_HMMSEARCH           } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { CSVTK_FILTER2 as CSVTK_FILTER2_HALLMARKS } from '../../../modules/local/csvtk/filter2/main'


workflow HCFILTER {

    take:
    hq_viruses_fna_gz
    classify_tsv_gz

    main:


    // //
    // // MODULE: Download hallmark HMMs and metadata from Genomad
    // //
    // GENOMAD_DOWNLOADHALLMARKS()

    //
    // MODULE: Identify uncertain viruses in UHVDB_CLASSIFY outputs
    //
    CSVTK_FILTER2(
        classify_tsv_gz,
    )

    // Create input for seqkit grep
    ch_seqkit_grep_input = hq_viruses_fna_gz
        .join(CSVTK_FILTER2.out.csv) // join with filtered completenesss by meta
        .multiMap { meta, fasta, tsv ->
            fasta: [ meta, fasta ] 
            tsv: [ tsv ]
        } 

    //
    // MODULE: Extract uncertain viruses from UVHDB_CLASSIFY outputs
    //
    SEQKIT_GREP(
        ch_seqkit_grep_input.fasta,
        ch_seqkit_grep_input.tsv
    )

    //
    // MODULE: Perform HMMsearch for Genomad hallmarks on uncertain viruses
    //
    HMMER_HMMSEARCH(
        SEQKIT_GREP.out.filter.combine()
    )

    //
    // MODULE: Combine HMMsearch output TSVs into a single file with a header
    //
    ch_catheader_input = GENOMAD_HMMSEARCH.out.tsv_gz.map { _meta, tsv_gz -> tsv_gz }.collect().map { tsv_gz -> [ [ id:'new_hcfilter' ], tsv_gz, 1, 'tsv.gz' ] }
    UHVDB_CATHEADER(
        ch_catheader_input,
        "${params.output_dir}/outputs/hcfilter"
    )

    //
    // MODULE: Combine certain and hmm-passing virus FASTA files into a single file
    //
    ch_catnoheader_input = GENOMAD_HMMSEARCH.out.fna_gz.mix( UHVDB_UNCERTAIN.out.certain_fna_gz ).map { _meta, fna_gz -> fna_gz }.collect().map { fna_gz -> [ [ id:'new_hq_hc_viruses' ], fna_gz, 'fna.gz' ] }
    UHVDB_CATNOHEADER(
        ch_catnoheader_input,
        "${params.output_dir}/outputs/hcfilter"
    )

    ch_hcfilter_tsv_gz     = UHVDB_CATHEADER.out.combined
    ch_hq_hc_viruses_fna_gz   = UHVDB_CATNOHEADER.out.combined

    emit:
    hcfilter_tsv_gz             = ch_hcfilter_tsv_gz
    hq_hc_viruses_fna_gz        = ch_hq_hc_viruses_fna_gz
}

