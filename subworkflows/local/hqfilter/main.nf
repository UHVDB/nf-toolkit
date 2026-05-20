include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { GENOMAD_DOWNLOAD          } from '../../../modules/nf-core/genomad/download/main'
include { VIRALVERIFY_DOWNLOAD      } from '../../../modules/local/viralverify/download/main'
include { GENOMAD_ENDTOEND          } from '../../../modules/nf-core/genomad/endtoend/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'
include { SEQKIT_REPLACE            } from '../../../modules/nf-core/seqkit/replace/main'
include { CSVTK_FILTER2             } from '../../../modules/local/csvtk/filter2/main'
include { VIRALVERIFY_VIRALVERIFY   } from '../../../modules/local/viralverify/viralverify/main'

workflow HQFILTER {
    take:
    fastas              // channel: [ val(meta), fna ]
    dtr_sequences_file  // string, DTR sequences file

    main:

    // Create channels for combining all geNomad results
    def ch_virus_summaries_tsv_gz   = channel.empty()
    def ch_virus_fna_gz             = channel.empty()
    def ch_genomad_genes_tsv_gz     = channel.empty()

    //
    // MODULE: Download genomad's database
    //
    if ( params.genomad_db ) {
        ch_genomad_db = channel.fromPath("${params.genomad_db}").first()
    } else {
        GENOMAD_DOWNLOAD()
        ch_genomad_db = GENOMAD_DOWNLOAD.out.genomad_db
    }

    //
    // MODULE: Download checkv's database
    //
    if ( params.checkv_db ) {
        ch_checkv_db = channel.fromPath("${params.checkv_db}").first()
    } else {
        CHECKV_DOWNLOADDATABASE()
        ch_checkv_db = CHECKV_DOWNLOADDATABASE.out.checkv_db
    }

    //
    // MODULE: Download viralverify's database
    //
    if ( params.viralverify_db ) {
        ch_viralverify_db = channel.fromPath("${params.viralverify_db}").first()
    } else {
        VIRALVERIFY_DOWNLOAD()
        ch_viralverify_db = VIRALVERIFY_DOWNLOAD.out.viralverify_db
    }

    //
    // MODULE: Run geNomad end-to-end
    //
    GENOMAD_ENDTOEND(
        fastas,
        ch_genomad_db
    )

    //
    // MODULE: Run CheckV end-to-end
    //
    CHECKV_ENDTOEND(
        GENOMAD_ENDTOEND.out.virus_fasta,
        ch_checkv_db
    )

    //
    // MODULE: Rename CheckV proviruses
    //
    SEQKIT_REPLACE(
        CHECKV_ENDTOEND.out.proviruses
    )

    //
    // MODULE: Filter CheckV's completeness output
    //
    CSVTK_FILTER2(
        CHECKV_ENDTOEND.out.completeness
    )

    //
    // MODULE: Filter CheckV output sequences
    //
    SEQKIT_GREP(
        CSVTK_FILTER2.out.csv.combine(CSVTK_FILTER2.out.csv, by:0)
    )

    //
    // MODULE: Run ViralVerify
    //
    VIRALVERIFY_VIRALVERIFY(
        CHECKV_ENDTOEND.out.viruses.mix(CHECKV_ENDTOEND.out.proviruses),
        ch_viralverify_db
    )

    //
    // MODULE: Classify viruses by combining results from geNomad, CheckV, and ViralVerify
    //
    // ch_virus_filter_input

    emit:
    virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz
    virus_fna_gz             = ch_virus_fna_gz
    genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz
}
