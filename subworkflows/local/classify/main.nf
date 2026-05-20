include { rmEmptyFastAs; rmEmptyTsvs; add_split; extractDigitBeforeExtension } from '../functions/main'
include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { GENOMAD_DOWNLOAD          } from '../../../modules/nf-core/genomad/download/main'
include { VIRALVERIFY_DOWNLOAD      } from '../../../modules/local/viralverify/download/main'
include { SEQKIT_SEQ                } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_REPLACE            } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_SPLIT2             } from '../../../modules/nf-core/seqkit/split2/main'
include { GENOMAD_ENDTOEND          } from '../../../modules/nf-core/genomad/endtoend/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_PROVIRUS } from '../../../modules/nf-core/seqkit/replace/main'
include { CAT_CAT as CAT_CAT_CHECKV } from '../../../modules/nf-core/cat/cat/main'
include { CSVTK_FILTER2             } from '../../../modules/local/csvtk/filter2/main'
include { SEQKIT_GREP               } from '../../../modules/nf-core/seqkit/grep/main'
include { VIRALVERIFY_VIRALVERIFY   } from '../../../modules/local/viralverify/viralverify/main'
include { UHVDB_CLASSIFY            } from '../../../modules/local/uhvdb/classify/main'
include { CAT_CAT                   } from '../../../modules/nf-core/cat/cat/main'

workflow CLASSIFY {
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
    // MODULE: Length filter fasta files
    //
    SEQKIT_SEQ(
        fastas
    )

    //
    // MODULE: Re-name Fasta sequences to include the sample_id as the prefix
    //
    SEQKIT_REPLACE(
        SEQKIT_SEQ.out.fastx
    )

    // Identify channels where the fasta has > params.split_size sequences
    ch_seqkit_replace_out_fastx = SEQKIT_REPLACE.out.fastx
        .branch { _meta, fastx ->
            split: file(fastx).countFasta( limit: params.split_size + 1 ) > params.split_size
            nosplit: true
        }

    //
    // MODULE: Split sequences into chunks of size params.chunk_size if they are larger than this
    //
    SEQKIT_SPLIT2(
        ch_seqkit_replace_out_fastx.split
    )
    ch_seqkit_split2_reads = SEQKIT_SPLIT2.out.reads
        .transpose()
        .map{ meta, fastx ->
            [add_split(meta, fastx.getName()), [fastx]]}

    //
    // MODULE: Run geNomad end-to-end
    //
    GENOMAD_ENDTOEND(
        ch_seqkit_replace_out_fastx.nosplit.mix(ch_seqkit_split2_reads),
        ch_genomad_db
    )

    //
    // MODULE: Run CheckV end-to-end
    //
    CHECKV_ENDTOEND(
        rmEmptyFastAs(GENOMAD_ENDTOEND.out.virus_fasta),
        ch_checkv_db
    )

    //
    // MODULE: Rename CheckV proviruses
    //
    SEQKIT_REPLACE_PROVIRUS(
        rmEmptyFastAs(CHECKV_ENDTOEND.out.proviruses)
    )

    //
    // MODULE: Combine CheckV viruses and proviruses
    //
    CAT_CAT_CHECKV(
        CHECKV_ENDTOEND.out.viruses,
        SEQKIT_REPLACE_PROVIRUS.out.fastx
    )

    //
    // MODULE: Filter CheckV's completeness output
    //
    CSVTK_FILTER2(
        CHECKV_ENDTOEND.out.completeness
    )

    // Get only unique sequences in hq fasta
    ch_mq_plus_viruses = CSVTK_FILTER2.out.csv
        .map { _meta, file -> file.text.readLines() }
        .flatten()
        .map { line -> line.split('\t')[0] }
        .unique()
        .collectFile(name: 'mq_plus_viruses.txt', newLine: true)

    //
    // MODULE: Filter CheckV output sequences
    //
    SEQKIT_GREP(
        rmEmptyFastAs(CHECKV_ENDTOEND.out.viruses.mix(SEQKIT_REPLACE.out.fastx)),
        ch_mq_plus_viruses
    )

    //
    // MODULE: Run ViralVerify
    //
    VIRALVERIFY_VIRALVERIFY(
        SEQKIT_GREP.out.filter,
        ch_viralverify_db
    )

    // Combine geNomad, CheckV, and viralverify outputs
    ch_uhvdb_classify_input = SEQKIT_GREP.out.filter
        .combine(SEQKIT_GREP.out.filter, by:0)
        .combine(GENOMAD_ENDTOEND.out.virus_summary, by:0)
        .combine(GENOMAD_ENDTOEND.out.virus_genes, by:0)
        .combine(CHECKV_ENDTOEND.out.quality_summary, by:0)
        .combine(rmEmptyTsvs(VIRALVERIFY_VIRALVERIFY.out.result_table), by:0)
        .view()

    //
    // MODULE: Classify viruses by combining results from geNomad, CheckV, and ViralVerify
    //
    UHVDB_CLASSIFY(
        ch_uhvdb_classify_input,
        dtr_sequences_file
    )

    // //
    // // MODULE: Combine filtered virus fastas
    // //
    // CAT_CAT(

    // )


    emit:
    virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz
    virus_fna_gz             = ch_virus_fna_gz
    genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz
}
