include { rmEmptyFastAs; rmEmptyTsvs; add_split; extractDigitBeforeExtension } from '../functions/main'
include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { GENOMAD_DOWNLOAD          } from '../../../modules/nf-core/genomad/download/main'
include { VIRALVERIFY_DOWNLOAD      } from '../../../modules/local/viralverify/download/main'
include { SEQKIT_SEQ_REPLACE_SPLIT2 } from '../../../modules/local/seqkit/seq_replace_split2/main'
include { GENOMAD_ENDTOEND          } from '../../../modules/nf-core/genomad/endtoend/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_PROVIRUS } from '../../../modules/nf-core/seqkit/replace/main'
include { FIND_CONCATENATE as FIND_CONCATENATE_CHECKV } from '../../../modules/nf-core/find/concatenate/main'
include { CSVTK_FILTER2             } from '../../../modules/local/csvtk/filter2/main'
include { SEQKIT_GREP               } from '../../../modules/nf-core/seqkit/grep/main'
include { VIRALVERIFY_VIRALVERIFY       } from '../../../modules/local/viralverify/viralverify/main'
include { UHVDB_CLASSIFY                } from '../../../modules/local/uhvdb/classify/main'
include { FIND_CONCATENATE          } from '../../../modules/nf-core/find/concatenate/main'

workflow CLASSIFY {

    take:
    fastas              // channel: [ val(meta), fna ]
    dtr_sequences_file  // string, DTR sequences file

    main:
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

    // Split fasta files based on source_type
    fastas = fastas.branch { meta, _fasta ->
        assembly: meta.source_type == 'Assembly'
        database: meta.source_type == 'Database'
    }

    //
    // MODULE: Length filter, add sample_id as prefix, and split database fasta files
    //
    SEQKIT_SEQ_REPLACE_SPLIT2(
        fastas.database
    )
    ch_split_database_fastas = SEQKIT_SEQ_REPLACE_SPLIT2.out.fastx
        .transpose()
        .map{ meta, fastx ->
            [add_split(meta, fastx.getName()), [fastx]]}

    //
    // MODULE: Run geNomad end-to-end on local database fasta files
    //
    GENOMAD_ENDTOEND(
        ch_split_database_fastas,
        ch_genomad_db
    )

    // Identify remote fasta files and split them into chunks of size params.url_split_size
    ch_assembly_fastas = fastas.assembly
        .map { meta, fasta ->
            def meta_new = [:]
            meta_new.body_site = meta.body_site    
            meta_new.source_type = meta.source_type
            meta_new.source_db = meta.source_db
            return [ meta_new, fasta ]
        }
        .groupTuple() // Combine URLs by source_db, source_type, and body_site
        .map { meta, fasta_list ->
            fasta_list
                .collate(params.url_split_size) // Split fasta files into chunks of size params.url_split_size
                .withIndex() // Add index to each chunk
                .collect { batch, idx ->
                    def batch_meta = meta + [ id: "${meta.source_db}_batch_${idx}" ]
                    [ batch_meta, batch ] // Create a new meta with the source_db and batch index
                }
        }
        .flatMap { batches -> batches }   // one [meta, batch] per emission
        .branch { _meta, fasta -> 
            remote: file(fasta[0]).toUri().toString().startsWith('https://')
            local: true
        }

    ch_assembly_fastas.remote.view()
    ch_assembly_fastas.local.view()

    

    // //
    // // MODULE: Run CheckV end-to-end
    // //
    // CHECKV_ENDTOEND(
    //     rmEmptyFastAs(GENOMAD_ENDTOEND.out.virus_fasta),
    //     ch_checkv_db
    // )

    // //
    // // MODULE: Rename CheckV proviruses
    // //
    // SEQKIT_REPLACE_PROVIRUS(
    //     rmEmptyFastAs(CHECKV_ENDTOEND.out.proviruses)
    // )

    // // Create a channel to join CheckV viruses and proviruses when necessary
    // ch_checkv_viruses = CHECKV_ENDTOEND.out.viruses
    //     .join(SEQKIT_REPLACE_PROVIRUS.out.fastx, remainder: true)
    //     .branch { _meta, _viruses, proviruses ->
    //         cat: proviruses != null
    //         no_cat: true
    //     }

    // //
    // // MODULE: Combine CheckV viruses and proviruses
    // //
    // FIND_CONCATENATE(
    //     ch_checkv_viruses.cat
    // )

    // //
    // // MODULE: Filter CheckV's completeness output
    // //
    // CSVTK_FILTER2(
    //     CHECKV_ENDTOEND.out.completeness
    // )

    // // Create input for seqkit grep
    // ch_seqkit_grep_input = ch_checkv_viruses.no_cat
    //     .map { meta, viruses, _proviruses -> [ meta, viruses ] } // Remove proviruses item when it is null
    //     .mix(FIND_CONCATENATE.out.file_out) // Combine with virus + provirus joined fasta
    //     .join(CSVTK_FILTER2.out.csv) // join with filtered completenesss by meta
    //     .multiMap { meta, fasta, tsv ->
    //         fasta: [ meta, fasta ]
    //         tsv: [ tsv ]
    //     }

    // //
    // // MODULE: Filter CheckV output sequences
    // //
    // SEQKIT_GREP(
    //     ch_seqkit_grep_input.fasta,
    //     ch_seqkit_grep_input.tsv
    // )

    // //
    // // MODULE: Run ViralVerify
    // //
    // VIRALVERIFY_VIRALVERIFY(
    //     SEQKIT_GREP.out.filter,
    //     ch_viralverify_db
    // )

    // // Combine geNomad, CheckV, and viralverify outputs
    // ch_uhvdb_classify_input = SEQKIT_GREP.out.filter
    //     .join(GENOMAD_ENDTOEND.out.virus_summary)
    //     .join(GENOMAD_ENDTOEND.out.virus_genes)
    //     .join(CHECKV_ENDTOEND.out.quality_summary)
    //     .join(rmEmptyTsvs(VIRALVERIFY_VIRALVERIFY.out.csv_gz))

    // //
    // // MODULE: Classify viruses by combining results from geNomad, CheckV, and ViralVerify
    // //
    // UHVDB_CLASSIFY(
    //     ch_uhvdb_classify_input,
    //     dtr_sequences_file
    // )

    // emit:
    // checkv_db       = ch_checkv_db
    // virus_fna_gz    = UHVDB_CLASSIFY.out.fna_gz
    // complete_fna_gz = UHVDB_CLASSIFY.out.complete_fna_gz
    // class_tsv_gz    = UHVDB_CLASSIFY.out.tsv_gz
}
