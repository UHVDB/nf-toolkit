/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL PLUGINS/FUNCTIONS/MODULES/SUBWORKFLOWS/WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// FUNCTIONS
def rmEmptyFastAs(ch_fastas) {
    def ch_nonempty_fastas = ch_fastas
        .filter { _meta, fasta ->
            try {
                file(fasta).countFasta( limit: 1 ) > 0
            } catch (java.util.zip.ZipException _e) {
                log.debug "[rmEmptyFastAs]: ${fasta} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                true
            } catch (_EOFException) {
                log.debug "[rmEmptyFastAs]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }
    return ch_nonempty_fastas
}

// MODULES
include { CHECKV_ENDTOEND               } from '../../../modules/local/checkv/endtoend'
include { ENA_GENOMAD                   } from '../../../modules/local/ena/genomad'
include { GENOMAD_ENDTOEND              } from '../../../modules/local/genomad/endtoend'
include { LOCAL_GENOMAD                 } from '../../../modules/local/local/genomad'
include { LOGAN_GENOMAD as ATB_GENOMAD  } from '../../../modules/local/logan/genomad'
include { LOGAN_GENOMAD                 } from '../../../modules/local/logan/genomad'
include { NCBI_GENOMAD                  } from '../../../modules/local/ncbi/genomad'
include { SEQKIT_SEQSPLIT2              } from '../../../modules/local/seqkit/seqsplit2'
include { SPIRE_GENOMAD                 } from '../../../modules/local/spire/genomad'
include { UHBDB_SPLIT                   } from '../../../modules/local/uhbdb/split'
include { UHVDB_VIRUSFILTER             } from '../../../modules/local/uhvdb/virusfilter'
include { UHVDB_CATHEADER               } from '../../../modules/local/uhvdb/catheader'
include { UHVDB_CATNOHEADER             } from '../../../modules/local/uhvdb/catnoheader'
include { VIRALVERIFY_DOWNLOAD          } from '../../../modules/local/viralverify/download'
include { VIRALVERIFY_VIRALVERIFY       } from '../../../modules/local/viralverify/viralverify'

//
// WORFKLOW: Classify viruses in input fasta files
//
workflow CLASSIFY {

    take:
    fna_gz      // channel: [ [ meta ], fna.gz ]
    genomad_db  // channel: [ genomad_db ]
    checkv_db   // channel: [ checkv_db ]

    main:

    // Create channel for DTR sequences file (sequences that had DTRs which were trimmed)
    if ( file("${params.dtr_sequences_file}").exists() ) {
        ch_dtr_sequences = channel.fromPath("${params.dtr_sequences_file}")
    } else {
        ch_dtr_sequences = channel.of([]) 
    }

    //
    // MODULE: Download and store viralverify's database
    //
    VIRALVERIFY_DOWNLOAD()

    // Create channels for combining all geNomad results
    def ch_virus_summaries_tsv_gz = channel.empty()
    def ch_virus_fna_gz = channel.empty()
    def ch_genomad_genes_tsv_gz = channel.empty()

    //
    // MODULE: Download ATB assemblies, run geNomad, and filter
    //
    def ch_atb_assembly_batches = fna_gz.filter { meta, _fasta -> meta.source_db == 'ATB' }
        .map { meta, fasta -> [ meta.id, fasta ] }
        .collate(params.url_split_size * 10)
        .toList()
        .flatMap{ id_fasta -> id_fasta.withIndex() }
        .map { id_fasta, index ->
            [ [ id: 'atb_batch_' + index, source_db: 'ATB' ], id_fasta ]
        }
    ATB_GENOMAD(
        ch_atb_assembly_batches,
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(ATB_GENOMAD.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(ATB_GENOMAD.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(ATB_GENOMAD.out.genes_tsv_gz)

    //
    // MODULE: Download ENA assemblies, run geNomad, and filter
    //
    def ch_ena_urls = fna_gz.filter { meta, _fasta -> meta.source_db == 'ENA' }
        .map { meta, fasta -> [ meta.id, fasta ] }
        .collate(params.url_split_size)
        .toList()
        .flatMap{ id_fasta -> id_fasta.withIndex() }
        .map { id_fasta, index ->
            [ [ id: 'ena_batch_' + index, source_db: 'ENA' ], id_fasta ]
        }
    ENA_GENOMAD(
        ch_ena_urls,
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(ENA_GENOMAD.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(ENA_GENOMAD.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(ENA_GENOMAD.out.genes_tsv_gz)

    //
    // MODULE: Download NCBI Virus assemblies, run geNomad, and filter
    //
    def ch_ncbi_virus_urls = fna_gz.filter { meta, _fasta -> meta.source_db == 'NCBI_VIRUS' }
        .map { meta, fasta -> [ meta.id, fasta ] }
        .collate(params.url_split_size * 10)
        .toList()
        .flatMap{ id_fasta -> id_fasta.withIndex() }
        .map { id_fasta, index ->
            [ [ id: 'ncbi_virus_batch_' + index, source_db: 'NCBI' ], id_fasta ]
        }
    NCBI_GENOMAD(
        ch_ncbi_virus_urls,
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(NCBI_GENOMAD.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(NCBI_GENOMAD.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(NCBI_GENOMAD.out.genes_tsv_gz)


    //
    // MODULE: Split local virus databases into batches
    //
    def ch_local_fastas = fna_gz
        .filter { meta, _fasta ->
            (
                meta.source_db != "ENA" &&
                meta.source_db != "NCBI_VIRUS" &&
                meta.source_db != "LOGAN" &&
                meta.source_db != "SPIRE" &&
                meta.source_db != "ATB" &&
                meta.source_db != "UHBDB" &&
                meta.source_db
            )
        }
        .map { meta, fasta ->
            meta.id = meta.source_db + "_" + meta.id
            [ meta, fasta ]
        }
    SEQKIT_SEQSPLIT2(
        ch_local_fastas,
        params.genomad_split_size
    )

    ch_split_fastas = SEQKIT_SEQSPLIT2.out.fastas_gz
        .filter { _meta, files -> files.size() > 0 }
        .map { _meta, file -> file }
        .flatten()
        .map { file ->
            [ [ id: file.getBaseName().replace(".fasta", "") ], file ]
        }

    //
    // MODULE: Run geNomad on UHBDB assemblies
    //
    def ch_uhbdb_dirs = fna_gz
        .filter { meta, _dir ->
            (
                meta.source_db == 'UHBDB'
            )
        }
    UHBDB_SPLIT(
        ch_uhbdb_dirs
    )
    ch_split_fastas = UHBDB_SPLIT.out.fastas_gz
        .filter { _meta, files -> files.size() > 0 }
        .map { _meta, file -> file }
        .flatten()
        .map { file ->
            [ [ id: file.getBaseName().replace(".fasta", "") ], file ]
        }
        .mix(ch_split_fastas)

    //
    // MODULE: Run geNomad on split local fastas, and filter
    //
    LOCAL_GENOMAD(
        ch_split_fastas,
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(LOCAL_GENOMAD.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(LOCAL_GENOMAD.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(LOCAL_GENOMAD.out.genes_tsv_gz)

    //
    // MODULE: Download LOGAN assemblies, run geNomad, and filter
    //
    def ch_logan_assembly_batches = fna_gz.filter { meta, _fasta -> meta.source_db == 'LOGAN' }
        .map { meta, fasta -> [ meta.id, fasta ] }
        .collate(params.url_split_size * 5)
        .toList()
        .flatMap{ id_fasta -> id_fasta.withIndex() }
        .map { id_fasta, index ->
            [ [ id: 'logan_batch_' + index, source_db: 'LOGAN' ], id_fasta ]
        }
    
    LOGAN_GENOMAD(
        ch_logan_assembly_batches,
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(LOGAN_GENOMAD.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(LOGAN_GENOMAD.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(LOGAN_GENOMAD.out.genes_tsv_gz)

    //
    // MODULE: Download SPIRE assemblies, run geNomad, and filter
    //
    def ch_spire_urls = fna_gz.filter { meta, _fasta -> meta.source_db == 'SPIRE' }
        .map { meta, fasta -> [ meta.id, fasta ] }
        .collate(params.url_split_size)
        .toList()
        .flatMap{ id_fasta -> id_fasta.withIndex() }
        .map { id_fasta, index ->
            [ [ id: 'spire_batch_' + index, source_db: 'SPIRE' ], id_fasta ]
        }
    SPIRE_GENOMAD(
        ch_spire_urls,
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(SPIRE_GENOMAD.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(SPIRE_GENOMAD.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(SPIRE_GENOMAD.out.genes_tsv_gz)

    //
    // MODULE: Run geNomad on local assemblies, and filter
    //
    def ch_no_db_fastas = fna_gz
        .filter { meta, _fasta ->
            (
                !meta.source_db
            )
        }
    GENOMAD_ENDTOEND(
        rmEmptyFastAs(ch_no_db_fastas),
        genomad_db
    )
    ch_virus_summaries_tsv_gz   = ch_virus_summaries_tsv_gz.mix(GENOMAD_ENDTOEND.out.summary_tsv_gz)
    ch_virus_fna_gz             = ch_virus_fna_gz.mix(GENOMAD_ENDTOEND.out.fna_gz)
    ch_genomad_genes_tsv_gz     = ch_genomad_genes_tsv_gz.mix(GENOMAD_ENDTOEND.out.genes_tsv_gz)

    //
    // MODULE: Run CheckV and filter
    //
    CHECKV_ENDTOEND(
        rmEmptyFastAs(ch_virus_fna_gz),
        checkv_db
    )

    //
    // MODULE: Run viralVerify and filter
    //
    VIRALVERIFY_VIRALVERIFY(
        rmEmptyFastAs(CHECKV_ENDTOEND.out.virus_fna_gz),
        VIRALVERIFY_DOWNLOAD.out.nbc_hmms
    )

    //
    // MODULE: Run UHVDB's virus filter
    //
    ch_virus_filter_input = rmEmptyFastAs(CHECKV_ENDTOEND.out.virus_fna_gz)
        .combine(ch_virus_summaries_tsv_gz, by:0)
        .combine(ch_genomad_genes_tsv_gz, by:0)
        .combine(CHECKV_ENDTOEND.out.quality_summary_tsv_gz, by:0)
        .combine(VIRALVERIFY_VIRALVERIFY.out.csv_gz, by:0)

    UHVDB_VIRUSFILTER(
        ch_virus_filter_input,
        ch_dtr_sequences.first()
    )

    //
    // MODULE: Combine UHVDB's virus filter TSVs
    //
    ch_catheader_input = UHVDB_VIRUSFILTER.out.tsv_gz.map { _meta, tsv_gz -> tsv_gz }.collect().map { tsv_gz -> [ [ id:"new_mq_plus_classify" ], tsv_gz, 1, 'tsv.gz' ] }
    UHVDB_CATHEADER(
        ch_catheader_input,
        "${params.output_dir}/outputs/classify/"
    )

    //
    // MODULE: Combine UHVDB's virus filter FASTAs viruses
    //
    ch_catnoheader_input = UHVDB_VIRUSFILTER.out.fna_gz.map { _meta, fasta_files -> fasta_files }.collect().map { fasta_files -> [ [ id:"new_mq_plus_viruses" ], fasta_files, 'fna.gz' ] }
    UHVDB_CATNOHEADER(
        ch_catnoheader_input,
        "${params.output_dir}/outputs/classify/"
    )

    ch_classify_tsv_gz = UHVDB_CATHEADER.out.combined
    ch_mq_viruses_fna_gz = UHVDB_CATNOHEADER.out.combined


    emit:
    classify_tsv_gz     = ch_classify_tsv_gz
    mq_viruses_fna_gz   = ch_mq_viruses_fna_gz
    sample_classify_tsv_gz = UHVDB_VIRUSFILTER.out.tsv_gz
}

