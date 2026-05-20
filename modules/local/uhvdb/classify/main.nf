process UHVDB_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff2c7e6a1d237f65056929cb69cf5301742780ff22a7f79d5b763c72874b8de6/data':
        'community.wave.seqera.io/library/biopython_polars:1c1d88559d24ac35' }"

    input:
    tuple val(meta), path(virus_fasta), path(provirus_fasta), path(virus_summary), path(genes), path(quality_summary), path(viralverify)
    path(dtr_sequences_txt)

    output:
    tuple val(meta), path("${meta.id}.uhvdb_viruses.fna.gz")    , emit: fna_gz
    tuple val(meta), path("${meta.id}.uhvdb_virus_class.tsv.gz"), emit: tsv_gz
    tuple val("${task.process}"), val('polars'), eval('python -c "import polars; print(polars.__version__)"'), topic: versions, emit: versions_polars
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), topic: versions, emit: versions_biopython
    tuple val("${task.process}"), val('uhvdb_classify'), eval('python -c "import uhvdb_classify; print(uhvdb_classify.__version__)"'), topic: versions, emit: versions_uhvdb_classify

    when:
    task.ext.when == null || task.ext.when

    script:
    def source_db = meta.source_db ?: 'no_source_db'
    def dtr_sequences = dtr_sequences_txt ? "--dtr_sequences ${dtr_sequences_txt}" : "--dtr_sequences ''"
    """
    ### Combine fasta files
    cat ${virus_fasta} ${provirus_fasta} > ${meta.id}.combined_viruses.fna

    ### Run uhvdb_classify
    uhvdb_virus_filter.py \\
        --fasta ${meta.id}.combined_viruses.fna \\
        --virus_summary ${virus_summary} \\
        --genes ${genes} \\
        --quality_summary ${quality_summary} \\
        --viralverify ${viralverify} \\
        ${dtr_sequences} \\
        --output_fasta ${meta.id}.uhvdb_viruses.fna \\
        --output_tsv ${meta.id}.uhvdb_virus_class.tsv \\
        --source_db ${source_db} \\
        --db_type ${meta.db_type} \\
        --body_site ${meta.body_site}

    ### Compress
    gzip ${meta.id}.uhvdb_viruses.fna
    gzip ${meta.id}.uhvdb_virus_class.tsv

    ### Cleanup
    rm -f ${meta.id}.combined_viruses.fna
    """

    stub:
    """
    ### Touch empty output files
    echo "" | gzip > ${meta.id}.uhvdb_viruses.fna.gz
    echo "" | gzip > ${meta.id}.uhvdb_virus_class.tsv.gz
    """
}