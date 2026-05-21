process SEQKIT_GENOMAD {
    tag "${meta.id}"
    label 'process_high'

    conda ( "${moduleDir}/environment.yml" )
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c39703d881069a69d14e68d12258fd1a93f2a9bd17870f6f6ab39a6b63e094f/data"
    // Singularity: https://wave.seqera.io/view/builds/bd-30e54ab816eb9c63_1?_gl=1*9mafpw*_gcl_au*MTI1MzgxOTA5MC4xNzY4MjM1MzM1
    // TODO: Add Docker

    input:
    tuple val(meta), val(id_files)
    path(genomad_db)

    output:
    tuple val(meta), path("*_virus.fna.gz")        , emit: fna_gz
    tuple val(meta), path("*_virus_summary.tsv.gz"), emit: summary_tsv_gz
    tuple val(meta), path("*_virus_genes.tsv.gz")  , emit: genes_tsv_gz
    tuple val("${task.process}"), val('genomad'), eval("genomad --version 2>&1 | sed 's/^.*geNomad, version //; s/ .*//'"), topic: versions, emit: versions_genomad

    // TODO: Add seqkit version

    script:
    def records   = id_files.collect { id, path -> "${id}\t${path}" }.join('\n')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp

    ### Remove short contigs
    while IFS=\$'\t' read -r sample_id file; do
        seqkit \\
            seq \\
            --threads ${task.cpus} \\
            --min-len 2000 \\
            "\${file}" \\
        | seqkit replace \\
            --threads ${task.cpus} \\
            -p ^ -r "\${sample_id}_" \\
            --out-file "tmp/\${sample_id}.fna.gz"
    done <<'RECORDS'
    ${records}
    RECORDS

    ### Run geNomad
    cat ./*.fna.gz > combined_filtered.fasta.gz

    genomad \\
        end-to-end \\
        combined_filtered.fasta.gz \\
        genomad_results \\
        ${genomad_db} \\
        --threads ${task.cpus} \\
        --splits 5 --relaxed

    ### Save virus outputs
    gzip -c genomad_results/*_summary/*_virus_summary.tsv > ${prefix}_virus_summary.tsv.gz
    gzip -c genomad_results/*_summary/*_virus_genes.tsv > ${prefix}_virus_genes.tsv.gz

    ### Cleanup
    rm -rf tmp/ genomad_results/ combined_filtered.fasta.gz ${prefix}_filtered_genomad.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_virus.fna.gz
    touch ${prefix}_virus_summary.tsv.gz
    touch ${prefix}_virus_genes.tsv.gz
    """
}
