process ARIA2C_SEQKIT_GENOMAD {
    tag "${meta.id}"
    label 'process_high'

    conda ( "${moduleDir}/environment.yml" )
    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c39703d881069a69d14e68d12258fd1a93f2a9bd17870f6f6ab39a6b63e094f/data"
    // Singularity: https://wave.seqera.io/view/builds/bd-30e54ab816eb9c63_1?_gl=1*9mafpw*_gcl_au*MTI1MzgxOTA5MC4xNzY4MjM1MzM1
    // TODO: Add Docker

    input:
    tuple val(meta), val(id_urls)
    path(genomad_db)

    output:
    tuple val(meta), path("*_virus.fna.gz")        , emit: fna_gz
    tuple val(meta), path("*_virus_summary.tsv.gz"), emit: summary_tsv_gz
    tuple val(meta), path("*_virus_genes.tsv.gz")  , emit: genes_tsv_gz

    tuple val("${task.process}"), val('genomad'), eval("genomad --version 2>&1 | sed 's/^.*geNomad, version //; s/ .*//'"), topic: versions, emit: versions_genomad
    // TODO: Add aria2c version
    // TODO: Add seqkit version

    script:
    def url_list   = id_urls.collect { id_url -> id_url[1].toString() + ',\sout=' + id_url[0].toString() + '.fna.gz' }.join(',')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ### Create arrays
    mkdir -p tmp
    IFS=',' read -r -a download_array <<< "${url_list}"
    printf '%s\\n' "\${download_array[@]}" > aria2_file.tsv

    ### Download assemblies
    for try in {1..6}; do
        aria2c \\
            --input=aria2_file.tsv \\
            --dir=tmp/ \\
            --max-connection-per-server=${task.cpus} \\
            --split=${task.cpus} \\
            --max-tries=10 \\
            --retry-wait=30 \\
            --max-concurrent-downloads=${task.cpus} && break || sleep \$((\$try^2*60))
    done

    ### Remove short contigs
    for file in tmp/*.fna.gz; do
        sample_id=\$(basename \${file} .fna.gz)

        seqkit \\
            seq \\
            --threads ${task.cpus} \\
            --min-len 2000 \\
            \$file \\
        | seqkit replace \\
            --threads ${task.cpus} \\
            -p ^ -r "\${sample_id}_" \\
            --out-file \${sample_id}.fna.gz

        rm \$file
    done

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
