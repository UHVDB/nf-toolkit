process TRTRIMMER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e26a6b2a7696825b05140168093374369b037a05433376788033210c777129d/data' :
        'community.wave.seqera.io/library/csvtk:0.37.0--113625988dd3285d' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tr-trimmer.fna.gz")    , emit: fna_gz

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ### Trim DTRs
    tr-trimmer \\
        ${fasta} \\
        ${args} \\
        > ${prefix}.fna

    ### Compress
    gzip ${prefix}.fna
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fna.gz
    """
}
