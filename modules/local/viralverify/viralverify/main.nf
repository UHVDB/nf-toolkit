process VIRALVERIFY_VIRALVERIFY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/472ca345aeaaaae70034c8bc04e1d62b226e188ec9edc16b7dfd43f83a25ed8e/data'
        : 'community.wave.seqera.io/library/viralverify:1.1--3527fad4ef7ee6f4'}"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path ("${meta.id}_viralverify.csv.gz") , emit: csv_gz
    tuple val("${task.process}"), val('viralverify'), val('1.1'), emit: versions_viralverify, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    script:
    """
    ### Uncompres
    gunzip -c ${fasta} > ${meta.id}.fasta

    ### Run ViralVerify
    viralverify \\
        -f ${meta.id}.fasta \\
        --hmm ${db} \\
        -o ${meta.id}_viralverify \\
        -t ${task.cpus}

    ### Compress
    mv ${meta.id}_viralverify/${file(meta.id + ".fasta").getBaseName()}_result_table.csv ${meta.id}_viralverify.csv
    gzip ${meta.id}_viralverify.csv
    """

    stub:
    """
    ### Touch empty output file
    echo "" | gzip > ${meta.id}_viralverify.csv.gz
    """
}
