process FASTP_DEACON_SPRING {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07f2814c3c48234ded08b489e77986b484886bf3a4b393b40017e7109d623d1b/data' :
        'community.wave.seqera.io/library/deacon_fastp_spring:cf5415ca9c1df5a8' }"
    // xsra is installed in the container via cargo (not available on bioconda)

    input:
    tuple val(meta) , path(fastq)
    path(index)

    output:
    tuple val(meta), path("*.spring*") , emit: spring
    tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | sed -e "s/fastp //g"'), emit: versions_fastp, topic: versions
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions
    tuple val("${task.process}"), val('spring'), val('1.1.1'), topic: versions, emit: versions_spring
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastp_reads_in      = meta.single_end ? "--in1 ${fastq}" : "--in1 ${fastq[0]} --in2 ${fastq[1]}"
    def fastp_reads_out     = meta.single_end ? "--out1 ${prefix}.fastp.fastq.gz" : "--out1 ${prefix}_1.fastp.fastq.gz --out2 ${prefix}_2.fastp.fastq.gz"
    def deacon_reads_in     = meta.single_end ? "${prefix}.fastp.fastq.gz" : "${prefix}_1.fastp.fastq.gz ${prefix}_2.fastp.fastq.gz"
    def deacon_reads_out    = meta.single_end ? "--output ${prefix}.deacon.fastq.gz" : "--output ${prefix}_1.deacon.fastq.gz --output2 ${prefix}_2.deacon.fastq.gz"
    def spring_input        = meta.single_end ? "${prefix}.deacon.fastq.gz" : "${prefix}_1.deacon.fastq.gz ${prefix}_2.deacon.fastq.gz"
    """
    ### Run fastp
    fastp \\
        ${fastp_reads_in} \\
        ${fastp_reads_out} \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe

    ### Run deacon
    deacon filter \\
        --deplete \\
        ${index} \\
        ${deacon_reads_in} \\
        ${deacon_reads_out} \\
        --threads ${task.cpus}

    rm -rf *.fastp.fastq.gz

    ### Run spring
    spring \\
        --compress \\
        --input ${spring_input} \\
        --num-threads ${task.cpus} \\
        --quality-opts ill_bin \\
        --gzipped-fastq \\
        --allow-read-reordering \\
        --output-file ${prefix}.spring

    ### Cleanup
    rm -rf ${prefix}*deacon*.fastq.gz *.fastp.html *.fastp.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.spring
    """
}
