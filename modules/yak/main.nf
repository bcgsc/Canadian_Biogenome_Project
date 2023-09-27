process YAK {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::yak=0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yak%3A0.1--he4a0461_4' :
        'quay.io/biocontainers/yak%3A0.1--he4a0461_4' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("*.yak"), emit: yak
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # for paired end: to provide two identical streams
    yak count -b37 -t32 -o ${prefix}.yak <(zcat $reads)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yak: \$(echo \$(yak --version 2>&1) | sed 's/^.*yak //; s/Using.*\$//')
    END_VERSIONS
    """
}
