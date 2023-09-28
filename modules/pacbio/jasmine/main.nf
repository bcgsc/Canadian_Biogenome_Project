process JASMINE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbjasmine"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbjasmine:2.0.0--h9ee0642_0':
        'quay.io/biocontainers/pbjasmine' }"

    input :
    tuple val(meta), path (unaligned_bam)

    output :
    tuple val(meta), path ('*_5mc.bam'), emit : cpg_bam
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    jasmine $unaligned_bam ${prefix}_5mc.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasmine: \$( jasmine --version | sed 's/jasmine //g'| head -n1  )
    END_VERSIONS
    """
}
