process THREENANOPLOT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nanoplot=1.41.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.0--pyhdfd78af_0' :
        'quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pacbio_bam), path(index)
    tuple val(meta2), path(pacbio_bam2), path(index2)
    tuple val(meta3), path(pacbio_bam3), path(index3)

    output:
    tuple val(meta), path("*.html"), optional: true                , emit: html
    tuple val(meta), path("*.png") , optional: true, emit: png
    tuple val(meta), path("*.txt") , optional: true                , emit: txt
    tuple val(meta), path("*.log")  , optional: true               , emit: log
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    NanoPlot \\
        $args \\
	--prefix ${meta.id} \\
        -t $task.cpus \\
        --ubam $pacbio_bam $pacbio_bam2 $pacbio_bam3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
