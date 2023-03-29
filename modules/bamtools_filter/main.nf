process BAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bamtools=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hd03093a_0' :
        'quay.io/biocontainers/bamtools:2.5.2--hd03093a_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam'), emit: filtered_bam
    path  "versions.yml"                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    bamtools \\
        filter \\
        -in $bam \\
        $args \\
        -out ${prefix}_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
