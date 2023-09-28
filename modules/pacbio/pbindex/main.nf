process PBINDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbtk==3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.0--h9ee0642_0':
        'quay.io/biocontainers/pbtk' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam.pbi'), emit: index
    path  "versions.yml"          , emit: versions

    script:
    """
    pbindex \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbindex: \$(pbindex --version | sed 's/pbindex //g'| head -n1 )
    END_VERSIONS
    """
}
