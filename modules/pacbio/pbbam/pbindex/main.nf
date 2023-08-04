process PBINDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbbam=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbbam%3A2.1.0--h3f0f298_2':
        'quay.io/biocontainers/pbbam' }"

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
