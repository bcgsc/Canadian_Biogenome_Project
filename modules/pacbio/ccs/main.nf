process CCS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pbccs"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbccs%3A6.4.0--h9ee0642_0':
        'quay.io/biocontainers/pbccs' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*_postccs.bam'), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    """
    ccs \\
        --report-file ${meta.id}_ccs_report.txt\\
	--all \\
	$bam \\
	${meta.id}_postccs.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ccs: \$(ccs --version | sed 's/ccs v//g')
    END_VERSIONS
    """
}
