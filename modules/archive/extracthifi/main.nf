process EXTRACTHIFI {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::extracthifi=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/extracthifi%3A1.0.0--0':
        'quay.io/biocontainers/extracthifi' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam'), emit: hifi_bam
 
    script:
    """
    extracthifi \\
	$bam \\
	${meta.id}_hifi.bam
    """
}
