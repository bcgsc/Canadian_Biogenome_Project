process BAM2FASTX {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bam2fastx=1.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bam2fastx%3A1.3.1--hf05d43a_1':
        'quay.io/biocontainers/bam2fastx' }"

    input:
    tuple val(meta), path(bam, stageAs: '??.bam'), path(index, stageAs: '??.bam.pbi')

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    bam2fastq \\
        $args \\
	-o ${prefix} \\
	$bam \\
        > ${prefix}.bam2fastx.log
    cat <<-END_VERSIONS > versions.yml
    """
}
