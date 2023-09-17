process FOURBAM2FASTX {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bam2fastx=1.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bam2fastx%3A1.3.1--hf05d43a_1':
        'quay.io/biocontainers/bam2fastx' }"

    input:
    tuple val(meta), path(bam), path(index)
    tuple val(meta2), path(bam2), path(index2)
    tuple val(meta3), path(bam3), path(index3)
    tuple val(meta4), path(bam4), path(index4)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path  "versions.yml"          , emit: versions
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    bam2fastq \\
        $args \\
        -o ${prefix} \\
        $bam \\
        $bam2 \\
        $bam3 \\
        $bam4 \\
        > ${prefix}.bam2fastx.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastq: \$(bam2fastq --version | sed 's/bam2fastq //g')
    END_VERSIONS
    """
}
