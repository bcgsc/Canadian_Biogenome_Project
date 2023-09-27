process PBMM2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pbmm2==1.12.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.12.0--h9ee0642_0':
        'quay.io/biocontainers/pbmm2:1.12.0--h9ee0642_0' }"

    input :
    tuple val(meta), path (unaligned_bam)
    tuple val(meta2), path (fasta)

    output :
    tuple val(meta), path ('*_aligned.bam'), emit : aligned_bam
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmm2 align  $unaligned_bam $fasta ${prefix}_aligned.bam --sort

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$(pbmm2 --version | sed 's/pbmm2 //g'| head -n1 )
    END_VERSIONS
    """
}
