process PRIMROSE {

    conda (params.enable_conda ? "bioconda::primrose" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/primrose:1.3.0--h9ee0642_0':
        'quay.io/biocontainers/primrose' }"

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
        primrose $unaligned_bam ${prefix}_5mc.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        primrose: \$(primrose --version | sed 's/primrose v//g')
    END_VERSIONS
        """
}
