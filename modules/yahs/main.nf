process YAHS {
    tag "$meta.id"
    label 'process_medium'

    container "https://depot.galaxyproject.org/singularity/yahs%3A1.2a.2--h7132678_0"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta3), path(index)
    tuple val(meta2), path(bam)

    output:
    tuple val(meta), path("*bin"), emit: bin
    tuple val(meta), path("*_scaffolds_final.agp"), emit: agp
    tuple val(meta), path("*_scaffolds_final.fa"), emit: fasta
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    yahs \\
        -o ${assembly.baseName} \\
        $assembly \\
        $bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version | sed 's/yahs v//g')
    END_VERSIONS
    """
}
