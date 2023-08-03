process GFASTATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gfastats=1.3.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.6--hdcf5f25_3':
        'biocontainers/gfastats:1.3.6--hdcf5f25_3' }"

    input:
    tuple val(meta), path(assembly)   // input.[fasta|fastq|gfa][.gz]

    output:
    tuple val(meta), path("*.assembly_summary"), emit: assembly_summary
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gfastats \\
        $args \\
        --threads $task.cpus \\
        $assembly > ${prefix}.assembly_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """
}
