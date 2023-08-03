process TIDK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::tidk=0.2.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.31--hdbdd923_2' :
        'quay.io/biocontainers/tidk:0.2.31' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")     , emit: tsv_telomere
    tuple val(meta), path("*.bedgraph"), emit: bedgraph_telomere
    tuple val(meta), path("*.svg")     , emit: plot_telomere
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tidk --version
    tidk search \\
        $args \\
        -o ${prefix} \\
        -e tsv \\
        --string ${params.string_telomere} \\
        --dir . \\
        ${fasta}

    tidk plot -t ${prefix}_telomeric_repeat_windows.tsv
    
    tidk search \\
        $args \\
        -o ${prefix} \\
        -e bedgraph \\
        --string ${params.string_telomere} \\
        --dir . \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidk: \$(tidk --version | sed -e "s/tidk//g")
    END_VERSIONS
    """
}
