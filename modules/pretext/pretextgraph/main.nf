process PRETEXTGRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextgraph=0.0.6"

    input:
    tuple val(meta), path(pretext_file)
    tuple val(meta2), path(bedgraph)
    val (graph_name)

    output:
    tuple val(meta), path("*_2.pretext"), emit: pretext
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -c $bedgraph > bedgraph.file.gz
    zcat bedgraph.file.gz | PretextGraph -i ${pretext_file} -n "$graph_name" -o ${prefix}_${graph_name}_2.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph: 0.0.6
    END_VERSIONS
    """
}
