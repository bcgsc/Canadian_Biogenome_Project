process FLYE_PACBIO_ONT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::flye=2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'quay.io/biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
    tuple val(meta), path(pacbio_reads)
    tuple val(meta2), path(ont_reads)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa.gz")  , emit: gfa
    tuple val(meta), path("*.gv.gz")   , emit: gv
    tuple val(meta), path("*.txt")     , emit: txt
    tuple val(meta), path("*.log")     , emit: log
    tuple val(meta), path("*.json")    , emit: json
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    flye \\
        --iterations 0 \\
        $args \\
        --out-dir . \\
        --threads $task.cpus \\
        --pacbio-raw $pacbio_reads $ont_reads

    flye \\
        $args \\
        --pacbio-raw $pacbio_reads \\
        --resume-from polishing \\
        --out-dir . \\
        --threads $task.cpus

    gzip -c assembly.fasta > ${prefix}.assembly.fasta.gz
    gzip -c assembly_graph.gfa > ${prefix}.assembly_graph.gfa.gz
    gzip -c assembly_graph.gv > ${prefix}.assembly_graph.gv.gz
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv flye.log ${prefix}.flye.log
    mv params.json ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}
