process HICANU {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::canu=2.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu%3A2.2--ha47f30e_0' :
        'quay.io/biocontainers/canu' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.report")       	, emit: report
    tuple val(meta), path("*.contigs.fasta")	, emit: assembly
    tuple val(meta), path("*.unassembled.fasta"), emit: contigs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
	"""
	canu \\
		-p ${meta.id} \\
		-pacbio-hifi \\
		$args \\
		$reads
	"""
}
