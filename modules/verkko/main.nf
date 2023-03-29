process VERKKO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::verkko=1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verkko:1.2--h64afbab_0' :
        'quay.io/biocontainers/verkko' }"

    input:
    tuple val(meta), path(pacbio_reads)
    tuple val(meta2), path(ont_reads)

    output:
    tuple val(meta), path("*_verkko_assembly.fasta")       , emit: assembly
    tuple val(meta), path("*_homopolymer-compressed.gfa")       , emit: gfa
//    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (params.assembly_secondary_mode == 'pacbio+ont'){
	"""
        verkko \\
            $args \\
            -d ${meta.id} \\
            --hifi $pacbio_reads \\
            --nano $ont_reads

	mv ${meta.id}/assembly.fasta ${meta.id}_verkko_assembly.fasta
	mv ${meta.id}/assembly.homopolymer-compressed.gfa ${meta.id}_homopolymer-compressed.gfa
        """
    } else if (params.assembly_secondary_mode == 'pacbio'){
        """
        verkko \\
            $args \\
            -d ${meta.id} \\
            --hifi $pacbio_reads 

        mv ${meta.id}/assembly.fasta ${meta.id}_verkko_assembly.fasta
        mv ${meta.id}/assembly.homopolymer-compressed.gfa ${meta.id}_homopolymer-compressed.gfa
        """
    } else if (params.assembly_secondary_mode == 'ont'){ 
        """
        verkko \\
            $args \\
            -d ${meta2.id} \\
            --nano $ont_reads

        mv ${meta.id}/assembly.fasta ${meta.id}_verkko_assembly.fasta
        mv ${meta.id}/assembly.homopolymer-compressed.gfa ${meta.id}_homopolymer-compressed.gfa
        """
    } else {
	error "Verkko needs a correct mode : 'pacbio', 'pacbio+ont' or 'ont'"
    }
}
