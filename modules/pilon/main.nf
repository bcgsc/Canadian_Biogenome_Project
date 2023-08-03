process PILON {
    tag "$meta.id"
    label 'process_medium'

//    conda "bioconda::pilon=1.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pilon:1.24--hdfd78af_0' :
        'quay.io/biocontainers/pilon' }"

    input:
    tuple val(meta), path(assembly)
    val pilon_mode
    tuple val(meta2), path (bam), path (bai)

    output:
    tuple val(meta), path("*.fasta")       , emit: improved_assembly
    tuple val(meta), path("*.change")       , emit: change_record, optional : true
    tuple val(meta), path("*.vcf")       , emit: vcf, optional : true
    tuple val(meta), path("*.bed")       , emit: tracks_bed, optional : true
    tuple val(meta), path("*.wig")       , emit: tracks_wig, optional : true
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["--frags", "--jumps", "--unpaired", "--bam"]
    if ( !valid_mode.contains(pilon_mode) )  { error "Unrecognised mode to run Pilon. Options: ${valid_mode.join(', ')}" }
    """
    pilon \\
	$args \\
	--genome $assembly \\
	--output ${meta.id} \\
	--threads $task.cpus \\
	$pilon_mode $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(pilon --version | sed 's/pilon v//g')
    END_VERSIONS
    """
}
