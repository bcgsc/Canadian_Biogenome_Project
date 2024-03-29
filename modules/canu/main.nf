process CANU {
    tag "$meta.id"
    label 'process_high'

//    conda "bioconda::canu=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'quay.io/biocontainers/canu:2.2--ha47f30e_0' }"

    input:
//As Canu doesn't take pacbio+ont, it can only be one or the other type
    tuple val(meta), path(reads)
    val(genome_size)

    output:
    tuple val(meta), path("*.report")                   , emit: report
    tuple val(meta), path("*.contigs.fasta.gz")         , emit: assembly                , optional: true
    tuple val(meta), path("*.unassembled.fasta.gz")     , emit: contigs
    tuple val(meta), path("*.correctedReads.fasta.gz")	, emit: corrected_reads         , optional: true
    tuple val(meta), path("*.trimmedReads.fasta.gz")	, emit: corrected_trimmed_reads , optional: true
    tuple val(meta), path("*.contigs.layout")           , emit: metadata                , optional: true
    tuple val(meta), path("*.contigs.layout.readToTig") , emit: contig_position         , optional: true
    tuple val(meta), path("*.contigs.layout.tigInfo")   , emit: contig_info             , optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genomesize = genome_size ? "genomesize=$genome_size" : "" 
    if (params.assembly_secondary_mode == 'hicanu'){
    """
    canu \\
        -p ${prefix} \\
        $args \\
        $genomesize \\
        maxThreads=$task.cpus \\
        -pacbio-hifi $reads
                
    gzip *.fasta
            
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """     
    } else if (params.assembly_secondary_mode == 'ont'){
    """
    canu \\
        -p ${prefix} \\
        $args \\
        $genomesize \\
        maxThreads=$task.cpus \\
        -nanopore $reads

    gzip *.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """
    } else if (params.assembly_secondary_mode == 'clr'){
    """
    canu \\
        -p ${prefix} \\
        $args \\
        $genomesize \\
        maxThreads=$task.cpus \\
        -pacbio $reads

    gzip *.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """
    } else {
	error "Canu needs a correct mode : 'hicanu', 'ont' or 'clr'"
    }
}
