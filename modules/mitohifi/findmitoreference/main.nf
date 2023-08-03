process FIND_MITO_REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/bgruening/mitohifi'

    input:
    tuple val(meta), path(reads_fasta)
    val specie

    output:
    tuple val (meta), path("*.fasta"), emit : reference_fasta
    tuple val (meta), path("*.gb"), emit : reference_gb
    path  "versions.yml"          , emit: versions
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    def VERSION = "3.0.0"
    """
    ${params.singularity_cache}/MitoHiFi/findMitoReference.py \\
    --species "$specie" \\
    --email ${params.email_adress} \\
    --outfolder . \\
    --min_length 16000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MitoHiFi: $VERSION
    END_VERSIONS
    """
}
