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
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    /projects/cbp/scratch/singularity/MitoHiFi/findMitoReference.py \\
    --species "$specie" \\
    --email $params.email \\
    --outfolder . \\
    --min_length 16000
    """
}
