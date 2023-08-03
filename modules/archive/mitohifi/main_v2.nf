process MITOHIFI {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/bgruening/mitohifi'
//	container 'ghcr.io/marcelauliano/MitoHiFi:main'

    input:
    tuple val(meta), path(assembly)
    val specie

    output:
	path('final_mitogenome.fasta', emit : final_mito_fasta)
	path('final_mitogenome.gb', emit : final_mito_gb)
	path('contigs_stats.tsv', emit : mito_contig_stat)
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    ${params.singularity_cache}/MitoHiFi/findMitoReference.py \\
    --species "$specie" \\
    --email ${params.email_adress} \\
    --outfolder . \\
    --min_length 16000

    python ${params.singularity_cache}/MitoHiFi/mitohifi.py \\
    $args \\
    -c "$assembly" \\
    -f *.fasta \\
    -g *.gb \\
    -t 10  
    """
}
