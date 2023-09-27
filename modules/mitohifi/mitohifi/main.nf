process MITOHIFI {
    tag "$meta.id"
    label 'process_medium'

    conda '/home/scorreard/miniconda3/envs/mitohifi_v3'

    input:
    tuple val(meta), path(reads_fasta)
    tuple val(meta2), path(reference_fasta)
    tuple val(meta3), path(reference_gb)

    output:
    tuple val (meta), path('final_mitogenome.fasta'), emit : final_mito_fasta
    tuple val (meta), path('final_mitogenome.gb'), emit : final_mito_gb
    tuple val (meta), path('contigs_stats.tsv'), emit : mito_contig_stat
    tuple val (meta), path('*.png'), emit : figures
    path  "versions.yml"          , emit: versions
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    python ${params.singularity_cache}/MitoHiFi/mitohifi.py \\
    $args \\
    -r "$reads_fasta" \\
    -f $reference_fasta \\
    -g $reference_gb \\
    -t 10 

    sed ' 1 s/.*/& [topology=circular] [location=mitochondrion]/' final_mitogenome.fasta > final_mitogenome_tagged.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        mitohifi : \$(python ${params.singularity_cache}/MitoHiFi/mitohifi.py --version | sed 's/MitoHiFi//g')
    END_VERSIONS
    """
}
