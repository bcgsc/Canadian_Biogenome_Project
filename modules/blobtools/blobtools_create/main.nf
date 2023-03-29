process BLOBTOOLS_CREATE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(config)

    output:
    tuple val(meta), path("${meta.id}"), emit: blobtools_folder
    tuple val(meta), path('*/*.json'), emit: json
    tuple val(meta), path('*/meta.json'), emit:meta_json
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """  
    singularity exec -B /projects blobtoolkit-blobtools_latest.sif blobtools create \
    	--fasta ${assembly} \
    	--meta ${config} \
	--taxid ${params.taxon_taxid} \
	--taxdump BlobtoolkitDatabase/taxdump \
    	${meta.id}
    """
}
