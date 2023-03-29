process BLOBTOOLS_CREATE {
    tag "$meta.id"
    label 'process_high'

//    container '/projects/cbp/scratch/singularity/blobtoolkit_latest.sif'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(config)

    output:
    tuple val(meta), path("${meta.id}"), emit: blobtools_folder
//    tuple val(meta), path('blobtools_folder'), emit: blobtools_folder
    tuple val(meta), path('*/*.json'), emit: json
    tuple val(meta), path('*/meta.json'), emit:meta_json
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    #source activate blobtools 

    #cp -r /projects/cbp/scratch/Arctic_surfclam_005/hifiasm_purgedups_yahs/blobtoolkit-2.0.0 .

    #./blobtoolkit-2.0.0/blobtools create \\
#	--fasta ${assembly} \\
#	datasets/${meta.id}
  
    singularity exec -B /projects /projects/cbp/scratch/singularity/blobtoolkit-blobtools_latest.sif blobtools create \
    	--fasta ${assembly} \
    	--meta ${config} \
	--taxid ${params.taxon_taxid} \
	--taxdump /projects/CanSeq/BlobtoolkitDatabase/taxdump \
    	${meta.id}
    """
}
