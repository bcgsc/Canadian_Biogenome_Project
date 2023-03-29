process BLOBTOOLS_VIEW_BLOB {
    tag "$meta.id"
    label 'process_high'
    time '10minutes'

    input:
    tuple val(meta), path(blobtools_folder)

    output:
    tuple val(meta), path('*.png'), emit: png
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    singularity exec -B /projects blobtoolkit-blobtools_latest.sif blobtools view \\
        --host http://localhost \\
	--timeout 60 \\
	--ports 8010-8099 \\
	--view blob \\
	--param largeFonts=true \\
	--format png \\
	--out ${meta.id} \\
	${meta.id}

    cp */*.png .
    """
}
