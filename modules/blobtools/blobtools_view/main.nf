process BLOBTOOLS_VIEW {
    tag "$meta.id"
    label 'process_high'
    time '10minutes'
    errorStrategy 'retry'
    maxRetries 5

    input:
    tuple val(meta), path(blobtools_folder)

    output:
    tuple val(meta), path('*.png'), emit: png
     path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    singularity exec -B /projects ${params.singularity_cache}/blobtoolkit-blobtools_latest.sif blobtools view \\
        --host http://localhost \\
	--timeout 60 \\
	--ports 8010-8099 \\
	--view blob \\
	--param largeFonts=true \\
	--format png \\
	--out ${meta.id} \\
	${meta.id}

    singularity exec -B /projects ${params.singularity_cache}/blobtoolkit-blobtools_latest.sif blobtools view \\
	--host http://localhost \\
	--timeout 600 \\
	--ports 8010-8099 \\
	--view cumulative \\
	--param largeFonts=true \\
	--format png \\
	--out ${meta.id} \\
	${meta.id}


    singularity exec -B /projects ${params.singularity_cache}/blobtoolkit-blobtools_latest.sif blobtools view \\
        --host http://localhost \\
	--timeout 60 \\
	--ports 8010-8099 \\
	--view snail \\
	--param largeFonts=true \\
	--format png \\
	--out ${meta.id} \\
	${meta.id}

    cp */*.png .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools : \$(singularity exec -B /projects ${params.singularity_cache}/blobtoolkit-blobtools_latest.sif blobtools --version | sed -e "s/blobtoolkit v//g")
    END_VERSIONS
    """
}
