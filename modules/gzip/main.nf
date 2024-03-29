process GZIP {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"


    input:
    tuple val(meta), path(file_to_compress)

    output:
    tuple val(meta), path('*.gz'), emit: gz
    path  "versions.yml"          , emit: versions
 
    script:
    """
    bgzip -c $file_to_compress > ${meta.id}.scaffolds_FINAL.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	bgzip: \$(bgzip --version | sed -n 1p |sed 's/bgzip (htslib) //g')
    END_VERSIONS
    """
}
