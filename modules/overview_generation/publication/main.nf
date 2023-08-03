process PUBLICATION {
	tag "$meta.id"
    label 'process_low'

//container = 'https://depot.galaxyproject.org/singularity/r-rjson%3A0.2.15--r3.3.2_0'
        conda "r-officer"



    input:
    tuple val(meta), path(overview_sample) //LonQC

    output :
    tuple val(meta), path('*.docx'), emit: genome_note

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    Rscript ${params.modules_path}/overview_generation/publication/publication_officer.R \\
	$overview_sample
    """
}
