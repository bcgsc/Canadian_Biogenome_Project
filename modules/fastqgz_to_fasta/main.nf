process FASTQGZ_TO_FASTA {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path  "versions.yml"          , emit: versions
 
    script:
    """
    bgzip -cd $reads | paste - - - -  | cut -f 1,2 | sed 's/^/>/'  | tr "\t" "\n" | sed 's/@//' > ${meta.id}_pacbio.fasta

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^bgzip (htslib) //; s/ Copyright.*//')
END_VERSIONS
    """
}
