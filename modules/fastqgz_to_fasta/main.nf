process FASTQGZ_TO_FASTA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
 
    script:
    """
    bgzip -cd $reads | paste - - - -  | cut -f 1,2 | sed 's/^/>/'  | tr "\t" "\n" | sed 's/@//' > ${meta.id}_pacbio.fasta
    """
}
