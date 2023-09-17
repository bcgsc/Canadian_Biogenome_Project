process COVERAGE_CALCULATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq)
    val (genome_size)

    output:
    path('*.txt'), emit: coverage
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    gzip -cd $fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c | awk '{ print "Genome coverage = "\$0/${genome_size}}' > coverage.txt
    """
}
