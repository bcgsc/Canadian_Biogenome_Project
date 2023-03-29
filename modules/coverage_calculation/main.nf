process COVERAGE_CALCULATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq)

    output:
    path('*.txt'), emit: coverage
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    gzip -cd $fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c | awk '{ print "Genome coverage = "\$0/(${params.hap_gen_size_Gb}*1000000000)}' > coverage.txt
    """
}
