process MATLOCK {
    tag "$meta.id"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/matlock%3A20181227--hd02f7ee_0"

    input:
    tuple val(meta), path(hic_bam)

    output:
    tuple val(meta), path("*sorted.links.txt"), emit: links_txt

    script:
    def args = task.ext.args ?: ''
    """
    matlock bam2 juicer $hic_bam out.links.txt  # this step sometimes crashes on memory
    sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
    """
}
