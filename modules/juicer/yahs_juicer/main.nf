process YAHS_JUICER {
    tag "$meta.id"
    label 'process_medium'

    container "https://depot.galaxyproject.org/singularity/yahs%3A1.2a.2--h7132678_0"

    input:
    tuple val(meta), path(agp), path (bin)
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("*_scaffolds_final.chrom.sizes"), emit: chrom_sizes
    tuple val(meta), path('*_alignments_sorted.txt'), emit: alignments_sorted_txt
    tuple val(meta), path('*.assembly'), emit: assembly
    tuple val(meta), path('*.txt'), emit: juicer_txt
    tuple val(meta), path('*.log'), emit: juicer_log

    script:
    def args = task.ext.args ?: ''
    """
    juicer pre \\
	 -a -o ${meta.id} \\
	$bin \\
	$agp \\
	$index

     juicer pre \\
	$bin \\
	$agp \\
	$index 2>tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -S32G | awk 'NF' > ${meta.id}_alignments_sorted.txt

    cat tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- > ${meta.id}_scaffolds_final.chrom.sizes
    """
}
