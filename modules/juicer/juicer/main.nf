process JUICER {
    tag "$meta.id"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/java-jdk%3A8.0.92--1"

    input:
    tuple val(meta), path(chrom_sizes)
    tuple val(meta), path(alignments_sorted_txt)

    output:
    tuple val(meta), path("*.hic"), emit: hic

    script:
    """
    java -Xmx400G -jar ${params.JUICER_JAR} pre \\
        $alignments_sorted_txt \\
        ${meta.id}.hic \\
        $chrom_sizes
    """
}
