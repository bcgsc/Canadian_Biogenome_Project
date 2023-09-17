process JUICER {
    tag "$meta.id"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/java-jdk%3A8.0.92--1"

    input:
    tuple val(meta), path(chrom_sizes)
    tuple val(meta), path(alignments_sorted_txt)

    output:
    tuple val(meta), path("*.hic"), emit: hic
    path  "versions.yml"          , emit: versions

    script:
    """
    java -Xmx400G -jar ${params.JUICER_JAR} pre \\
        $alignments_sorted_txt \\
        ${meta.id}.hic \\
        $chrom_sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openjdk: \$(echo \$(java -version 2>&1) | grep version | sed 's/\"//g' | cut -f3 -d ' ')
        juicer : \$(echo \$(java -jar ${params.JUICER_JAR} --version | sed -n 2p | sed 's/Juicer Tools Version //g'))
    END_VERSIONS
    """
}
