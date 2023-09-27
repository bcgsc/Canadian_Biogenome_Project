process PREPROCESS_MERGED {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::pbtk==3.1.0 bioconda::bamtools=2.5.2"

    input:
    tuple val(meta), path(bam, stageAs: "?/*")

    output:
    tuple val(meta), path('*_filtered.bam'), emit: filtered_bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args_pbmerge = task.ext.args_pbmerge ?: ''
    def args_bamtools_filter = task.ext.args_bamtools_filter ?: ''
    def args_bam2fastq = task.ext.args_bam2fastq ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmerge \\
        -o ${prefix}_merged.bam \\
        $args_pbmerge \\
        */*.bam

    bamtools \\
        filter \\
        -in ${prefix}_merged.bam \\
        $args_bamtools_filter \\
        -out ${prefix}_filtered.bam

    rm *_merged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbbam: \$( pbmerge --version | head -n1 | sed 's/pbmerge //' | sed -E 's/ .+//' )
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
