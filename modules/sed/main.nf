process SED_NONE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(file1)

    output:
    tuple val(meta), path("*_sed.fa"), emit: assembly

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed 's/None-None//' $file1 > ${prefix}_sed.fa
    """
}
