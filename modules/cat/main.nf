process CAT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(file1)
    tuple val(meta2), path (file2)

    output:
    tuple val(meta), path("*_alternate_contigs_full.fa"), emit: alternate_contigs_full

    script:
    """
    cat $file1 $file2 > ${meta.id}_alternate_contigs_full.fa
    """
}
