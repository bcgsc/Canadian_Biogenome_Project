process CAT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(file1)
    tuple val(meta2), path (file2)

    output:
    tuple val(meta), path("*_alternate_contigs_full.fa"), emit: alternate_contigs_full
    path  "versions.yml"          , emit: versions

    script:
    """
    cat $file1 $file2 > ${meta.id}_alternate_contigs_full.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | sed 's/cat (GNU coreutils) //g' | sed -n 1p)
    END_VERSIONS
    """
}
