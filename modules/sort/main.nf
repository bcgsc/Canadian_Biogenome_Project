process SORT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*_sorted.bed"), emit: sorted_bed
    path  "versions.yml"          , emit: versions

    script:
    """
    sort --parallel=8 --buffer-size=80%  --temporary-directory=$projectDir --output=${bed.simpleName}_sorted.bed $bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed 's/sort (GNU coreutils) //g' | sed -n 1p)
    END_VERSIONS
    """
}
