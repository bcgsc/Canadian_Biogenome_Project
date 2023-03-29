process SORT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*_sorted.bed"), emit: sorted_bed

    script:
    """
    sort --parallel=8 --buffer-size=80%  --temporary-directory=$projectDir --output=${bed.simpleName}_sorted.bed $bed
    """
}
