process BED_PROCESSING {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_sorted.bed"), emit: sorted_bed

    script:
    """
    bedtools bamtobed -i $bam > merged_bed
    sort --parallel=8 --buffer-size=80%  --temporary-directory=$projectDir --output=bed_sorted.bed \$merged_bed
    """
}
