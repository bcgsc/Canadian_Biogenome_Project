process HICPRO_SAMPLESHEET {
    tag "$meta.id"
    label 'process_LOW'

    input:
    tuple val(meta), path(hic)

    output:
    tuple val(meta), path("*.csv"), emit: hicpro_csv

    script:
    """
    echo "sample,fastq_1,fastq_2\n${params.id},${params.hic_read1},${params.hic_read2}" > hicpro_samplesheet.csv
    """
}
