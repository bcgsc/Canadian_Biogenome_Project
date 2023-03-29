process LONGQC {
    tag "$meta.id"
    label 'process_medium'

    container 'longqc.sif'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.html'), emit: report
    tuple val(meta), path('*.json'), emit: report_json

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python \\
	$args \\
	longQC.py \\
	sampleqc \\
	-x pb-sequel \\
	-o LongQC \\
	--sample_name ${meta.id} \\
	-p 32 \\
	${reads} \\
        > LongQC.log

    mv LongQC/*.html .
    mv LongQC/*.json .

    rm LongQC/analysis/*.fastq
    """
}
