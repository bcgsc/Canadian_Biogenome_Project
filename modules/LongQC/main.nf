process LONGQC {
    tag "$meta.id"
    label 'process_medium'

    container '/projects/cbp/scratch/singularity/longqc.sif'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.html'), emit: report
    tuple val(meta), path('*.json'), emit: report_json
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python \\
        ${params.singularity_cache}/LongQC/longQC.py \\
        sampleqc \\
        -o LongQC \\
        --sample_name ${meta.id} \\
        -p 32 \\
        $args \\
        ${reads} \\
        > LongQC.log

    mv LongQC/*.html .
    mv LongQC/*.json .

    rm LongQC/analysis/*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        longqc : \$(python ${params.singularity_cache}/LongQC/longQC.py --version | sed 's/LongQC //g')
    END_VERSIONS
    """
}
