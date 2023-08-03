process FCS_FCSGX_CLEAN {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(fcsgx_report)

    output:
    tuple val(meta), path("*.cleaned.fasta.gz")     , emit: cleaned_fasta
    tuple val(meta), path("*.contam.fasta.gz")     , emit: cantam_fasta
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.4.0'
    """
    cp $fcsgx_report local_report.txt

    zcat $assembly | python3 ${params.singularity_cache}/fcs.py \\
        --image=${params.singularity_cache}/ftp.ncbi.nlm.nih.gov-genomes-TOOLS-FCS-releases-0.4.0-fcs-gx.sif \\
        clean genome \\
        --action-report=local_report.txt \\
        --output=cleaned.fasta \\
        --contam-fasta-out=contam.fasta

    gzip -c cleaned.fasta > ${assembly.baseName}.cleaned.fasta.gz
    gzip -c contam.fasta > ${assembly.baseName}.contam.fasta.gz
    rm cleaned.fasta
    rm contam.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """
}
