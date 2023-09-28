process QUAST {
    tag "$meta.id" 
    label 'process_medium'

    conda 'bioconda::quast=5.2.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    tuple val(meta), path (consensus)
    path fasta
    path gff
    val use_fasta
    val use_gff
    val genome_size

    output:
    path "*_quast_report.tsv"    , emit: renamed_tsv
    path 'report.tsv'        , emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def est_ref_size = genome_size ? "--est-ref-size $genome_size" : ""
    prefix   = task.ext.prefix ?: 'quast'
    def features  = use_gff ? "--features $gff" : ''
    def reference = use_fasta ? "-r $fasta" : ''
    """
    quast.py \\
        --output-dir $prefix \\
        $reference \\
        $features \\
        --threads $task.cpus \\
        $args \\
        $est_ref_size \\
        ${consensus.join(' ')}

    mv ${prefix}/report.tsv report.tsv
    cp report.tsv ${consensus.baseName}_quast_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
