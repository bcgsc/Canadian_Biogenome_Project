process LONGSTITCH {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::longstitch=1.0.4"
//    container "https://depot.galaxyproject.org/singularity/longstitch%3A1.0.3--hdfd78af_0"

    container "docker://quay.io/biocontainers/longstitch:1.0.4--hdfd78af_0"

    input:
    tuple val(meta2), path(reads)
    tuple val(meta), path(assembly)
    val(genome_size)

    output:
    tuple val(meta), path('*.ntLink.scaffolds.fa') , emit: assembly
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def G = genome_size ? "G=$genome_size" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cd  ${assembly} > ${assembly.simpleName}.fa

    #gunzip -c ${reads} > ${reads.simpleName}.fq.gz

    cp $reads ${reads.simpleName}.fq.gz

    longstitch tigmint-ntLink-arks \\
    draft=${assembly.simpleName} \\
    reads=${reads.simpleName} \\
    $args \\
    $G \\
    t=$task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longstitch : \$( longstitch  --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
