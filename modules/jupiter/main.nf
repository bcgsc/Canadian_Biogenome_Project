process JUPITER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.png")     , emit: plot
    tuple val(meta), path("*.txt")     , emit: txt
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    source \$(conda info --json | awk '/conda_prefix/ { gsub(/"|,/, "", \$2); print \$2 }')/bin/activate  /home/scorreard/miniconda3/envs/circos

    gzip -cd $reference > jupiter_reference.fa

    ${params.singularity_cache}/Jupiter/./jupiter \\
        t=$task.cpus \\
        name=${prefix}_${reference.simpleName} \\
        ref=jupiter_reference.fa \\
        fa=$assembly \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jupiter: \$( ${params.singularity_cache}/Jupiter/./jupiter --version | head -n1)
        circos : \$( circos --version | sed 's/circos | v //g'| sed 's/ | 15 Jun 2019 | Perl 5.032001//g') 
    END_VERSIONS
    """
}
