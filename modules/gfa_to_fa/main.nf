process GFA_TO_FA {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(GFA_assembly)

    output:
    tuple val(meta), path('*.fa'), emit: fa_assembly
    path  "versions.yml"          , emit: versions
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    awk '/^S/{print ">"\$2;print \$3}' $GFA_assembly > ${GFA_assembly.baseName}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed 's/GNU Awk //g'| sed -n 1p)
    END_VERSIONS
    """
}
