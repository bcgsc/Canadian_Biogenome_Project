process GFA_TO_FA {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(GFA_assembly)

    output:
    tuple val(meta), path('*.fa'), emit: fa_assembly
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    awk '/^S/{print ">"\$2;print \$3}' $GFA_assembly > ${GFA_assembly.baseName}.fa
    """
}
