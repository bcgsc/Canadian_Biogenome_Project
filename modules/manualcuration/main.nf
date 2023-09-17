process RAPIDCURATION_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqtk=1.3 conda-forge::perl=5.32.1"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
//        'biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.tpf")    , emit: split_tpf
    path "versions.yml"               , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    perl ${params.singularity_cache}/rapid-curation/rapid_split.pl -fa $assembly
    mv *.tpf ${prefix}.tpf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(v//g' | sed 's/)//g')
    END_VERSIONS
    """ 
} 
