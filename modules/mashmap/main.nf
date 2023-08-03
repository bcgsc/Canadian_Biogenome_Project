process MASHMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mashmap=3.0.4 bioconda::perl-bioperl=1.7.2 conda-forge::gsl=2.7 mkl=2023.1.0 conda-forge::gnuplot=5.4.5"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/mashmap%3A3.0.4--h97b747e_0' :
//        'biocontainers/mashmap:3.0.4--h97b747e_0' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.png")    , emit: mashmap_png
    tuple val(meta), path("*.txt"), emit: mashmap_txt
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    mashmap \\
        -r $reference \\
        -q $assembly \\
        $args \\
        -t $task.cpus \\
        -o ${prefix}_${reference.simpleName}.out

    perl ${params.singularity_cache}/MashMap/scripts/generateDotPlot png large ${prefix}_${reference.simpleName}.out

    awk '{print \$1"\t"\$6"\t"\$5}' ${prefix}_${reference.simpleName}.out | sort | uniq > mashmap_correspondance.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mashmap: \$(echo \$(mashmap --version 2>&1))
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(v//g' | sed 's/)//g')
    END_VERSIONS
    """ 
} 
