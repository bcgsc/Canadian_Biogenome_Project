process THREED_DNA {
    tag "$meta.id"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/3d-dna%3A201008--hdfd78af_0"

    input:
    tuple val(meta), path(links_txt)
    tuple val(meta2), path(assembly)

    output:
    tuple val(meta), path("*sorted.links.txt"), emit: links_txt

    script:
    def args = task.ext.args ?: ''
    """
    ./run-assembly-visualizer.sh -p false $assembly $links_txt
    """
}
