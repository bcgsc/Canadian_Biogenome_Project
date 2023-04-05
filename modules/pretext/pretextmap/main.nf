
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextmap=0.1.9 bioconda::samtools=1.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c6242a6c1a522137de7a9e9ff90779ede11cf5c5-0':
        'quay.io/biocontainers/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c6242a6c1a522137de7a9e9ff90779ede11cf5c5-0' }"

    input:
    tuple val(meta), path(chrom_sizes)
    tuple val(meta), path(alignments_sorted_txt)

    output:
    tuple val(meta), path("*.pretext"), emit: pretext

    script:
    """
    (awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"\$1"\t"\$2} END {print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' $chrom_sizes; awk '{print ".\t"\$2"\t"\$3"\t"\$6"\t"\$7"\t.\t."}' $alignments_sorted_txt) | PretextMap -o ${meta.id}.pretext
    """
}
