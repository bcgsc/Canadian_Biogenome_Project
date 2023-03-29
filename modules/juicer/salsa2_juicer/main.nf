process SALSA2_JUICER {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salsa2:2.3--py27hee3b9ab_0':
        'quay.io/biocontainers/salsa2:2.3--py27hee3b9ab_0' }"

    input:
    tuple val(meta), path(fasta), path(index), path(agp), path (alignment_iteration_1_bed), path (scaffold_length_iteration_1)

    output:
    tuple val(meta), path("*.hic"), emit: hic
    tuple val(meta), path('*.assembly'), emit: assembly


    script:
    def args = task.ext.args ?: ''
    """
    cut -f 1,2 $index > chromosome_sizes.tsv

    python script/alignments2txt.py -b ${alignment_iteration_1_bed} -a ${agp} -l $scaffold_length_iteration_1 > alignments.txt

    awk '{if (\$2 > \$6) {print \$1"\t"\$6"\t"\$7"\t"\$8"\t"\$5"\t"\$2"\t"\$3"\t"\$4} else {print}}' alignments.txt | sort -k2,2d -k6,6d -T $projectDir --parallel=8 | awk 'NF'  > alignments_sorted.txt

    java -jar ${params.JUICER_JAR} pre $args -o ${meta.id} alignments_sorted.txt salsa_scaffolds.hic chromosome_sizes.tsv
    """
}
