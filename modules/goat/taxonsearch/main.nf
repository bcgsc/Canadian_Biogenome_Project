process GOAT_TAXONSEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::goat==0.2.0 libgcc-ng"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goat:0.2.0--h92d785c_0':
        'biocontainers/goat:0.2.0--h92d785c_0' }"

    input:
    tuple val(meta), val(taxon), path(taxa_file)

    output:
    tuple val(meta), path("*.tsv"), emit: taxonsearch
    path "versions.yml"           , emit: versions
    env(ploidy)	 		  , emit: ploidy
    env(haploid_number)		  , emit: haploid_number
    env(scientific_name)	  , emit: scientific_name
    env(genome_size)		  , emit: genome_size

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input = taxa_file ? "-f ${taxa_file}" : "-t \"${taxon}\""
    if (!taxon && !taxa_file) error "No input. Valid input: single taxon identifier or a .txt file with identifiers"
    if (taxon && taxa_file ) error "Only one input is required: a single taxon identifier or a .txt file with identifiers"
    """
    goat-cli taxon search \\
        $args \\
        $input > ${prefix}.tsv

    sed 's/\t/,/g' ${prefix}.tsv > ${prefix}.csv

    ploidy=\$(awk -v get='^(ploidy)' 'BEGIN{FS=OFS=","}FNR==1{for(i=1;i<=NF;i++)if(\$i~get)cols[++c]=i}{for(i=1; i<=c; i++)printf "%s%s", \$(cols[i]), (i<c ? OFS : ORS)}' ${prefix}.csv | awk 'FNR == 2')
    haploid_number=\$(awk -v get='^(haploid_number)' 'BEGIN{FS=OFS=","}FNR==1{for(i=1;i<=NF;i++)if(\$i~get)cols[++c]=i}{for(i=1; i<=c; i++)printf "%s%s", \$(cols[i]), (i<c ? OFS : ORS)}' ${prefix}.csv | awk 'FNR == 2')
    scientific_name=\$(awk -v get='^(scientific_name)' 'BEGIN{FS=OFS=","}FNR==1{for(i=1;i<=NF;i++)if(\$i~get)cols[++c]=i}{for(i=1; i<=c; i++)printf "%s%s", \$(cols[i]), (i<c ? OFS : ORS)}' ${prefix}.csv | awk 'FNR == 2')
    genome_size=\$(awk -v get='^(genome_size)' 'BEGIN{FS=OFS=","}FNR==1{for(i=1;i<=NF;i++)if(\$i~get)cols[++c]=i}{for(i=1; i<=c; i++)printf "%s%s", \$(cols[i]), (i<c ? OFS : ORS)}' ${prefix}.csv | awk 'FNR == 2' | awk -F, '{print \$1}') 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goat: \$(goat-cli --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
