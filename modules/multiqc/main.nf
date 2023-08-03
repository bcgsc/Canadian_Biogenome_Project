process MULTIQC {
    label 'process_single'

    conda (params.enable_conda ? 'bioconda::multiqc=1.13' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(coverage_estimation)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    """
    echo ${coverage_estimation}
    coverage_estimation=\$(cat $coverage_estimation)
    echo \$coverage_estimation

    multiqc \\
        $args \\
        $config \\
        $extra_config \\
        -b "Specie common name : ${params.id}"  \\
        -b "Specie taxonomic ID : ${params.taxon_taxid}" \\
        -b "Specie scientific name : ${params.taxon_name}" \\
        -b "Specie ploidy : ${params.ploidy}" \\
        -b "Specie genome size : ${params.hap_gen_size_Gb}" \\
        -b "Specie number of chromosomes : ${params.chrom_num}" \\
        -b "\${coverage_estimation}" \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
