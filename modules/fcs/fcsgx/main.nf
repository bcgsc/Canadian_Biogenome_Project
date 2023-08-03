process FCS_FCSGX {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-gx.sif' :
        'ncbi/fcs-gx:0.4.0' }"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "FCS_FCSGX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.fcs_gx_report.txt"), emit: fcs_gx_report
    tuple val(meta), path("*.taxonomy.rpt")     , emit: taxonomy_report
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSGX_VERSION = '0.4.0'

    """
    python3 /app/bin/run_gx \\
        --fasta $assembly \\
        --out-dir . \\
        --gx-db ${params.fcs_gx_database} \\
        --tax-id ${params.taxon_taxid} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed -e "s/Python //g")
        FCS-GX: $FCSGX_VERSION
    END_VERSIONS
    """
}
