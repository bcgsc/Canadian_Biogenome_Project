process FCS_FCSADAPTOR {
    stageInMode 'copy'
    tag "$meta.id"
    label 'process_low'

    // Use local Singularity image (if provided)
    container "/projects/cbp/scratch/singularity/fcs-adaptor-with-bash.sif"

    // Exit if running this module with -profile conda / mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "FCS_FCSADAPTOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.cleaned_sequences.fa.gz"), emit: cleaned_assembly
    tuple val(meta), path("*.fcs_adaptor_report.txt") , emit: adaptor_report
    tuple val(meta), path("*.fcs_adaptor.log")        , emit: log
    tuple val(meta), path("*.pipeline_args.yaml")     , emit: pipeline_args
    tuple val(meta), path("*.skipped_trims.jsonl")    , emit: skipped_trims
    path "versions.yml"                               ,emit: versions

    script:
    def args = task.ext.args ?: '--euk'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FCSADAPTOR_VERSION = '0.5.5'
    """
    mkdir -p output

    /app/fcs/bin/av_screen_x \\
        ./${assembly} \\
        -o ./output \\
        $args

    gzip -cf output/cleaned_sequences/* > "${assembly.baseName}.cleaned_sequences.fa.gz"
    cp output/fcs_adaptor_report.txt     "${assembly.baseName}.fcs_adaptor_report.txt"
    cp output/fcs_adaptor.log            "${assembly.baseName}.fcs_adaptor.log"
    cp output/pipeline_args.yaml         "${assembly.baseName}.pipeline_args.yaml"
    cp output/skipped_trims.jsonl        "${assembly.baseName}.skipped_trims.jsonl"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FCS-adaptor: $FCSADAPTOR_VERSION
    END_VERSIONS
    """
}
