process CPG {
        input :
        tuple val(meta), path (aligned_bam)
        tuple val(meta2), path (bam_index)
        tuple val(meta3), path (ref)

        output :
        tuple val(meta), path ('*.bed'), emit : cbp_bed
        tuple val(meta), path ('*.bw') , emit : cbp_bigwig
        tuple val(meta), path ('*.log'), emit : log
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
	"""
	/projects/cbp/scratch/methylation_SC/pb-CpG-tools/pb-CpG-tools-v2.1.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
	  --bam $aligned_bam \
	  --output-prefix ${prefix} \
	  --model /projects/cbp/scratch/methylation_SC/pb-CpG-tools/pb-CpG-tools-v2.1.0-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
	  --threads 16 \
	  --ref $ref \
	  --pileup-mode model


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbp: \$(/projects/cbp/scratch/methylation_SC/pb-CpG-tools/pb-CpG-tools-v2.1.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores --version | sed 's/cpg v//g')
    END_VERSIONS
	"""
}

