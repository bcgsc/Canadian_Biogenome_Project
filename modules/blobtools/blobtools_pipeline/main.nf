process BLOBTOOLS_PIPELINE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(config)
    tuple val(meta), path(fa_gz)

    output:
    tuple val(meta), path('*.blastn.nt.out'), emit: blast_out
    tuple val(meta), path('*.diamond.reference_proteomes.out'), emit: diamond_proteome_out
    tuple val(meta), path('*.diamond.busco_genes.out'), emit: diamond_busco_out
    tuple val(meta), path('assembly_minimap.bam'), emit: assembly_minimap_bam
    tuple val(meta), path('hic_minimap.bam'), emit:hic_minimap_bam
    tuple val(meta), path('lineage1_full_table.tsv.gz'), emit: lineage1_full_table_tsv
    tuple val(meta), path('lineage2_full_table.tsv.gz'), emit: lineage2_full_table_tsv, optional: true
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    # Load the Conda environment containing snakemake
    source activate blobtools

    # Run the snakemake pipeline
    #Directory is where the log files are going to be created
    snakemake -p \
          -j 60 \
          --directory . \
          --configfile $config \
          --latency-wait 60 \
          --stats blobtoolkit.stats \
          -s ${params.blobtoolkit_path}/insdc-pipeline/blobtoolkit.smk

    cp blastn/*.blastn.nt.out .
    cp diamond/*.diamond.reference_proteomes.out .
    cp diamond_blastp/*.diamond.busco_genes.out .
    cp minimap/*.bam .
    mv *.${meta.id}.bam assembly_minimap.bam
    mv ${meta.id}.*.bam hic_minimap.bam
    cp busco/${meta.id}.busco.*/full_table.tsv.gz lineage1_full_table.tsv.gz
    #cp busco/${meta.id}.busco.${params.lineage2}/full_table.tsv.gz lineage2_full_table.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snakemake : \$(snakemake --version)
        minimap2: \$(minimap2 --version 2>&1)
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
        busco: \$(singularity run -B /projects ${params.singularity_cache}/busco5.sif busco --version 2>&1 | sed 's/^BUSCO //' ) 
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//' | head -1)
    END_VERSIONS
    """
}
