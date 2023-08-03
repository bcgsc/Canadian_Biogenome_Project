process OVERVIEW_GENERATION_SAMPLE {
	tag "$meta.id"
    label 'process_low'

        //container = 'https://depot.galaxyproject.org/singularity/r-stringr%3A1.1.0--r3.3.1_0'
	container = 'https://depot.galaxyproject.org/singularity/r-rjson%3A0.2.15--r3.3.2_0'

    input:
    tuple val(meta), path(longqc) //LonQC
    tuple val(meta1), path(kraken_pacbio) //kraken_pacbio
    tuple val(meta2), path(kraken_hic) //kraken_hic
    path(quast_contig) //quast_contig
    path(quast_contig_purged) //quast_contig_purged
    path(quast_scaffold) //quast_scaffold
    tuple val(meta6), path(busco_lin1) //busco_lineage1
    tuple val(meta7), path(busco_lin2) //busco_lineage2
    tuple val(meta8), path(busco_lin3) //busco_lineage3
    tuple val(meta9), path(busco_lin4) //busco_lineage4
    tuple val(meta10), path(chrom_size)
    val (ploidy)
    val(haploid_number)
    val(scientific_name)
    val(genome_size)

    output :
    tuple val(meta), path('*.tsv'), emit: overview

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    """
    Rscript ${params.modules_path}/overview_generation/sample/overview_generation_sample.R \\
	$params.id \\
        $params.pipeline_version \\
	$params.outdir \\
	$params.pacbio_input_type \\
	$params.bam_cell1 \\
	$params.bam_cell2 \\
	$params.bam_cell3 \\
	$params.bam_cell4 \\
	$params.ont_fastq_1 \\
	$params.hic_read1 \\
	$params.hic_read2 \\
	$params.illumina_SR_read1 \\
	$params.illumina_SR_read2 \\
	$params.pacbio_rq \\
	$params.assembly_method \\
	$params.assembly_secondary_mode \\
	$params.polishing_method \\
        $params.purging_method \\
	$params.scaffolding_method \\
        $params.manual_curation \\
        $params.mitohifi \\
	$params.taxon_taxid \\
	$ploidy \\
	$genome_size \\
	$haploid_number \\
	$params.lineage \\
	$params.lineage2 \\
	$params.lineage3 \\
	$params.lineage4 \\
	$longqc \\
	$kraken_pacbio \\
	$kraken_hic \\
	$quast_contig \\
	$quast_contig_purged \\
	$quast_scaffold \\
	$busco_lin1 \\
	$busco_lin2 \\
	$busco_lin3 \\
	$busco_lin4 \\
        $chrom_size \\
	$scientific_name
    """
}
