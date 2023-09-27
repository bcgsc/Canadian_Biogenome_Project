process BLOBTOOLS_CONFIG_1LINEAGE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(pacbio_fastq)

    output:
    tuple val(meta), path('config.yaml'), emit: config
 
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
//The assembly file must be compressed (fa.gz) bgzip -c 
    """
    echo "assembly:
      accession: ${meta.id}
      file: ${params.outdir}/QC/blobtools/${assembly}
      level: scaffold
      prefix: ${meta.id}
    busco:
      download_dir: ${params.busco_lineages_path}
      lineages:
        - ${params.lineage}
      basal_lineages:
        - ${params.lineage}
    reads:
      paired:
        - prefix: ${params.Illumina_prefix}
          platform: ILLUMINA
          file: ${params.hic_read1};${params.hic_read2}
      single:
        - prefix: ${meta.id}
          platform: PACBIO_SMRT
          file: ${params.outdir}/preprocessing/bam2fastx/${pacbio_fastq}
    revision: 0
    settings:
      blast_chunk: 100000
      blast_max_chunks: 10
      blast_overlap: 0
      blast_min_length: 1000
      taxdump: ${params.Blobtoolkit_db}/taxdump 
      tmp: /tmp
    similarity:
      defaults:
        evalue: 1.0e-10
        import_evalue: 1.0e-25
        max_target_seqs: 10
        taxrule: bestdistorder
      diamond_blastx:
        name: reference_proteomes
        path: ${params.Blobtoolkit_db}/uniprot
      diamond_blastp:
        name: reference_proteomes
        path: ${params.Blobtoolkit_db}/uniprot
        import_max_target_seqs: 100000
      blastn:
        name: nt
        path: ${params.Blobtoolkit_db}/nt
    taxon:
      name: ${params.taxon_name}
      taxid: '${params.taxon_taxid}'
    version: 1" > config.yaml    
    """
}
