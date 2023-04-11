#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
CBP pipeline - Solenne Correard - Jones lab
=============================================
Specie id			: ${params.id}
Taxon				: ${params.taxon_name}
Taxon number			: ${params.taxon_taxid}
PacBio input type		: ${params.pacbio_input_type}
PacBio reads cell 1		: ${params.bam_cell1}
PacBio reads cell 2		: ${params.bam_cell2}
PacBio reads cell 3     	: ${params.bam_cell3}
PacBio reads cell 4             : ${params.bam_cell4}
PacBio read quality threashold  : ${params.pacbio_rq}
ONT reads			: ${params.ont_fastq_1}
Hi-C reads F			: ${params.hic_read1}
Hi-C reads R     		: ${params.hic_read2}
Output path			: ${params.outdir}
Pipeline version		: ${params.pipeline_version}
Assembly method			: ${params.assembly_method}
Assembly mode			: ${params.assembly_secondary_mode}
Polishing method		: ${params.polishing_method}
Purging method			: ${params.purging_method}
Scaffolding method		: ${params.scaffolding_method}
Manual curation			: ${params.manual_curation}
Mitochondrial assembly		: ${params.mitohifi}
"""

//Pre-processing
include { CCS as CCS_PACBIO_CELL1; CCS as CCS_PACBIO_CELL2; CCS as CCS_PACBIO_CELL3; CCS as CCS_PACBIO_CELL4 } from './modules/pacbio/ccs/main.nf'
include { BAMTOOLS_FILTER as BAMTOOLS_FILTER_PACBIO_CELL1; BAMTOOLS_FILTER as BAMTOOLS_FILTER_PACBIO_CELL2; BAMTOOLS_FILTER as BAMTOOLS_FILTER_PACBIO_CELL3; BAMTOOLS_FILTER as BAMTOOLS_FILTER_PACBIO_CELL4 } from './modules/bamtools_filter/main.nf'
include { PBINDEX as PBINDEX_FILTERED_PACBIO_CELL1; PBINDEX as PBINDEX_FILTERED_PACBIO_CELL2; PBINDEX as PBINDEX_FILTERED_PACBIO_CELL3; PBINDEX as PBINDEX_FILTERED_PACBIO_CELL4 } from './modules/pacbio/pbbam/pbindex/main.nf'
include { BAM2FASTX } from './modules/pacbio/bam2fastx/main.nf'
include { TWOBAM2FASTX } from './modules/pacbio/bam2fastx/2bam2fastx/main.nf'
include { THREEBAM2FASTX } from './modules/pacbio/bam2fastx/3bam2fastx/main.nf'
include { FOURBAM2FASTX } from './modules/pacbio/bam2fastx/4bam2fastx/main.nf'

include { CUTADAPT }  from './modules/cutadapt/main.nf'

//QC Input data
include { LONGQC } from './modules/LongQC/main.nf'
include { MERYL_COUNT } from './modules/meryl/count/main.nf'
include { MERYL_UNIONSUM } from './modules/meryl/unionsum/main.nf'
include { MERYL_HISTOGRAM } from './modules/meryl/histogram/main.nf'
include { GENOMESCOPE2 } from './modules/genomescope2/main.nf'
include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_PACBIO_BAM; KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_HIC_READS; KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_SR_READS; KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ONT_READS } from './modules/kraken2/main.nf'
include { COVERAGE_CALCULATION } from './modules/coverage_calculation/main.nf'

//Assembly
//HifiASM
include { HIFIASM } from './modules/hifiasm/main.nf'
include { GFA_TO_FA; GFA_TO_FA as GFA_TO_FA2 } from './modules/gfa_to_fa/main.nf'

//Canu
include { CANU } from './modules/canu/main.nf'

//Flye
include { FLYE } from './modules/flye/main.nf'
include { FLYE_PACBIO_ONT } from './modules/flye/flye_pacbio_ont/main.nf'
include { MINIMAP2_ALIGN as MINIMAP_ALIGN_FLYE } from './modules/minimap2/align/main.nf'
include { RACON } from './modules/racon/main.nf'
include { LONGSTITCH } from './modules/longstitch/main.nf'

//Verkko
include { VERKKO } from './modules/verkko/main.nf'

//Polishing
//Pilon
include { PILON } from './modules/pilon/main.nf'
include { BWAMEM2_INDEX } from './modules/bwamem2/index/main.nf'
include { BWAMEM2_MEM } from './modules/bwamem2/mem/main.nf'
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

//PurgeDups
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_CONTIG; MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_SELF; MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_CONTIG_ALT; MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_SELF_ALT } from './modules/minimap2/align/main.nf'
include { PURGEDUPS_SPLITFA; PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_ALT } from './modules/purgedups/splitfa/main.nf'
include { PURGEDUPS_PBCSTAT; PURGEDUPS_PBCSTAT as PURGEDUPS_PBCSTAT_ALT } from './modules/purgedups/pbcstat/main.nf'
include { PURGEDUPS_CALCUTS; PURGEDUPS_CALCUTS as PURGEDUPS_CALCUTS_ALT } from './modules/purgedups/calcuts/main.nf'
include { PURGEDUPS_PURGEDUPS; PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_ALT } from './modules/purgedups/purgedups/main.nf'
include { PURGEDUPS_GETSEQS; PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_ALT } from './modules/purgedups/getseqs/main.nf'

//HIC scaffolding-SALSA2
include { PREPARE_GENOME } from './modules/nfcore_hic/subworkflows/local/prepare_genome.nf'
include { FASTQC } from './modules/nfcore_hic/modules/nf-core/fastqc/main.nf'
include { HICPRO } from './modules/nfcore_hic/subworkflows/local/hicpro.nf'
include { COOLER } from './modules/nfcore_hic/subworkflows/local/cooler.nf'
include { HIC_PLOT_DIST_VS_COUNTS } from './modules/nfcore_hic/modules/local/hicexplorer/hicPlotDistVsCounts.nf'
include { COMPARTMENTS } from './modules/nfcore_hic/subworkflows/local/compartments.nf'
include { TADS } from './modules/nfcore_hic/subworkflows/local/tads.nf'

include { MINIMAP2_ALIGN as MINIMAP_ALIGN_HIC_F_GETSEQS; MINIMAP2_ALIGN as MINIMAP_ALIGN_HIC_R_GETSEQS } from './modules/minimap2/align/main.nf'
include { BEDTOOLS_BAMTOBED as BEDTOOLS_BAMTOBED_HIC_F_GETSEQS; BEDTOOLS_BAMTOBED as BEDTOOLS_BAMTOBED_HIC_R_GETSEQS;} from './modules/bedtools/bamtobed/main.nf'
include { BED_PROCESSING } from './modules/bed_processing/main.nf'
include { SALSA2 } from './modules/salsa2/main.nf'
include { SALSA2_JUICER } from './modules/juicer/salsa2_juicer/main.nf'

//HIC scaffolding-YAHS
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX1; SAMTOOLS_FAIDX as SAMTOOLS_FAIDX2; SAMTOOLS_FAIDX as SAMTOOLS_FAIDX1_ALT; SAMTOOLS_FAIDX as SAMTOOLS_FAIDX2_ALT } from './modules/samtools/faidx/main.nf'
include { CHROMAP_INDEX; CHROMAP_INDEX as CHROMAP_INDEX_ALT } from './modules/chromap/index/main.nf'
include { CHROMAP_CHROMAP; CHROMAP_CHROMAP as CHROMAP_CHROMAP_ALT } from './modules/chromap/chromap/main.nf'
include { YAHS; YAHS as YAHS_ALT } from './modules/yahs/main.nf'

//Assembly QC
include { YAHS_JUICER } from './modules/juicer/yahs_juicer/main.nf'
include { JUICER } from './modules/juicer/juicer/main.nf'
include { PRETEXTMAP } from './modules/pretext/pretextmap/main.nf'
include { PRETEXTSNAPSHOT } from './modules/pretext/pretextsnapshot/main.nf'

include { CAT } from './modules/cat/main.nf'
include { BUSCO ; BUSCO as BUSCO_lin2; BUSCO as BUSCO_lin3; BUSCO as BUSCO_lin4; BUSCO as BUSCO_ALT; BUSCO as BUSCO_lin2ALT; BUSCO as BUSCO_lin3ALT; BUSCO as BUSCO_lin4ALT } from './modules/busco/main.nf'
include { MERQURY as MERQURY1; MERQURY as MERQURY2; MERQURY as MERQURY3 } from './modules/merqury/main.nf'
include { MERQURY_DOUBLE as MERQURY1_DOUBLE; MERQURY_DOUBLE as MERQURY2_DOUBLE; MERQURY_DOUBLE as MERQURY3_DOUBLE } from './modules/merqury/merqury_double/main.nf'
include { QUAST as QUAST1; QUAST as QUAST2; QUAST as QUAST3; QUAST as QUAST_PILON } from './modules/quast/main.nf'
include { QUAST_DOUBLE as QUAST1_DOUBLE; QUAST_DOUBLE as QUAST2_DOUBLE; QUAST_DOUBLE as QUAST3_DOUBLE } from './modules/quast/quast_double/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

include { GZIP } from './modules/gzip/main.nf'
include { BLOBTOOLS_CONFIG } from './modules/blobtools/blobtools_config/main.nf'
include { BLOBTOOLS_PIPELINE } from './modules/blobtools/blobtools_pipeline/main.nf'
include { BLOBTOOLS_CREATE } from './modules/blobtools/blobtools_create/main.nf'
include { BLOBTOOLS_ADD } from './modules/blobtools/blobtools_add/main.nf'
include { BLOBTOOLS_VIEW_SNAIL } from './modules/blobtools/blobtools_view_snail/main.nf'
include { BLOBTOOLS_VIEW_BLOB } from './modules/blobtools/blobtools_view_blob/main.nf'
include { BLOBTOOLS_VIEW_CUMULATIVE } from './modules/blobtools/blobtools_view_cumulative/main.nf'

//Mitochondrial assembly
include { FASTQGZ_TO_FASTA } from './modules/fastqgz_to_fasta/main.nf'
include { FIND_MITO_REFERENCE } from './modules/mitohifi/findmitoreference/main.nf'
include { MITOHIFI } from './modules/mitohifi/mitohifi/main.nf'

workflow {

//////////////////////////////////////////////////  INPUT   //////////////////////////////////

//PacBio data
	input_pacbio_cell1 = [
                [ id:params.id, single_end: true], // meta map
                [ file(params.bam_cell1, checkIfExists: true) ]
           ]
	if( params.bam_cell2 ){
		input_pacbio_cell2 = [
                	[ id:'pacbio_cell2', single_end: true], // meta map
                	[ file(params.bam_cell2, checkIfExists: true) ]
           	]
	}
	if( params.bam_cell3 ){
		input_pacbio_cell3 = [
                	[ id:'pacbio_cell3', single_end: true], // meta map
                	[ file(params.bam_cell3, checkIfExists: true) ]
           	]
	}
        if( params.bam_cell4 ){
                input_pacbio_cell4 = [
                        [ id:'pacbio_cell4', single_end: true], // meta map
                        [ file(params.bam_cell4, checkIfExists: true) ]
                ]
        }


//ONT data
	if (params.ont_fastq_1) {
		input_ont_fastq_1 = [
                        [ id:'input_ont_fastq_1', single_end: true], // meta map
                        [ file(params.ont_fastq_1, checkIfExists: true) ]
                ]
	}

//Hi-C data (Only if hic defined)
        if (( params.hic_read1 ) && ( params.hic_read2 )) {
		input_hic_R1_R2 = [
        	[ id:'hic_R1R2', single_end: false], // meta map
        	[
                	file(params.hic_read1, checkIfExists: true),
                	file(params.hic_read2, checkIfExists: true) 
        	]
    	]
	}

//Illumina SR data (Only if illumina_SR defined)
        if (( params.illumina_SR_read1 ) && ( params.illumina_SR_read2 )) {
                input_illumina_SR_R1_R2 = [
                [ id:'input_illumina_SR_R1_R2', single_end: false], // meta map
                [
                        file(params.illumina_SR_read1, checkIfExists: true),
                        file(params.illumina_SR_read2, checkIfExists: true)
                ]
        ]
        }

//////////////////////////////////////////////////  DUMMY FILES   //////////////////////////////////
        quast_fasta             = file('fasta_dummy')
        quast_gff                = file('gff_dummy')
	dummy_file		= file('dummy')
	input_dummy = [
                [ id:'dummy', single_end: true], // meta map
                [ file('dummy')]
	]

//////////////////////////////////////////////////  WORKFLOW   //////////////////////////////////

	if (params.pacbio_input_type == 'subreads') {
		CCS_PACBIO_CELL1(input_pacbio_cell1)
		if( params.bam_cell2 ) {
			CCS_PACBIO_CELL2(input_pacbio_cell2)
		}
                if( params.bam_cell3 ) {
                        CCS_PACBIO_CELL3(input_pacbio_cell3)
                }
                if( params.bam_cell4 ) {
                        CCS_PACBIO_CELL4(input_pacbio_cell4)
                }

                //Pre-processing
                BAMTOOLS_FILTER_PACBIO_CELL1 (CCS_PACBIO_CELL1.out.bam)
                PBINDEX_FILTERED_PACBIO_CELL1 (BAMTOOLS_FILTER_PACBIO_CELL1.out.filtered_bam)

                if( params.bam_cell2 ) {
                        BAMTOOLS_FILTER_PACBIO_CELL2 (CCS_PACBIO_CELL2.out.bam)
                        PBINDEX_FILTERED_PACBIO_CELL2 (BAMTOOLS_FILTER_PACBIO_CELL2.out.filtered_bam)
                }
                if( params.bam_cell3 ){
                        BAMTOOLS_FILTER_PACBIO_CELL3 (CCS_PACBIO_CELL3.out.bam)
                        PBINDEX_FILTERED_PACBIO_CELL3 (BAMTOOLS_FILTER_PACBIO_CELL3.out.filtered_bam)
                }
		if( params.bam_cell4 ){
                        BAMTOOLS_FILTER_PACBIO_CELL4 (CCS_PACBIO_CELL4.out.bam)
                        PBINDEX_FILTERED_PACBIO_CELL4 (BAMTOOLS_FILTER_PACBIO_CELL4.out.filtered_bam)
                }
	} else {
		//Pre-processing
		BAMTOOLS_FILTER_PACBIO_CELL1 (input_pacbio_cell1)
		PBINDEX_FILTERED_PACBIO_CELL1 (BAMTOOLS_FILTER_PACBIO_CELL1.out.filtered_bam)

		if( params.bam_cell2 ) {
			BAMTOOLS_FILTER_PACBIO_CELL2 (input_pacbio_cell2)
        		PBINDEX_FILTERED_PACBIO_CELL2 (BAMTOOLS_FILTER_PACBIO_CELL2.out.filtered_bam)
		}
		if( params.bam_cell3 ){
        	        BAMTOOLS_FILTER_PACBIO_CELL3 (input_pacbio_cell3)
        	        PBINDEX_FILTERED_PACBIO_CELL3 (BAMTOOLS_FILTER_PACBIO_CELL3.out.filtered_bam)
		}
		if( params.bam_cell4 ){
                        BAMTOOLS_FILTER_PACBIO_CELL4 (input_pacbio_cell4)
                        PBINDEX_FILTERED_PACBIO_CELL4 (BAMTOOLS_FILTER_PACBIO_CELL4.out.filtered_bam)
                }
	}

	//Merge the multiple pacbio bam files if multiple are generated and generate fastq files
	
	if( params.bam_cell4 ) {
                FOURBAM2FASTX(BAMTOOLS_FILTER_PACBIO_CELL1.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL1.out.index), BAMTOOLS_FILTER_PACBIO_CELL2.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL2.out.index), BAMTOOLS_FILTER_PACBIO_CELL3.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL3.out.index), BAMTOOLS_FILTER_PACBIO_CELL4.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL4.out.index))
                CUTADAPT (FOURBAM2FASTX.out.reads)
        } else if( params.bam_cell3 ) {
		THREEBAM2FASTX(BAMTOOLS_FILTER_PACBIO_CELL1.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL1.out.index), BAMTOOLS_FILTER_PACBIO_CELL2.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL2.out.index), BAMTOOLS_FILTER_PACBIO_CELL3.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL3.out.index))
		CUTADAPT (THREEBAM2FASTX.out.reads)
	} else if( params.bam_cell2 ){
		TWOBAM2FASTX(BAMTOOLS_FILTER_PACBIO_CELL1.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL1.out.index), BAMTOOLS_FILTER_PACBIO_CELL2.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL2.out.index))
		CUTADAPT (TWOBAM2FASTX.out.reads)
	} else {
		BAM2FASTX (BAMTOOLS_FILTER_PACBIO_CELL1.out.filtered_bam.join(PBINDEX_FILTERED_PACBIO_CELL1.out.index))
	        CUTADAPT (BAM2FASTX.out.reads)
        }

        //QC Input data
        mqc_input = Channel.empty()
        MERYL_COUNT (CUTADAPT.out.reads)
        MERYL_HISTOGRAM (MERYL_COUNT.out.meryl_db)
        GENOMESCOPE2 (MERYL_HISTOGRAM.out.hist)
        COVERAGE_CALCULATION(CUTADAPT.out.reads)

//All the following steps are commented out as they require some local installation to work
/*
        LONGQC (CUTADAPT.out.reads)
        KRAKEN2_KRAKEN2_PACBIO_BAM (CUTADAPT.out.reads, params.kraken_db, false, false )
	mqc_input = mqc_input.mix(KRAKEN2_KRAKEN2_PACBIO_BAM.out.report.collect{it[1]})
        COVERAGE_CALCULATION(CUTADAPT.out.reads)

	//ONLY if Hi-C data available
	if (( params.hic_read1 ) && (params.hic_read2 )) {
        	KRAKEN2_KRAKEN2_HIC_READS (input_hic_R1_R2, params.kraken_db, false, false )
        	mqc_input = mqc_input.mix(KRAKEN2_KRAKEN2_HIC_READS.out.report.collect{it[1]})
		kraken_hic = KRAKEN2_KRAKEN2_HIC_READS.out.report
	} else {
                kraken_hic = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('kraken_hic_dummy')]
                ]   
	}

        //ONLY if Illumina SR data available
        if ((params.illumina_SR_read1 ) && (params.illumina_SR_read2)) {
                KRAKEN2_KRAKEN2_SR_READS (input_illumina_SR_R1_R2, params.kraken_db, false, false )
                mqc_input = mqc_input.mix(KRAKEN2_KRAKEN2_SR_READS.out.report.collect{it[1]})
        }

        //ONLY if ONT data available
        if (params.ont_fastq_1) {
                KRAKEN2_KRAKEN2_ONT_READS (input_ont_fastq_1, params.kraken_db, false, false )
                mqc_input = mqc_input.mix(KRAKEN2_KRAKEN2_ONT_READS.out.report.collect{it[1]})
        }

	if (params.mitohifi == 'yes') {
	        //Mitochondrial assembly
		FASTQGZ_TO_FASTA(CUTADAPT.out.reads)
	        FIND_MITO_REFERENCE(FASTQGZ_TO_FASTA.out.fasta, params.taxon_name)
	        MITOHIFI(FASTQGZ_TO_FASTA.out.fasta, FIND_MITO_REFERENCE.out.reference_fasta, FIND_MITO_REFERENCE.out.reference_gb)
	}
*/
	//Assembly : The method is selected in the parameters : 'hifiasm' or 'flye' or 'canu' or 'verkko'
	if ( params.assembly_method == 'hifiasm') {
        	//HifiASM : Need to select a secondary mode : 'pacbio' or 'pacbio+hic' or 'pacbio+ont' or 'pacbio+ont+hic'
        	if (params.assembly_secondary_mode == 'pacbio+hic') {
			HIFIASM (CUTADAPT.out.reads, [], [], params.hic_read1, params.hic_read2, [] )
        	} else if (params.assembly_secondary_mode == 'pacbio+ont') {
                        HIFIASM (CUTADAPT.out.reads, [], [], [], [], params.ont_fastq_1 )
                } else if (params.assembly_secondary_mode == 'pacbio') {
			HIFIASM (CUTADAPT.out.reads, [], [], [], [], [] )
                } else if (params.assembly_secondary_mode == 'pacbio+ont+hic') {
                        HIFIASM (CUTADAPT.out.reads, [], [], params.hic_read1, params.hic_read2,params.ont_fastq_1 )
		} else {
			error "Invalid hifiasm mode: params.assembly_secondary_mode. These modes are currently supported : 'pacbio' or 'pacbio+hic' or 'pacbio+ont' or 'pacbio+ont+hic'"
		}
        	GFA_TO_FA (HIFIASM.out.hap1_contigs)
		assembly_primary = GFA_TO_FA.out.fa_assembly	
		GFA_TO_FA2 (HIFIASM.out.hap2_contigs)
		assembly_alternate = GFA_TO_FA2.out.fa_assembly
	} else if ( params.assembly_method == 'canu') {
		//CANU
		if (params.assembly_secondary_mode == 'hicanu') {
                        CANU(CUTADAPT.out.reads)
                } else if (params.assembly_secondary_mode == 'ont') {
                        CANU(input_ont_fastq_1)
                } else {
			error "Invalid canu mode: params.assembly_secondary_mode. These modes are currently supported : 'hicanu' or 'ont'"
		}
                assembly_primary = CANU.out.assembly
	} else if ( params.assembly_method == 'flye') {
        	//FLYE
		if (params.assembly_secondary_mode== 'hifi') {
			mode = "--pacbio-hifi"
        		FLYE (CUTADAPT.out.reads, mode)
			MINIMAP_ALIGN_FLYE (CUTADAPT.out.reads, FLYE.out.fasta.collect{it[1]}, false, false, false)
			RACON (CUTADAPT.out.reads, FLYE.out.fasta.join (MINIMAP_ALIGN_FLYE.out.paf))
                	LONGSTITCH (CUTADAPT.out.reads, RACON.out.improved_assembly)
        	} else if (params.assembly_secondary_mode== 'ont') {
                        mode = "--nano-raw"
                        FLYE (input_ont_fastq_1, mode)
			MINIMAP_ALIGN_FLYE (input_ont_fastq_1, FLYE.out.fasta.collect{it[1]}, false, false, false)
                        RACON (input_ont_fastq_1, FLYE.out.fasta.join (MINIMAP_ALIGN_FLYE.out.paf))
                        LONGSTITCH (input_ont_fastq_1, RACON.out.improved_assembly)
                } else if (params.assembly_secondary_mode== 'pacbio+ont') {
                        FLYE_PACBIO_ONT (CUTADAPT.out.reads, input_ont_fastq_1)
                        MINIMAP_ALIGN_FLYE (input_ont_fastq_1, FLYE_PACBIO_ONT.out.fasta.collect{it[1]}, false, false, false)
                        RACON (input_ont_fastq_1, FLYE_PACBIO_ONT.out.fasta.join (MINIMAP_ALIGN_FLYE.out.paf))
                        LONGSTITCH (input_ont_fastq_1, RACON.out.improved_assembly)
                } else {
                        error "Invalid flye mode: params.assembly_secondary_mode. These modes are currently supported : 'hifi' or 'ont' or 'pacbio+ont'"
		}

		assembly_primary = LONGSTITCH.out.assembly
        } else if ( params.assembly_method == 'verkko') {
                //VERKKO
		if (params.assembly_secondary_mode == 'pacbio+ont') {
			VERKKO(CUTADAPT.out.reads, input_ont_fastq_1)
		} else if (params.assembly_secondary_mode == 'pacbio') {
			VERKKO(CUTADAPT.out.reads, input_dummy)
		} else if (params.assembly_secondary_mode == 'ont') {
			VERKKO(input_dummy, input_ont_fastq_1)
		} else {
                        error "Invalid verkko mode: params.assembly_secondary_mode. These modes are currently supported : 'pacbio', 'ont', 'pacbio+ont'"
		}
                assembly_primary = VERKKO.out.assembly
	} else {
		error "Invalid alignment method: params.assembly_method. These methods are currently supported : 'hifiasm', 'canu', 'flye', 'verkko'. "
	}

	//QC post assembly
	if ((params.assembly_method == 'hifiasm')&& (params.ploidy != '1')) {
	        QUAST1_DOUBLE (assembly_primary, assembly_alternate, quast_fasta, quast_gff, false, false )
	        mqc_input = mqc_input.mix(QUAST1_DOUBLE.out.tsv)
		quast_contig = QUAST1_DOUBLE.out.tsv
		MERQURY1_DOUBLE (MERYL_COUNT.out.meryl_db, assembly_primary, assembly_alternate)
	} else {
                QUAST1 (assembly_primary, quast_fasta, quast_gff, false, false )
                mqc_input = mqc_input.mix(QUAST1.out.tsv)
                quast_contig = QUAST1.out.tsv
		MERQURY1 (MERYL_COUNT.out.meryl_db.join(assembly_primary))
	}

	//Polishing (likely not going to happen, only for primary assembly for now)
	if (params.polishing_method == 'pilon') {
		BWAMEM2_INDEX(assembly_primary)
		BWAMEM2_MEM(input_illumina_SR_R1_R2, BWAMEM2_INDEX.out.index, true)
		SAMTOOLS_INDEX(BWAMEM2_MEM.out.bam)
		pilon_mode="--frags"
		PILON(assembly_primary, pilon_mode, BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai))
		assembly_polished = PILON.out.improved_assembly
		QUAST_PILON(assembly_polished, quast_fasta, quast_gff, false, false )
		mqc_input = mqc_input.mix(QUAST_PILON.out.tsv)
		assembly_unpurged = assembly_polished
	} else {
		assembly_unpurged = assembly_primary
	}
		

	//PurgeDups for primary assembly
        PURGEDUPS_SPLITFA (assembly_unpurged)
        MINIMAP2_ALIGN_TO_CONTIG (CUTADAPT.out.reads, assembly_unpurged.collect{it[1]}, false, false, false)
        MINIMAP2_ALIGN_TO_SELF (PURGEDUPS_SPLITFA.out.split_fasta, [], false, false, false)
        PURGEDUPS_PBCSTAT (MINIMAP2_ALIGN_TO_CONTIG.out.paf)
        PURGEDUPS_CALCUTS (PURGEDUPS_PBCSTAT.out.stat) 
        PURGEDUPS_PURGEDUPS (
                PURGEDUPS_PBCSTAT.out.basecov
                        .join (PURGEDUPS_CALCUTS.out.cutoff )
                        .join (MINIMAP2_ALIGN_TO_SELF.out.paf )
                )
        PURGEDUPS_GETSEQS (assembly_unpurged.join(PURGEDUPS_PURGEDUPS.out.bed))
        SAMTOOLS_FAIDX1 (PURGEDUPS_GETSEQS.out.purged)
	purged_primary = PURGEDUPS_GETSEQS.out.purged

        if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
                //Merge haplotig from purge_dups and alternate assembly from hifiasm
                CAT (assembly_alternate, PURGEDUPS_GETSEQS.out.haplotigs)
                PURGEDUPS_SPLITFA_ALT (CAT.out.alternate_contigs_full)
                MINIMAP2_ALIGN_TO_CONTIG_ALT (CUTADAPT.out.reads, CAT.out.alternate_contigs_full.collect{it[1]}, false, false, false)
                MINIMAP2_ALIGN_TO_SELF_ALT (PURGEDUPS_SPLITFA_ALT.out.split_fasta, [], false, false, false)
                PURGEDUPS_PBCSTAT_ALT (MINIMAP2_ALIGN_TO_CONTIG_ALT.out.paf)
                PURGEDUPS_CALCUTS_ALT (PURGEDUPS_PBCSTAT_ALT.out.stat)
                PURGEDUPS_PURGEDUPS_ALT (
                        PURGEDUPS_PBCSTAT_ALT.out.basecov
                                .join (PURGEDUPS_CALCUTS_ALT.out.cutoff )
                                .join (MINIMAP2_ALIGN_TO_SELF_ALT.out.paf )
                        )
                PURGEDUPS_GETSEQS_ALT (assembly_alternate.join(PURGEDUPS_PURGEDUPS_ALT.out.bed))
                SAMTOOLS_FAIDX1_ALT (PURGEDUPS_GETSEQS_ALT.out.purged)
                purged_alternate = PURGEDUPS_GETSEQS_ALT.out.purged

                QUAST2_DOUBLE (purged_primary, purged_alternate, quast_fasta, quast_gff, false, false )
                mqc_input = mqc_input.mix(QUAST2_DOUBLE.out.tsv)
		quast_contig_purged=QUAST2_DOUBLE.out.tsv
                MERQURY2_DOUBLE(MERYL_COUNT.out.meryl_db, purged_primary, purged_alternate)
        } else {
		QUAST2 (purged_primary, quast_fasta, quast_gff, false, false  )
		mqc_input = mqc_input.mix(QUAST2.out.tsv)
		quast_contig_purged = QUAST2.out.tsv
                MERQURY2(MERYL_COUNT.out.meryl_db.join(purged_primary))
        }


	//Only if HiC data is available
	if (( params.hic_read1 ) && ( params.hic_read2 )) {
		//HIC scaffolding: The method is selected in the parameters : salsa or yahs
		if ( params.scaffolding_method == "salsa") {
			// For SALSA2, need to run nf-core/hic
			PREPARE_GENOME (PURGEDUPS_GETSEQS.out.purged, params.restriction_site)
			FASTQC (input_hic_R1_R2)
			ch_map_res = Channel.from( params.bin_size ).splitCsv().flatten().toInteger()
			HICPRO (input_hic_R1_R2, PREPARE_GENOME.out.index, PREPARE_GENOME.out.res_frag, PREPARE_GENOME.out.chromosome_size, params.ligation_site, ch_map_res)
			COOLER (HICPRO.out.pairs, PREPARE_GENOME.out.chromosome_size, ch_map_res)
			ch_ddecay_res = Channel.empty()
			COOLER.out.cool
				.combine(ch_ddecay_res)
				.filter{ it[0].resolution == it[2] }
				.map { it -> [it[0], it[1]]}
				.set{ ch_distdecay }
			HIC_PLOT_DIST_VS_COUNTS (ch_distdecay)
			COOLER.out.cool
				.combine(ch_comp_res)
				.filter{ it[0].resolution == it[2] }
				.map { it -> [it[0], it[1], it[2]]}
				.set{ ch_cool_compartments }
			COMPARTMENTS (ch_cool_compartments, PURGEDUPS_GETSEQS.out.purged, PREPARE_GENOME.out.chromosome_size)
		 	COOLER.out.cool
				.combine(ch_tads_res)
				.filter{ it[0].resolution == it[2] }
				.map { it -> [it[0], it[1]]}
				.set{ ch_cool_tads }
			TADS (ch_cool_tads)
			mqc_input = mqc_input.mix(ch_multiqc_config)
                        mqc_input = mqc_input.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
                        mqc_input = mqc_input.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
                        mqc_input = mqc_input.mix(FASTQC.out.zip.map{it->it[1]})
                        mqc_input = mqc_input.mix(HICPRO.out.mqc)
			// ln -s PATH/TO/HIC_READ1 hic_R1.fastq.gz
			// ln -s PATH/TO/HIC_READ2 hic_R2.fastq.gz
			// nextflow run nf-core/hic -profile singularity --input 'hic_R{1,2}.fastq.gz' --fasta purge_dups/*.purged.fa --ligation_site 'GATCGATC,GANTGATC,GANTANTC,GATCANTC'       --restriction_site '^GATC,G^ANTC,C^TNAG,T^TAA'

			//Then start the pipeline again.
        		BED_PROCESSING(HICPRO.out.bam)
        		SALSA2 (PURGEDUPS_GETSEQS.out.purged.join(SAMTOOLS_FAIDX1.out.fai), BED_PROCESSING.out.sorted_bed.collect{it[1]}, [], [], [] )
			SAMTOOLS_FAIDX2 (SALSA2.out.fasta)
        		scaffold                = SALSA2.out.fasta
        		scaffold_agp            = SALSA2.out.agp
        		scaffold_index          = SAMTOOLS_FAIDX2.out.fai
	
			SALSA2_JUICER (
                		scaffold
					.join(scaffold_index)
                	        	.join(scaffold_agp)
                	        	.join(SALSA2.out.alignment_iteration_1_bed)
                	        	.join(SALSA2.out.scaffold_length_iteration_1)
                		)	
		} else if ( params.scaffolding_method == "yahs") {
        		CHROMAP_INDEX(PURGEDUPS_GETSEQS.out.purged)
        		CHROMAP_CHROMAP(input_hic_R1_R2, PURGEDUPS_GETSEQS.out.purged, CHROMAP_INDEX.out.index, [],[],[],[])
        		YAHS(PURGEDUPS_GETSEQS.out.purged, SAMTOOLS_FAIDX1.out.fai, CHROMAP_CHROMAP.out.bam)
        		SAMTOOLS_FAIDX2(YAHS.out.fasta)
	
        		scaffold 		= YAHS.out.fasta
			scaffold_agp 		= YAHS.out.agp
			scaffold_bin 		= YAHS.out.bin
        		scaffold_index		= SAMTOOLS_FAIDX2.out.fai

		        if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
				CHROMAP_INDEX_ALT(purged_alternate)
				CHROMAP_CHROMAP_ALT(input_hic_R1_R2, purged_alternate, CHROMAP_INDEX_ALT.out.index, [],[],[],[])
                        	YAHS_ALT(purged_alternate, SAMTOOLS_FAIDX1_ALT.out.fai, CHROMAP_CHROMAP_ALT.out.bam)
                        	SAMTOOLS_FAIDX2_ALT(YAHS_ALT.out.fasta)

	                        scaffold_alt                = YAHS_ALT.out.fasta
	                        scaffold_agp_alt            = YAHS_ALT.out.agp
	                        scaffold_bin_alt            = YAHS_ALT.out.bin
	                        scaffold_index_alt          = SAMTOOLS_FAIDX2_ALT.out.fai
			}
		} else {
			error "Invalid alignment mode: params.scaffolding_method "
		}

		//Scaffold QC
        	if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
                	QUAST3_DOUBLE (scaffold, scaffold_alt, quast_fasta, quast_gff, false, false )
                	mqc_input = mqc_input.mix(QUAST3_DOUBLE.out.tsv)
			quast_scaffold = QUAST3_DOUBLE.out.tsv
        		MERQURY3_DOUBLE(MERYL_COUNT.out.meryl_db, scaffold, scaffold_alt)
		} else {
	                QUAST3 (scaffold, quast_fasta, quast_gff, false, false )
	                mqc_input = mqc_input.mix(QUAST3.out.tsv)
			quast_scaffold = QUAST3.out.tsv
			MERQURY3(MERYL_COUNT.out.meryl_db.join(scaffold))
	        }
	
		// JUICER must have contig fai for scaffold assembly
	        YAHS_JUICER (scaffold_agp.join(scaffold_bin), SAMTOOLS_FAIDX1.out.fai)
		chrom_size = YAHS_JUICER.out.chrom_sizes
//        	JUICER(YAHS_JUICER.out.chrom_sizes, YAHS_JUICER.out.alignments_sorted_txt)
//		PRETEXTMAP(YAHS_JUICER.out.chrom_sizes, YAHS_JUICER.out.alignments_sorted_txt)
//        	PRETEXTSNAPSHOT (PRETEXTMAP.out.pretext)
	
//Blobtoolkit is commented out as it requires some local installation
/*
        	//BLOBTOOLSKIT
		GZIP(scaffold)
	        if( params.bam_cell4 ) {
                        BLOBTOOLS_CONFIG(GZIP.out.gz, FOURBAM2FASTX.out.reads)
                } else if( params.bam_cell3 ) {
			BLOBTOOLS_CONFIG(GZIP.out.gz, THREEBAM2FASTX.out.reads)
	        } else if( params.bam_cell2 ){
			BLOBTOOLS_CONFIG(GZIP.out.gz, TWOBAM2FASTX.out.reads)
		} else {
			BLOBTOOLS_CONFIG(GZIP.out.gz, BAM2FASTX.out.reads)
		}
		BLOBTOOLS_PIPELINE(BLOBTOOLS_CONFIG.out.config, GZIP.out.gz)
		BLOBTOOLS_CREATE(scaffold, BLOBTOOLS_CONFIG.out.config)
		BLOBTOOLS_ADD(BLOBTOOLS_PIPELINE.out.blast_out, BLOBTOOLS_PIPELINE.out.diamond_proteome_out, BLOBTOOLS_PIPELINE.out.diamond_busco_out, BLOBTOOLS_PIPELINE.out.assembly_minimap_bam, BLOBTOOLS_PIPELINE.out.hic_minimap_bam , BLOBTOOLS_PIPELINE.out.lineage1_full_table_tsv , BLOBTOOLS_PIPELINE.out.lineage2_full_table_tsv, BLOBTOOLS_CREATE.out.blobtools_folder)
		BLOBTOOLS_VIEW_SNAIL(BLOBTOOLS_ADD.out.blobtools_folder)
		BLOBTOOLS_VIEW_BLOB(BLOBTOOLS_ADD.out.blobtools_folder)
		BLOBTOOLS_VIEW_CUMULATIVE(BLOBTOOLS_ADD.out.blobtools_folder)
*/
	} else {
		quast_scaffold = file('quast_scaffold_dummy')
                chrom_size = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('chrom_size_dummy')]
                ]
	}


//BUSCO is commented out as it requires local database
/*
	//Busco ran only once on the most final assembly (scaffold > purged > contig)
	//Define which file to run Busco on
	if (( params.hic_read1 ) && ( params.hic_read2 )) {
		busco_assembly = scaffold
		if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
			busco_assembly_alt = scaffold_alt
		}
	} else {
		busco_assembly = purged_primary
		if (( params.assembly_method == 'hifiasm' ) && (params.ploidy != '1')) {
			busco_assembly_alt = purged_alternate
		}
	}

        BUSCO (busco_assembly, params.lineage, params.busco_lineages_path, [])
        mqc_input = mqc_input.mix(BUSCO.out.short_summaries_txt.collect{it[1]})
        if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
                BUSCO_ALT (busco_assembly_alt, params.lineage, params.busco_lineages_path, [])
                mqc_input = mqc_input.mix(BUSCO_ALT.out.short_summaries_txt.collect{it[1]})
        }

        if (params.lineage2) {
                BUSCO_lin2(busco_assembly, params.lineage2, params.busco_lineages_path, [])
                mqc_input = mqc_input.mix(BUSCO_lin2.out.short_summaries_txt.collect{it[1]})
                busco_lin2_json = BUSCO_lin2.out.short_summaries_json
//		if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
//                        BUSCO_lin2ALT (busco_assembly_alt, params.lineage2, params.busco_lineages_path, [])
//                        mqc_input = mqc_input.mix(BUSCO_lin2ALT.out.short_summaries_txt.collect{it[1]})
//                }
        } else {
                busco_lin2_json = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('busco_lin2_json_dummy')]
                ]
	}
        if (params.lineage3) {
                BUSCO_lin3(busco_assembly, params.lineage3, params.busco_lineages_path, [])
                mqc_input = mqc_input.mix(BUSCO_lin3.out.short_summaries_txt.collect{it[1]})
		busco_lin3_json = BUSCO_lin3.out.short_summaries_json
//                if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
//                        BUSCO_lin3ALT (busco_assembly_alt, params.lineage3, params.busco_lineages_path, [])
//                        mqc_input = mqc_input.mix(BUSCO_lin3ALT.out.short_summaries_txt.collect{it[1]})
//                }
        } else {
                busco_lin3_json = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('busco_lin3_json_dummy')]
                ]
	}
        if (params.lineage4) {
                BUSCO_lin4(busco_assembly, params.lineage4, params.busco_lineages_path, [])
                mqc_input = mqc_input.mix(BUSCO_lin4.out.short_summaries_txt.collect{it[1]})
		busco_lin4_json = BUSCO_lin4.out.short_summaries_json
//                if ((params.assembly_method == 'hifiasm') && (params.ploidy != '1')) {
//                        BUSCO_lin4ALT (busco_assembly_alt, params.lineage4, params.busco_lineages_path, [])
//                        mqc_input = mqc_input.mix(BUSCO_lin4ALT.out.short_summaries_txt.collect{it[1]})
//                }
        } else {
                busco_lin4_json = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('busco_lin4_json_dummy')]
                ]
	}
*/

	//MultiQC report
	MULTIQC (mqc_input.collect(), [], [], [], COVERAGE_CALCULATION.out.coverage)
}


