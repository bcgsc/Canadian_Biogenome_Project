#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
CBP pipeline - Solenne Correard - Jones lab
=============================================
Specie id			: ${params.id}
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
FCS (Foreign Conta Screen)	: ${params.fcs}
Polishing method		: ${params.polishing_method}
Purging method			: ${params.purging_method}
Scaffolding method		: ${params.scaffolding_method}
Mitochondrial assembly		: ${params.mitohifi}
Pretext Hi-C map		: ${params.pretext}
Juicer Hi-C map			: ${params.juicer}
Methylation calling		: ${params.methylation_calling}
Comparison to related genome	: ${params.genome_comparison}
Blobtools			: ${params.blobtools}
Manual curation                 : ${params.manual_curation}
"""

include { GOAT_TAXONSEARCH } from './modules/goat/taxonsearch/main.nf'

//Pre-processing
include { CCS as CCS_PACBIO } from './modules/pacbio/ccs/main.nf'
include { BAMTOOLS_FILTER as BAMTOOLS_FILTER_PACBIO } from './modules/bamtools_filter/main.nf'
include { PBINDEX as PBINDEX_FILTERED_PACBIO } from './modules/pacbio/pbbam/pbindex/main.nf'
include { PBBAM_PBMERGE } from './modules/pacbio/pbbam/pbmerge/main.nf'
include { BAM2FASTX } from './modules/pacbio/bam2fastx/main.nf'

include { PREPROCESS_MERGED } from './modules/pacbio/preprocess_merged/main.nf'

include { CUTADAPT }  from './modules/cutadapt/main.nf'

//QC Input data
include { LONGQC as LONGQC_PACBIO; LONGQC as LONGQC_ONT } from './modules/LongQC/main.nf'
include { MERYL_COUNT } from './modules/meryl/count/main.nf'
include { MERYL_UNIONSUM } from './modules/meryl/unionsum/main.nf'
include { MERYL_HISTOGRAM } from './modules/meryl/histogram/main.nf'
include { GENOMESCOPE2 } from './modules/genomescope2/main.nf'
include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_PACBIO_BAM; KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_HIC_READS; KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_SR_READS; KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ONT_READS } from './modules/kraken2/main.nf'
include { COVERAGE_CALCULATION } from './modules/coverage_calculation/main.nf'

//Mitochondrial assembly
include { FASTQGZ_TO_FASTA } from './modules/fastqgz_to_fasta/main.nf'
include { FIND_MITO_REFERENCE } from './modules/mitohifi/findmitoreference/main.nf'
include { MITOHIFI } from './modules/mitohifi/mitohifi/main.nf'


//Assembly
//HifiASM
include { HIFIASM } from './modules/hifiasm/main.nf'
include { GFA_TO_FA as GFA_TO_FA_hap1; GFA_TO_FA as GFA_TO_FA_hap2 } from './modules/gfa_to_fa/main.nf'

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
include { BWAMEM2_INDEX as BWAMEM2_INDEX_PILON } from './modules/bwamem2/index/main.nf'
include { BWAMEM2_MEM as BWAMEM2_MEM_PILON } from './modules/bwamem2/mem/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PILON } from './modules/samtools/index/main.nf'

//NCBI cleaning sequence
include { FCS_FCSADAPTOR as FCS_FCSADAPTOR_hap1; FCS_FCSADAPTOR as FCS_FCSADAPTOR_ALT } from './modules/fcs/fcsadaptor/'
include { FCS_FCSGX as FCS_FCSGX_hap1; FCS_FCSGX as FCS_FCSGX_ALT } from './modules/fcs/fcsgx'
include { FCS_FCSGX_CLEAN as FCS_FCSGX_CLEAN_hap1; FCS_FCSGX_CLEAN as FCS_FCSGX_CLEAN_ALT } from './modules/fcs/fcsgx_clean'


//PurgeDups
include { CAT } from './modules/cat/main.nf'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_CONTIG; MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_SELF; MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_CONTIG_ALT; MINIMAP2_ALIGN as MINIMAP2_ALIGN_TO_SELF_ALT } from './modules/minimap2/align/main.nf'
include { PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_hap1; PURGEDUPS_SPLITFA as PURGEDUPS_SPLITFA_ALT } from './modules/purgedups/splitfa/main.nf'
include { PURGEDUPS_PBCSTAT as PURGEDUPS_PBCSTAT_hap1; PURGEDUPS_PBCSTAT as PURGEDUPS_PBCSTAT_ALT } from './modules/purgedups/pbcstat/main.nf'
include { PURGEDUPS_CALCUTS as PURGEDUPS_CALCUTS_hap1; PURGEDUPS_CALCUTS as PURGEDUPS_CALCUTS_ALT } from './modules/purgedups/calcuts/main.nf'
include { PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_hap1; PURGEDUPS_PURGEDUPS as PURGEDUPS_PURGEDUPS_ALT } from './modules/purgedups/purgedups/main.nf'
include { PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_hap1; PURGEDUPS_GETSEQS as PURGEDUPS_GETSEQS_ALT } from './modules/purgedups/getseqs/main.nf'

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
include { CHROMAP_INDEX as CHROMAP_INDEX_hap1; CHROMAP_INDEX as CHROMAP_INDEX_ALT } from './modules/chromap/index/main.nf'
include { CHROMAP_CHROMAP as CHROMAP_CHROMAP_hap1; CHROMAP_CHROMAP as CHROMAP_CHROMAP_ALT } from './modules/chromap/chromap/main.nf'
include { YAHS as YAHS_hap1; YAHS as YAHS_ALT } from './modules/yahs/main.nf'

//Map PacBio data against newly genevated assembly
include { JASMINE } from './modules/pacbio/jasmine/main.nf'
include { PBMM2 } from './modules/pacbio/pbmm2/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PBMM2 } from './modules/samtools/index/main.nf'

//Assembly QC
include { YAHS_JUICER } from './modules/juicer/yahs_juicer/main.nf'
include { JUICER } from './modules/juicer/juicer/main.nf'
include { PRETEXTMAP } from './modules/pretext/pretextmap/main.nf'
include { PRETEXTSNAPSHOT } from './modules/pretext/pretextsnapshot/main.nf'
include { BEDTOOLS_GENOMECOV } from './modules/bedtools/genomecov/main.nf'
include { GFASTATS } from './modules/gfastats/main.nf'
include { TIDK } from './modules/tidk/main.nf'
include { PRETEXTGRAPH as PRETEXTGRAPH_TELO; PRETEXTGRAPH as PRETEXTGRAPH_TELO_COV } from './modules/pretext/pretextgraph/main.nf'

//Genome comparison
include { NCBIGENOMEDOWNLOAD } from './modules/ncbigenomedownload/main.nf'
include { JUPITER } from './modules/jupiter/main.nf'
include { MASHMAP } from './modules/mashmap/main.nf'

//QC of assemblies
include { BUSCO as BUSCO_lin1_PRIM; BUSCO as BUSCO_lin1_cleaned; BUSCO as BUSCO_lin1_purged; BUSCO as BUSCO_lin1_SCAFF; BUSCO as BUSCO_lin2; BUSCO as BUSCO_lin3; BUSCO as BUSCO_lin4; BUSCO as BUSCO_ALT } from './modules/busco/main.nf'
include { MERQURY as MERQURY_ASS; MERQURY as MERQURY_PURGED; MERQURY as MERQURY_SCAFF } from './modules/merqury/main.nf'
include { MERQURY_DOUBLE as MERQURY_ASS_DOUBLE; MERQURY_DOUBLE as MERQURY_PURGED_DOUBLE; MERQURY_DOUBLE as MERQURY_SCAFF_DOUBLE } from './modules/merqury/merqury_double/main.nf'
include { QUAST as QUAST_ASS; QUAST as QUAST_PILON; QUAST as QUAST_CLEAN; QUAST as QUAST_PURGED; QUAST as QUAST_SCAFF } from './modules/quast/main.nf'
include { QUAST_DOUBLE as QUAST_ASS_DOUBLE; QUAST_DOUBLE as QUAST_CLEAN_DOUBLE; QUAST_DOUBLE as QUAST_PURGED_DOUBLE; QUAST_DOUBLE as QUAST_SCAFF_DOUBLE } from './modules/quast/quast_double/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

//Blobtoolskit
include { GZIP } from './modules/gzip/main.nf'
include { BLOBTOOLS_CONFIG } from './modules/blobtools/blobtools_config/main.nf'
include { BLOBTOOLS_PIPELINE } from './modules/blobtools/blobtools_pipeline/main.nf'
include { BLOBTOOLS_CREATE } from './modules/blobtools/blobtools_create/main.nf'
include { BLOBTOOLS_ADD } from './modules/blobtools/blobtools_add/main.nf'
include { BLOBTOOLS_VIEW } from './modules/blobtools/blobtools_view/main.nf'

include { OVERVIEW_GENERATION_SAMPLE } from './modules/overview_generation/sample/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from './modules/dumpsoftwareversions/main'

include { RAPIDCURATION_SPLIT } from './modules/manualcuration/main.nf'

workflow {

//////////////////////////////////////////////////  INPUT   //////////////////////////////////

taxon  = [
      [ id:params.id ], // meta map
      taxon = params.taxon_taxid,
      []
  ]

//PacBio data

        if( params.bam_cell4 ){
                input_pacbio = [
                        [ id:params.id, single_end: true], // meta map
                        [
                                file(params.bam_cell1, checkIfExists: true),
                                file(params.bam_cell2, checkIfExists: true),
                                file(params.bam_cell3, checkIfExists: true),
                                file(params.bam_cell4, checkIfExists: true)

                         ]
                ]
        } else if( params.bam_cell3 ){
                input_pacbio = [
                        [ id:params.id, single_end: true], // meta map
                        [
                                file(params.bam_cell1, checkIfExists: true),
                                file(params.bam_cell2, checkIfExists: true),
                                file(params.bam_cell3, checkIfExists: true)
                         ]
                ]
        } else if( params.bam_cell2 ){
                input_pacbio = [
                        [ id:params.id, single_end: true], // meta map
                        [
                                file(params.bam_cell1, checkIfExists: true),
                                file(params.bam_cell2, checkIfExists: true)
                         ]
                ]
        } else {
	        input_pacbio = [
	                [ id:params.id, single_end: true], // meta map
	                [ file(params.bam_cell1, checkIfExists: true) ]
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

    // To gather all QC reports for MultiQC
    mqc_input  = Channel.empty()
    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()

	GOAT_TAXONSEARCH(taxon)

//PacBio data is very large
//When there is several HIFI SMRT cells, 
//doing the merging, filtering and bam2fastx steps in different modules is space consumming
//Limiting the capacity to run several genomes in parrallel
//Merging the steps in one module and deleting the intermediate files was a solution

        if ((params.bam_cell2) && (params.pacbio_input_type == 'ccs')){
		//MERGED STEPS : PBBAM_PBMERGE + BAMTOOLS_FILTER_PACBIO 
                PREPROCESS_MERGED(input_pacbio)
                bamtools_filter_output= PREPROCESS_MERGED.out.filtered_bam
                PBINDEX_FILTERED_PACBIO (bamtools_filter_output)
		ch_versions = ch_versions.mix(PREPROCESS_MERGED.out.versions)
                ch_versions = ch_versions.mix(PBINDEX_FILTERED_PACBIO.out.versions)
	} else if ((params.bam_cell2) && (params.pacbio_input_type == 'subreads')) {
                PBBAM_PBMERGE(input_pacbio)
		CCS_PACBIO(PBBAM_PBMERGE.out.bam)
		BAMTOOLS_FILTER_PACBIO (CCS_PACBIO_CELL1.out.bam)
	        PBINDEX_FILTERED_PACBIO (BAMTOOLS_FILTER_PACBIO.out.filtered_bam)
		bamtools_filter_output=BAMTOOLS_FILTER_PACBIO.out.filtered_bam
                ch_versions = ch_versions.mix(PBBAM_PBMERGE.out.versions)
                ch_versions = ch_versions.mix(CCS_PACBIO_CELL1.out.versions)
                ch_versions = ch_versions.mix(BAMTOOLS_FILTER_PACBIO.out.versions)
                ch_versions = ch_versions.mix(PBINDEX_FILTERED_PACBIO.out.versions)
	} else {
		//If only one SMRT cell
		BAMTOOLS_FILTER_PACBIO (input_pacbio)
                PBINDEX_FILTERED_PACBIO (BAMTOOLS_FILTER_PACBIO.out.filtered_bam)
		bamtools_filter_output=BAMTOOLS_FILTER_PACBIO.out.filtered_bam
                ch_versions = ch_versions.mix(BAMTOOLS_FILTER_PACBIO.out.versions)
                ch_versions = ch_versions.mix(PBINDEX_FILTERED_PACBIO.out.versions)
	}
	BAM2FASTX (bamtools_filter_output.join(PBINDEX_FILTERED_PACBIO.out.index))
        bam2fastx_output=BAM2FASTX.out.reads
	CUTADAPT (bam2fastx_output)
        // Gather versions of all tools used
	ch_versions = ch_versions.mix(CUTADAPT.out.versions)
	ch_versions = ch_versions.mix(BAM2FASTX.out.versions)

        //QC Input data
        LONGQC_PACBIO (CUTADAPT.out.reads)
        MERYL_COUNT (CUTADAPT.out.reads)
        MERYL_HISTOGRAM (MERYL_COUNT.out.meryl_db)
        GENOMESCOPE2 (MERYL_HISTOGRAM.out.hist, GOAT_TAXONSEARCH.out.ploidy)
        COVERAGE_CALCULATION(CUTADAPT.out.reads, GOAT_TAXONSEARCH.out.genome_size)

	if (params.execute_kraken == 'yes') {
	        KRAKEN2_KRAKEN2_PACBIO_BAM (CUTADAPT.out.reads, params.kraken_db, false, false )
		mqc_input = mqc_input.mix(KRAKEN2_KRAKEN2_PACBIO_BAM.out.report.collect{it[1]})
		kraken_pacbio = KRAKEN2_KRAKEN2_PACBIO_BAM.out.report
	        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2_PACBIO_BAM.out.versions)
	} else {
		kraken_pacbio = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('kraken_pacbio_dummy')]
                ]
        }


        // Gather versions of all tools used
        ch_versions = ch_versions.mix(LONGQC_PACBIO.out.versions)
        ch_versions = ch_versions.mix(MERYL_COUNT.out.versions)
        ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions)

	//ONLY if Hi-C data available
	if (( params.hic_read1 ) && (params.hic_read2 ) && (params.execute_kraken == 'yes')) {
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
		LONGQC_ONT(input_ont_fastq_1)
        }

	if (params.mitohifi == 'yes') {
	        //Mitochondrial assembly
		FASTQGZ_TO_FASTA(CUTADAPT.out.reads)
	        FIND_MITO_REFERENCE(FASTQGZ_TO_FASTA.out.fasta, GOAT_TAXONSEARCH.out.scientific_name)
	        MITOHIFI(FASTQGZ_TO_FASTA.out.fasta, FIND_MITO_REFERENCE.out.reference_fasta, FIND_MITO_REFERENCE.out.reference_gb)

	        // Gather versions of all tools used
	        ch_versions = ch_versions.mix(FASTQGZ_TO_FASTA.out.versions)
                ch_versions = ch_versions.mix(FIND_MITO_REFERENCE.out.versions)
                ch_versions = ch_versions.mix(MITOHIFI.out.versions)
	}

	//Assembly : The method is selected in the parameters : 'hifiasm' or 'flye' or 'canu' or 'verkko'
	if ( params.assembly_method == 'hifiasm') {
        	//HifiASM : Need to select a secondary mode : 'pacbio' or 'pacbio+hic' or 'pacbio+ont' or 'pacbio+ont+hic'
        	if (params.assembly_secondary_mode == 'pacbio+hic') {
			HIFIASM (CUTADAPT.out.reads, [], [], params.hic_read1, params.hic_read2, [], GOAT_TAXONSEARCH.out.ploidy, GOAT_TAXONSEARCH.out.genome_size )
        	} else if (params.assembly_secondary_mode == 'pacbio+ont') {
                        HIFIASM (CUTADAPT.out.reads, [], [], [], [], params.ont_fastq_1, GOAT_TAXONSEARCH.out.ploidy, GOAT_TAXONSEARCH.out.genome_size )
                } else if (params.assembly_secondary_mode == 'pacbio') {
			HIFIASM (CUTADAPT.out.reads, [], [], [], [], [], GOAT_TAXONSEARCH.out.ploidy, GOAT_TAXONSEARCH.out.genome_size )
                } else if (params.assembly_secondary_mode == 'pacbio+ont+hic') {
                        HIFIASM (CUTADAPT.out.reads, [], [], params.hic_read1, params.hic_read2,params.ont_fastq_1, GOAT_TAXONSEARCH.out.ploidy, GOAT_TAXONSEARCH.out.genome_size )
		} else {
			error "Invalid hifiasm mode: params.assembly_secondary_mode. These modes are currently supported : 'pacbio' or 'pacbio+hic' or 'pacbio+ont' or 'pacbio+ont+hic'"
		}
        	GFA_TO_FA_hap1 (HIFIASM.out.hap1_contigs)
		assembly_primary = GFA_TO_FA_hap1.out.fa_assembly	
		GFA_TO_FA_hap2 (HIFIASM.out.hap2_contigs)
		assembly_alternate = GFA_TO_FA_hap2.out.fa_assembly

                // Gather versions of all tools used
                ch_versions = ch_versions.mix(HIFIASM.out.versions)
                ch_versions = ch_versions.mix(GFA_TO_FA_hap1.out.versions)
	} else if ( params.assembly_method == 'canu') {
		//CANU
		if (params.assembly_secondary_mode == 'hicanu') {
                        CANU(CUTADAPT.out.reads, GOAT_TAXONSEARCH.out.genome_size)
                } else if (params.assembly_secondary_mode == 'ont') {
                        CANU(input_ont_fastq_1, GOAT_TAXONSEARCH.out.genome_size)
                } else if (params.assembly_secondary_mode == 'clr') {
                        CANU(CUTADAPT.out.reads, GOAT_TAXONSEARCH.out.genome_size)
                } else {
			error "Invalid canu mode: params.assembly_secondary_mode. These modes are currently supported : 'hicanu', 'ont' or 'clr'"
		}
                assembly_primary = CANU.out.assembly

                // Gather versions of all tools used
                ch_versions = ch_versions.mix(CANU.out.versions)
	} else if ( params.assembly_method == 'flye') {
        	//FLYE
		if (params.assembly_secondary_mode== 'hifi') {
			mode = "--pacbio-hifi"
        		FLYE (CUTADAPT.out.reads, mode)
			MINIMAP_ALIGN_FLYE (CUTADAPT.out.reads, FLYE.out.fasta.collect{it[1]}, false, false, false)
			RACON (CUTADAPT.out.reads, FLYE.out.fasta, MINIMAP_ALIGN_FLYE.out.paf)
                	LONGSTITCH (CUTADAPT.out.reads, RACON.out.improved_assembly, GOAT_TAXONSEARCH.out.genome_size)
	                // Gather versions of all tools used
	                ch_versions = ch_versions.mix(FLYE.out.versions)
        	} else if (params.assembly_secondary_mode== 'ont') {
                        mode = "--nano-raw"
                        FLYE (input_ont_fastq_1, mode)
			MINIMAP_ALIGN_FLYE (input_ont_fastq_1, FLYE.out.fasta.collect{it[1]}, false, false, false)
                        RACON (input_ont_fastq_1, FLYE.out.fasta, MINIMAP_ALIGN_FLYE.out.paf)
                        LONGSTITCH (input_ont_fastq_1, RACON.out.improved_assembly, GOAT_TAXONSEARCH.out.genome_size)
                        // Gather versions of all tools used
                        ch_versions = ch_versions.mix(FLYE.out.versions)
                } else if (params.assembly_secondary_mode== 'pacbio+ont') {
                        FLYE_PACBIO_ONT (CUTADAPT.out.reads, input_ont_fastq_1)
                        MINIMAP_ALIGN_FLYE (input_ont_fastq_1, FLYE_PACBIO_ONT.out.fasta.collect{it[1]}, false, false, false)
                        RACON (input_ont_fastq_1, FLYE_PACBIO_ONT.out.fasta, MINIMAP_ALIGN_FLYE.out.paf)
                        LONGSTITCH (input_ont_fastq_1, RACON.out.improved_assembly, GOAT_TAXONSEARCH.out.genome_size)
                        // Gather versions of all tools used
                        ch_versions = ch_versions.mix(FLYE_PACBIO_ONT.out.versions)
                } else if (params.assembly_secondary_mode== 'clr') {
                        mode = "--pacbio-raw"
                        FLYE (CUTADAPT.out.reads, mode)
                        MINIMAP_ALIGN_FLYE (CUTADAPT.out.reads, FLYE.out.fasta.collect{it[1]}, false, false, false)
                        RACON (CUTADAPT.out.reads, FLYE.out.fasta, MINIMAP_ALIGN_FLYE.out.paf)
                        LONGSTITCH (CUTADAPT.out.reads, RACON.out.improved_assembly, GOAT_TAXONSEARCH.out.genome_size)
                        // Gather versions of all tools used
                        ch_versions = ch_versions.mix(FLYE.out.versions)
                } else {
                        error "Invalid flye mode: params.assembly_secondary_mode. These modes are currently supported : 'hifi' or 'ont' or 'pacbio+ont' or 'clr'"
		}

		assembly_primary = LONGSTITCH.out.assembly
                // Gather versions of all tools used
                ch_versions = ch_versions.mix(MINIMAP_ALIGN_FLYE.out.versions)
                ch_versions = ch_versions.mix(RACON.out.versions)
                ch_versions = ch_versions.mix(LONGSTITCH.out.versions)
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
                // Gather versions of all tools used
                ch_versions = ch_versions.mix(VERKKO.out.versions)
	} else {
		error "Invalid alignment method: params.assembly_method. These methods are currently supported : 'hifiasm', 'canu', 'flye', 'verkko'. "
	}

	//QC post assembly
	if (params.lineage) {
		BUSCO_lin1_PRIM(assembly_primary, params.lineage, params.busco_lineages_path, [])
	} else {
		BUSCO_lin1_PRIM(assembly_primary, 'auto', [], [])
	}
	mqc_input = mqc_input.mix(BUSCO_lin1_PRIM.out.short_summaries_txt.collect{it[1]})
	// Gather versions of all tools used
	ch_versions = ch_versions.mix(BUSCO_lin1_PRIM.out.versions)
	if ((params.assembly_method == 'hifiasm') && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
	        QUAST_ASS_DOUBLE (assembly_primary, assembly_alternate, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
	        mqc_input = mqc_input.mix(QUAST_ASS_DOUBLE.out.tsv)
		quast_contig = QUAST_ASS_DOUBLE.out.renamed_tsv
		MERQURY_ASS_DOUBLE (MERYL_COUNT.out.meryl_db, assembly_primary, assembly_alternate)
       		// Gather versions of all tools used
	        ch_versions = ch_versions.mix(QUAST_ASS_DOUBLE.out.versions)
                ch_versions = ch_versions.mix(MERQURY_ASS_DOUBLE.out.versions)
	} else {
                QUAST_ASS (assembly_primary, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
                mqc_input = mqc_input.mix(QUAST_ASS.out.tsv)
                quast_contig = QUAST_ASS.out.renamed_tsv
		MERQURY_ASS (MERYL_COUNT.out.meryl_db.join(assembly_primary))
                // Gather versions of all tools used
                ch_versions = ch_versions.mix(QUAST_ASS.out.versions)
                ch_versions = ch_versions.mix(MERQURY_ASS.out.versions)
	}

	//Polishing (likely not going to happen, only for primary assembly for now)
	if (params.polishing_method == 'pilon') {
		BWAMEM2_INDEX_PILON(assembly_primary)
		BWAMEM2_MEM_PILON(input_illumina_SR_R1_R2, BWAMEM2_INDEX_PILON.out.index, true)
		SAMTOOLS_INDEX_PILON(BWAMEM2_MEM_PILON.out.bam)
		pilon_mode="--frags"
		PILON(assembly_primary, pilon_mode, BWAMEM2_MEM_PILON.out.bam.join(SAMTOOLS_INDEX_PILON.out.bai))
		assembly_polished = PILON.out.improved_assembly
		QUAST_PILON(assembly_polished, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
		mqc_input = mqc_input.mix(QUAST_PILON.out.tsv)
		assembly_unpurged = assembly_polished
                // Gather versions of all tools used
                ch_versions = ch_versions.mix(BWAMEM2_INDEX_PILON.out.versions)
                ch_versions = ch_versions.mix(BWAMEM2_MEM_PILON.out.versions)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PILON.out.versions)
                ch_versions = ch_versions.mix(PILON.out.versions)
	} else {
		assembly_unpurged = assembly_primary
	}

	if (params.fcs == 'yes') {		
		//Assembly cleaning
		FCS_FCSADAPTOR_hap1(assembly_unpurged)
		FCS_FCSGX_hap1(FCS_FCSADAPTOR_hap1.out.cleaned_assembly)
		FCS_FCSGX_CLEAN_hap1(FCS_FCSADAPTOR_hap1.out.cleaned_assembly, FCS_FCSGX_hap1.out.fcs_gx_report)
	        cleaned_hap1 = FCS_FCSGX_CLEAN_hap1.out.cleaned_fasta
		// Gather versions of all tools used
		ch_versions = ch_versions.mix(FCS_FCSADAPTOR_hap1.out.versions)
	        ch_versions = ch_versions.mix(FCS_FCSGX_hap1.out.versions)
	        ch_versions = ch_versions.mix(FCS_FCSGX_CLEAN_hap1.out.versions)
		if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
			FCS_FCSADAPTOR_ALT(assembly_alternate)
		        FCS_FCSGX_ALT(FCS_FCSADAPTOR_ALT.out.cleaned_assembly)
		        FCS_FCSGX_CLEAN_ALT(FCS_FCSADAPTOR_ALT.out.cleaned_assembly, FCS_FCSGX_ALT.out.fcs_gx_report)
               	        cleaned_hap2 = FCS_FCSGX_CLEAN_ALT.out.cleaned_fasta
		}
	} else {
		cleaned_hap1 = assembly_unpurged
		if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
			cleaned_hap2 = assembly_alternate
		}
	}

	//QC on cleaned contig assemblies
	if (params.lineage) {
		BUSCO_lin1_cleaned(cleaned_hap1, params.lineage, params.busco_lineages_path, [])
	} else {
                BUSCO_lin1_cleaned(cleaned_hap1, 'auto', [], [])
	}
        mqc_input = mqc_input.mix(BUSCO_lin1_cleaned.out.short_summaries_txt.collect{it[1]})
        if ((params.assembly_method == 'hifiasm') && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
                QUAST_CLEAN_DOUBLE (cleaned_hap1, cleaned_hap2, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
                mqc_input = mqc_input.mix(QUAST_CLEAN_DOUBLE.out.tsv)
        } else {
                QUAST_CLEAN (cleaned_hap1, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
                mqc_input = mqc_input.mix(QUAST_CLEAN.out.tsv)
        }

	//PurgeDups for primary assembly
        PURGEDUPS_SPLITFA_hap1 (cleaned_hap1)
        MINIMAP2_ALIGN_TO_CONTIG (CUTADAPT.out.reads, cleaned_hap1.collect{it[1]}, false, false, false)
        MINIMAP2_ALIGN_TO_SELF (PURGEDUPS_SPLITFA_hap1.out.split_fasta, [], false, false, false)
        PURGEDUPS_PBCSTAT_hap1 (MINIMAP2_ALIGN_TO_CONTIG.out.paf)
        PURGEDUPS_CALCUTS_hap1 (PURGEDUPS_PBCSTAT_hap1.out.stat) 
        PURGEDUPS_PURGEDUPS_hap1 (PURGEDUPS_PBCSTAT_hap1.out.basecov.join (PURGEDUPS_CALCUTS_hap1.out.cutoff), MINIMAP2_ALIGN_TO_SELF.out.paf )
        PURGEDUPS_GETSEQS_hap1 (cleaned_hap1, PURGEDUPS_PURGEDUPS_hap1.out.bed)
        SAMTOOLS_FAIDX1 (PURGEDUPS_GETSEQS_hap1.out.purged)
	purged_primary = PURGEDUPS_GETSEQS_hap1.out.purged
	if (params.lineage) {
		BUSCO_lin1_purged(purged_primary, params.lineage, params.busco_lineages_path, [])
	} else {
                BUSCO_lin1_purged(purged_primary, 'auto', [], [])
	}
        mqc_input = mqc_input.mix(BUSCO_lin1_purged.out.short_summaries_txt.collect{it[1]})

        // Gather versions of all tools used
        ch_versions = ch_versions.mix(PURGEDUPS_SPLITFA_hap1.out.versions)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_TO_CONTIG.out.versions)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_TO_SELF.out.versions)
        ch_versions = ch_versions.mix(PURGEDUPS_PBCSTAT_hap1.out.versions)
        ch_versions = ch_versions.mix(PURGEDUPS_CALCUTS_hap1.out.versions)
        ch_versions = ch_versions.mix(PURGEDUPS_PURGEDUPS_hap1.out.versions)
        ch_versions = ch_versions.mix(PURGEDUPS_GETSEQS_hap1.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX1.out.versions)

        if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
                //Merge haplotig from purge_dups and alternate assembly from hifiasm
                CAT (cleaned_hap2, PURGEDUPS_GETSEQS_hap1.out.haplotigs)
                PURGEDUPS_SPLITFA_ALT (CAT.out.alternate_contigs_full)
                MINIMAP2_ALIGN_TO_CONTIG_ALT (CUTADAPT.out.reads, CAT.out.alternate_contigs_full.collect{it[1]}, false, false, false)
                MINIMAP2_ALIGN_TO_SELF_ALT (PURGEDUPS_SPLITFA_ALT.out.split_fasta, [], false, false, false)
                PURGEDUPS_PBCSTAT_ALT (MINIMAP2_ALIGN_TO_CONTIG_ALT.out.paf)
                PURGEDUPS_CALCUTS_ALT (PURGEDUPS_PBCSTAT_ALT.out.stat)
                PURGEDUPS_PURGEDUPS_ALT (PURGEDUPS_PBCSTAT_ALT.out.basecov.join (PURGEDUPS_CALCUTS_ALT.out.cutoff), MINIMAP2_ALIGN_TO_SELF_ALT.out.paf )
                PURGEDUPS_GETSEQS_ALT (CAT.out.alternate_contigs_full, PURGEDUPS_PURGEDUPS_ALT.out.bed)
                SAMTOOLS_FAIDX1_ALT (PURGEDUPS_GETSEQS_ALT.out.purged)
                purged_alternate = PURGEDUPS_GETSEQS_ALT.out.purged

                QUAST_PURGED_DOUBLE (purged_primary, purged_alternate, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
                mqc_input = mqc_input.mix(QUAST_PURGED_DOUBLE.out.tsv)
		quast_contig_purged=QUAST_PURGED_DOUBLE.out.renamed_tsv
                MERQURY_PURGED_DOUBLE(MERYL_COUNT.out.meryl_db, purged_primary, purged_alternate)
        } else {
		QUAST_PURGED (purged_primary, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
		mqc_input = mqc_input.mix(QUAST_PURGED.out.tsv)
		quast_contig_purged = QUAST_PURGED.out.renamed_tsv
                MERQURY_PURGED(MERYL_COUNT.out.meryl_db.join(purged_primary))
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
        		CHROMAP_INDEX_hap1(PURGEDUPS_GETSEQS_hap1.out.purged)
        		CHROMAP_CHROMAP_hap1(input_hic_R1_R2, PURGEDUPS_GETSEQS_hap1.out.purged, CHROMAP_INDEX_hap1.out.index, [],[],[],[])
        		YAHS_hap1(PURGEDUPS_GETSEQS_hap1.out.purged, SAMTOOLS_FAIDX1.out.fai, CHROMAP_CHROMAP_hap1.out.bam)
        		SAMTOOLS_FAIDX2(YAHS_hap1.out.fasta)

		        // Gather versions of all tools used
		        ch_versions = ch_versions.mix(CHROMAP_INDEX_hap1.out.versions)
                        ch_versions = ch_versions.mix(CHROMAP_CHROMAP_hap1.out.versions)
                        ch_versions = ch_versions.mix(YAHS_hap1.out.versions)
                        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX2.out.versions)

        		scaffold 		= YAHS_hap1.out.fasta
			scaffold_agp 		= YAHS_hap1.out.agp
			scaffold_bin 		= YAHS_hap1.out.bin
        		scaffold_index		= SAMTOOLS_FAIDX2.out.fai

		        if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
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
		if (params.lineage) {
			BUSCO_lin1_SCAFF(scaffold, params.lineage, params.busco_lineages_path, [])
		} else {
                        BUSCO_lin1_SCAFF(scaffold, 'auto', [], [])
		}
        	mqc_input = mqc_input.mix(BUSCO_lin1_SCAFF.out.short_summaries_txt.collect{it[1]})
        	if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
                	QUAST_SCAFF_DOUBLE (scaffold, scaffold_alt, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
                	mqc_input = mqc_input.mix(QUAST_SCAFF_DOUBLE.out.tsv)
			quast_scaffold = QUAST_SCAFF_DOUBLE.out.renamed_tsv
        		MERQURY_SCAFF_DOUBLE(MERYL_COUNT.out.meryl_db, scaffold, scaffold_alt)
		} else {
	                QUAST_SCAFF (scaffold, quast_fasta, quast_gff, false, false, GOAT_TAXONSEARCH.out.genome_size)
	                mqc_input = mqc_input.mix(QUAST_SCAFF.out.tsv)
			quast_scaffold = QUAST_SCAFF.out.renamed_tsv
			MERQURY_SCAFF(MERYL_COUNT.out.meryl_db.join(scaffold))
	        }

		if (params.methylation_calling == 'yes') {
			//Map PacBio data against newly genevated assembly
		        JASMINE (bamtools_filter_output)
 	               ch_versions = ch_versions.mix(JASMINE.out.versions)
	                PBMM2 (JASMINE.out.cpg_bam, scaffold)
		} else {
			PBMM2 (bamtools_filter_output, scaffold)
		}
                SAMTOOLS_INDEX_PBMM2 (PBMM2.out.aligned_bam)
	        // Gather versions of all tools used
	        ch_versions = ch_versions.mix(PBMM2.out.versions)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PBMM2.out.versions)
         

                // JUICER must have contig fai for scaffold assembly
	        YAHS_JUICER (scaffold_agp, scaffold_bin, SAMTOOLS_FAIDX1.out.fai)
		chrom_size = YAHS_JUICER.out.chrom_sizes

// Gather versions of all tools used
		ch_versions = ch_versions.mix(YAHS_JUICER.out.versions)
		
                GFASTATS(scaffold)
                //Identify telomere sequences
                TIDK(scaffold)
                //Calculate Pacbio coverage and output a bedgraph
                BEDTOOLS_GENOMECOV(PBMM2.out.aligned_bam, '1', [], 'bedgraph')
                ch_versions = ch_versions.mix(GFASTATS.out.versions)
                ch_versions = ch_versions.mix(TIDK.out.versions)

		if (params.pretext == 'yes'){
			//PRETEXT
	                PRETEXTMAP(YAHS_JUICER.out.chrom_sizes, YAHS_JUICER.out.alignments_sorted_txt)
	                PRETEXTSNAPSHOT (PRETEXTMAP.out.pretext)
			//Add the telomere bedgraph to pretextgraph
		        PRETEXTGRAPH_TELO(PRETEXTMAP.out.pretext, TIDK.out.bedgraph_telomere, 'telomere')
	                //Add the coverage bedgraph to pretextgrap
		        PRETEXTGRAPH_TELO_COV(PRETEXTGRAPH_TELO.out.pretext, BEDTOOLS_GENOMECOV.out.genomecov, 'coverage')

	                ch_versions = ch_versions.mix(PRETEXTMAP.out.versions)
	                ch_versions = ch_versions.mix(PRETEXTSNAPSHOT.out.versions)
	                ch_versions = ch_versions.mix(PRETEXTGRAPH_TELO.out.versions)
		}

                if ((params.juicer == 'yes')) {
                        // JUICER must have contig fai for scaffold assembly
                        JUICER(YAHS_JUICER.out.chrom_sizes, YAHS_JUICER.out.alignments_sorted_txt)
                        // Gather versions of all tools used
                        ch_versions = ch_versions.mix(JUICER.out.versions)
                }


	        //Genome comparison
	        if (params.genome_comparison == 'yes') {
			input_ncbi = [ [ id:params.related_genome, single_end:true ] ]
			NCBIGENOMEDOWNLOAD(input_ncbi, [])
	                JUPITER(scaffold, NCBIGENOMEDOWNLOAD.out.fna)
			MASHMAP(scaffold, NCBIGENOMEDOWNLOAD.out.fna)
                        ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions)
                        ch_versions = ch_versions.mix(JUPITER.out.versions)
			ch_versions = ch_versions.mix(MASHMAP.out.versions)
	        }
		

		GZIP(scaffold)

		if (params.blobtools == 'yes'){
			BLOBTOOLS_CONFIG(GZIP.out.gz, bam2fastx_output)
	                BLOBTOOLS_PIPELINE(BLOBTOOLS_CONFIG.out.config, GZIP.out.gz)
	                BLOBTOOLS_CREATE(scaffold, BLOBTOOLS_CONFIG.out.config)
	                BLOBTOOLS_ADD(BLOBTOOLS_PIPELINE.out.blast_out, BLOBTOOLS_PIPELINE.out.diamond_proteome_out, BLOBTOOLS_PIPELINE.out.diamond_busco_out, BLOBTOOLS_PIPELINE.out.assembly_minimap_bam, BLOBTOOLS_PIPELINE.out.hic_minimap_bam , BLOBTOOLS_PIPELINE.out.lineage1_full_table_tsv , BLOBTOOLS_PIPELINE.out.lineage2_full_table_tsv, BLOBTOOLS_CREATE.out.blobtools_folder)
	                BLOBTOOLS_VIEW(BLOBTOOLS_ADD.out.blobtools_folder)

	                ch_versions = ch_versions.mix(BLOBTOOLS_PIPELINE.out.versions)
		}

		RAPIDCURATION_SPLIT(scaffold)
                ch_versions = ch_versions.mix(RAPIDCURATION_SPLIT.out.versions)
	} else {
		quast_scaffold = file('quast_scaffold_dummy')
                chrom_size = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('chrom_size_dummy')]
                ]
	}


	//Busco is ran on all lineages for hap 1 only for the most final assembly (scaffold > purged > contig)
	//BUSCO is ran on the principal lineage for hap2 only for the most final assembly (scaffold > purged > contig)
	//Define which file to run Busco on
	if (( params.hic_read1 ) && ( params.hic_read2 )) {
		busco_assembly = scaffold
		if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
			busco_assembly_alt = scaffold_alt
		}
	} else {
		busco_assembly = purged_primary
		if (( params.assembly_method == 'hifiasm' )  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
			busco_assembly_alt = purged_alternate
		}
	}

       if (( params.hic_read1 ) && ( params.hic_read2 )) {
		busco_lin1_json=BUSCO_lin1_SCAFF.out.short_summaries_json
	} else {
		busco_lin1_json=BUSCO_lin1_purged.out.short_summaries_json
	}

/*
        if ((params.assembly_method == 'hifiasm')  && (GOAT_TAXONSEARCH.out.ploidy != '1')) {
		if (params.lineage) {		
	                BUSCO_ALT (busco_assembly_alt, params.lineage, params.busco_lineages_path, [])
		} else {
                        BUSCO_ALT (busco_assembly_alt, 'auto', [], [])
		}
                mqc_input = mqc_input.mix(BUSCO_ALT.out.short_summaries_txt.collect{it[1]})
        }
*/

        if (params.lineage2) {
		BUSCO_lin2(busco_assembly, params.lineage2, params.busco_lineages_path, [])         
		mqc_input = mqc_input.mix(BUSCO_lin2.out.short_summaries_txt.collect{it[1]})
                busco_lin2_json = BUSCO_lin2.out.short_summaries_json
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
        } else {
                busco_lin4_json = [
                        [ id:'dummy', single_end: true], // meta map
                        [ file('busco_lin4_json_dummy')]
                ]
	}

	// Gather versions of all tools used
	ch_version_yaml = Channel.empty()
	CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
	ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

	//MultiQC report
        mqc_input = mqc_input.mix(ch_version_yaml)
	MULTIQC (mqc_input.collect(), [], [], [], COVERAGE_CALCULATION.out.coverage)
	ch_versions = ch_versions.mix(MULTIQC.out.versions)

//	OVERVIEW_GENERATION_SAMPLE(LONGQC_PACBIO.out.report_json, kraken_pacbio, kraken_hic, quast_contig, quast_contig_purged, quast_scaffold, busco_lin1_json, busco_lin2_json, busco_lin3_json, busco_lin4_json, chrom_size, GOAT_TAXONSEARCH.out.ploidy, GOAT_TAXONSEARCH.out.haploid_number, GOAT_TAXONSEARCH.out.scientific_name, GOAT_TAXONSEARCH.out.genome_size)


}


