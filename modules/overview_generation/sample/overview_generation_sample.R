#rjson_url = "https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.20.tar.gz"

#install.packages(rjson_url, repos=NULL, type="source", lib="${params.singularity_cache}")
library("rjson")
#library(dplyr)
#library(stringr)


args <- commandArgs(trailingOnly = TRUE)

xxid =(args[1])
pipeline_version = (args[2])
outdir =(args[3])
pacbio_input_type =(args[4])
bam_cell1 =(args[5])
bam_cell2 =(args[6])
bam_cell3 =(args[7])
bam_cell4 =(args[8])
ont_fastq_1 =(args[9])
hic_read1 =(args[10])
hic_read2 =(args[11])
illumina_SR_read1 =(args[12])
illumina_SR_read2 =(args[13])
pacbio_rq =(args[14])
assembly_method =(args[15])
assembly_secondary_mode =(args[16])
polishing_method =(args[17])
purging_method = (args[18])
scaffolding_method =(args[19])
manual_curation=(args[20])
mitohifi=(args[21])
taxon_taxid =(args[22])
ploidy =(args[23])
hap_gen_size_Gb =(args[24])
chrom_num =(args[25])
lineage =(args[26])
lineage2 =(args[27])
lineage3 =(args[28])
lineage4 =(args[29])

taxon_name =paste0((args[41]), "_", (args[42]))

#Depending on the version, the genome size is in bp or in Gb, need to make it all in Gb
if (hap_gen_size_Gb >  1000) {
	hap_gen_size_Gb = as.numeric(hap_gen_size_Gb)/1000000000
}

#Count the number of files
#number of pacbio bam files (including both barcoded and unbarcoded)
pacbio_concat = paste(bam_cell1, bam_cell2, bam_cell3, bam_cell4, sep = "; ")
#pacbio_n_files = 4-str_count(pacbio_concat, "null")
pacbio_n_files = 4 - lengths(regmatches(pacbio_concat, gregexpr("null", pacbio_concat)))


#hic
#Number of PAIRED filed
hic_concat = paste(hic_read1, hic_read2, sep = "; ")
#hic_n_paired_files = 2-str_count(hic_concat, "null")
hic_n_paired_files = (2 - lengths(regmatches(hic_concat, gregexpr("null", hic_concat))))/2

#pe150 (illumina short reads)
#Number of PAIRED filed
pe150_concat = paste(illumina_SR_read1, illumina_SR_read2, sep = "; ")
#pe150_n_paired_files = 2-str_count(pe150_concat, "null")
pe150_n_paired_files = (2 - lengths(regmatches(pe150_concat, gregexpr("null", pe150_concat))))/2

#ONT
ont_concat = paste(ont_fastq_1, sep = "; ")
#ont_n_files = 1-str_count(ont_concat, "null")
ont_n_files = 1- lengths(regmatches(ont_concat, gregexpr("null", ont_concat)))

#Extract the information from the different files to generate aggregated tsv and figures

#LongQC data
if (file.exists(args[30])) {
	#longqc_report= fromJSON(file="~/Downloads/Greenland_cockle_R/QC_vals_longQC_sampleqc_Greenland_cockle_004.json")
	longqc_report= fromJSON(file=(args[30]))
	pacbio_n_reads_longqc = longqc_report[["Num_of_reads"]]
	pacbio_longest_read_longqc_kb = longqc_report[["Longest_read"]]/1000
	pacbio_mean_read_length_longqc_kb = longqc_report[["Length_stats"]][["Mean_read_length"]]/1000
	pacbio_coverage_x=(pacbio_n_reads_longqc*pacbio_mean_read_length_longqc_kb)/(as.numeric(hap_gen_size_Gb)*1000000)	
} else {
	pacbio_n_reads_longqc = "NA"
	pacbio_longest_read_longqc_kb = "NA"
	pacbio_mean_read_length_longqc_kb = "NA"
	pacbio_coverage_xpacbio_coverage_xpacbio_coverage_x = "NA"
	pacbio_coverage_x = "NA"
}
lonqc_overview = cbind(pacbio_n_reads_longqc, pacbio_longest_read_longqc_kb, pacbio_mean_read_length_longqc_kb, pacbio_coverage_x)
#head(lonqc_overview)


#Pacbio Kraken results
#kraken_pacbio_results = read.table("~/Downloads/Greenland_cockle_R/Greenland_cockle_004.kraken2.report.txt", sep = "\t", blank.lines.skip = FALSE, quote="", comment.char="")
kraken_pacbio_results = read.table(args[31], sep = "\t", blank.lines.skip = FALSE, quote="", comment.char="")
colnames(kraken_pacbio_results)=c("perc", "n_reads_clade", "n_reads_taxon", "rank_code", "NCBI_tax", "scientific_name")
#Remove spaces in scientific_name column
kraken_pacbio_results$scientific_name = gsub('\\s+', '', kraken_pacbio_results$scientific_name)
kraken_pacbio_actinopteri_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Actinopteri", 1]
kraken_pacbio_amphibia_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Amphibia", 1]
kraken_pacbio_aves_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Aves", 1]
kraken_pacbio_bivalvia_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Bivalvia", 1]
kraken_pacbio_insecta_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Insecta", 1]
kraken_pacbio_magnoliopsida_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Magnoliopsida", 1]
kraken_pacbio_mammallia_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "Mammalia", 1]
kraken_pacbio_unclassified_perc = kraken_pacbio_results[kraken_pacbio_results$scientific_name == "unclassified", 1]
kraken_pacbio_other_perc = 100 - (kraken_pacbio_actinopteri_perc + kraken_pacbio_amphibia_perc + kraken_pacbio_aves_perc + kraken_pacbio_bivalvia_perc + kraken_pacbio_insecta_perc + kraken_pacbio_magnoliopsida_perc + kraken_pacbio_mammallia_perc + kraken_pacbio_unclassified_perc)

kraken_pacbio_overview = cbind(kraken_pacbio_actinopteri_perc, kraken_pacbio_amphibia_perc, kraken_pacbio_aves_perc, kraken_pacbio_bivalvia_perc, kraken_pacbio_insecta_perc, kraken_pacbio_magnoliopsida_perc, kraken_pacbio_mammallia_perc, kraken_pacbio_unclassified_perc, kraken_pacbio_other_perc)

#Kraken hic results
if (hic_n_paired_files > 0){
	#kraken_hic_results = read.table("~/Downloads/Greenland_cockle_R/Greenland_cockle_004.kraken2.report.txt", sep = "\t", blank.lines.skip = FALSE, quote="", comment.char="")
	kraken_hic_results = read.table(args[32], sep = "\t", blank.lines.skip = FALSE, quote="", comment.char="")
	colnames(kraken_hic_results)=c("perc", "n_reads_clade", "n_reads_taxon", "rank_code", "NCBI_tax", "scientific_name")
	#Remove spaces in scientific_name column
	kraken_hic_results$scientific_name = gsub('\\s+', '', kraken_hic_results$scientific_name)
	kraken_hic_actinopteri_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Actinopteri", 1]
	kraken_hic_amphibia_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Amphibia", 1]
	kraken_hic_aves_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Aves", 1]
	kraken_hic_bivalvia_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Bivalvia", 1]
	kraken_hic_insecta_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Insecta", 1]
	kraken_hic_magnoliopsida_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Magnoliopsida", 1]
	kraken_hic_mammallia_perc = kraken_hic_results[kraken_hic_results$scientific_name == "Mammalia", 1]
	kraken_hic_unclassified_perc = kraken_hic_results[kraken_hic_results$scientific_name == "unclassified", 1]
	kraken_hic_other_perc = 100 - (kraken_hic_actinopteri_perc + kraken_hic_amphibia_perc + kraken_hic_aves_perc + kraken_hic_bivalvia_perc + kraken_hic_insecta_perc + kraken_hic_magnoliopsida_perc + kraken_hic_mammallia_perc + kraken_hic_unclassified_perc)
	kraken_hic_root_numer_reads = kraken_hic_results[kraken_hic_results$scientific_name == "root", 2]
        kraken_hic_unclassified_numer_reads = kraken_hic_results[kraken_hic_results$scientific_name == "unclassified", 2]
	kraken_hic_root_unclassified_numer_reads = kraken_hic_root_numer_reads + kraken_hic_unclassified_numer_reads
} else {
	kraken_hic_actinopteri_perc = "NA"
	kraken_hic_amphibia_perc = "NA"
	kraken_hic_aves_perc = "NA"
	kraken_hic_bivalvia_perc = "NA"
	kraken_hic_insecta_perc = "NA"
	kraken_hic_magnoliopsida_perc = "NA"
	kraken_hic_mammallia_perc = "NA"
	kraken_hic_unclassified_perc = "NA"
	kraken_hic_other_perc = "NA"
	kraken_hic_root_numer_reads = "NA"
	kraken_hic_unclassified_numer_reads = "NA"
	kraken_hic_root_unclassified_numer_reads = "NA"
}

kraken_hic_overview = cbind(kraken_hic_actinopteri_perc, kraken_hic_amphibia_perc, kraken_hic_aves_perc, kraken_hic_bivalvia_perc, kraken_hic_insecta_perc, kraken_hic_magnoliopsida_perc, kraken_hic_mammallia_perc, kraken_hic_unclassified_perc, kraken_hic_other_perc, kraken_hic_root_numer_reads, kraken_hic_unclassified_numer_reads, kraken_hic_root_unclassified_numer_reads)

#Quast
#Differs if hifiasm is used as there is both haplotypes (Quast_double)
#from quast_report.tsv
#quast_table_contig = read.table("~/Downloads/Greenland_cockle_R/Greenland_cockle_004_report.tsv", sep = "\t", header=T, comment.char="@")
quast_table_contig = read.table(args[33], sep = "\t", header=T, comment.char="@")

hap1_contig_n50_quast_mb = quast_table_contig[quast_table_contig$Assembly=="N50",2]/1000000
hap1_contig_l50_quast = quast_table_contig[quast_table_contig$Assembly=="L50",2]
hap1_contig_n90_quast_mb = quast_table_contig[quast_table_contig$Assembly=="N90",2]/1000000
hap1_contig_l90_quast = quast_table_contig[quast_table_contig$Assembly=="L90",2]
hap1_contig_assembly_length_quast_gb = quast_table_contig[quast_table_contig$Assembly =="Total length",2]/1000000000
hap1_contig_largest_contig_quast_mb = quast_table_contig[quast_table_contig$Assembly =="Largest contig",2]/1000000
hap1_contig_number_quast = quast_table_contig[quast_table_contig$Assembly =="# contigs",2]
hap1_contig_GC_quast = quast_table_contig[quast_table_contig$Assembly == "GC (%)", 2]

if (assembly_method == "hifiasm") {
	hap2_contig_n50_quast_mb = quast_table_contig[quast_table_contig$Assembly=="N50",3]/1000000
	hap2_contig_l50_quast = quast_table_contig[quast_table_contig$Assembly=="L50",3]
	hap2_contig_n90_quast_mb = quast_table_contig[quast_table_contig$Assembly=="N90",3]/1000000
	hap2_contig_l90_quast = quast_table_contig[quast_table_contig$Assembly=="L90",3]
	hap2_contig_assembly_length_quast_gb = quast_table_contig[quast_table_contig$Assembly =="Total length",3]/1000000000
	hap2_contig_largest_contig_quast_mb = quast_table_contig[quast_table_contig$Assembly =="Largest contig",3]/1000000
	hap2_contig_number_quast = quast_table_contig[quast_table_contig$Assembly =="# contigs",3]
	hap2_contig_GC_quast = quast_table_contig[quast_table_contig$Assembly == "GC (%)", 3]
} else {
	hap2_contig_n50_quast_mb = "NA"
	hap2_contig_l50_quast = "NA"
	hap2_contig_n90_quast_mb = "NA"
	hap2_contig_l90_quast = "NA"
	hap2_contig_assembly_length_quast_gb = "NA"
	hap2_contig_largest_contig_quast_mb = "NA"
	hap2_contig_number_quast = "NA"
	hap2_contig_GC_quast = "NA"
}

quast_contig_overview = cbind (hap1_contig_n50_quast_mb, hap1_contig_l50_quast, hap1_contig_n90_quast_mb, hap1_contig_l90_quast, hap1_contig_assembly_length_quast_gb, hap1_contig_largest_contig_quast_mb, hap1_contig_number_quast, hap1_contig_GC_quast, hap2_contig_n50_quast_mb, hap2_contig_l50_quast, hap2_contig_n90_quast_mb, hap2_contig_l90_quast, hap2_contig_assembly_length_quast_gb, hap2_contig_largest_contig_quast_mb, hap2_contig_number_quast, hap2_contig_GC_quast)

#If assembly is purged
if (file.exists(args[34])){
	#quast_table_contig_purged = read.table("~/Downloads/Greenland_cockle_R/Greenland_cockle_004_report.tsv", sep = "\t", header=T, comment.char="@")
	quast_table_contig_purged = read.table(args[34], sep = "\t", header=T, comment.char="@")
	
	hap1_contig_purged_n50_quast_mb = quast_table_contig_purged[quast_table_contig_purged$Assembly=="N50",2]/1000000
	hap1_contig_purged_l50_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly=="L50",2]
	hap1_contig_purged_n90_quast_mb = quast_table_contig_purged[quast_table_contig_purged$Assembly=="N90",2]/1000000
	hap1_contig_purged_l90_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly=="L90",2]
	hap1_contig_purged_assembly_length_quast_gb = quast_table_contig_purged[quast_table_contig_purged$Assembly =="Total length",2]/1000000000
	hap1_contig_purged_largest_contig_quast_mb = quast_table_contig_purged[quast_table_contig_purged$Assembly =="Largest contig",2]/1000000
	hap1_contig_purged_number_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly =="# contigs",2]
        hap1_contig_purged_GC_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly =="GC (%)", 2]
} else {
	hap1_contig_purged_n50_quast_mb = "NA"
	hap1_contig_purged_l50_quast = "NA"
	hap1_contig_purged_n90_quast_mb = "NA"
	hap1_contig_purged_l90_quast = "NA"
	hap1_contig_purged_assembly_length_quast_gb = "NA"
	hap1_contig_purged_largest_contig_quast_mb = "NA"
	hap1_contig_purged_number_quast = "NA"
        hap1_contig_purged_GC_quast = "NA"
}

if ((file.exists(args[34])) & assembly_method == "hifiasm") {
	#quast_table_contig_purged = read.table("~/Downloads/Greenland_cockle_R/Greenland_cockle_004_report.tsv", sep = "\t", header=T, comment.char="@")
	quast_table_contig_purged = read.table(args[34], sep = "\t", header=T, comment.char="@")

	hap2_contig_purged_n50_quast_mb = quast_table_contig_purged[quast_table_contig_purged$Assembly=="N50",3]/1000000
	hap2_contig_purged_l50_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly=="L50",3]
	hap2_contig_purged_n90_quast_mb = quast_table_contig_purged[quast_table_contig_purged$Assembly=="N90",3]/1000000
	hap2_contig_purged_l90_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly=="L90",3]
	hap2_contig_purged_assembly_length_quast_gb = quast_table_contig_purged[quast_table_contig_purged$Assembly =="Total length",3]/1000000000
	hap2_contig_purged_largest_contig_quast_mb = quast_table_contig_purged[quast_table_contig_purged$Assembly =="Largest contig",3]/1000000
	hap2_contig_purged_number_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly =="# contigs",3]
        hap2_contig_purged_GC_quast = quast_table_contig_purged[quast_table_contig_purged$Assembly =="GC (%)", 3]
} else {
	hap2_contig_purged_n50_quast_mb = "NA"
	hap2_contig_purged_l50_quast = "NA"
	hap2_contig_purged_n90_quast_mb = "NA"
	hap2_contig_purged_l90_quast = "NA"
	hap2_contig_purged_assembly_length_quast_gb = "NA"
	hap2_contig_purged_largest_contig_quast_mb = "NA"
	hap2_contig_purged_number_quast = "NA"
        hap2_contig_purged_GC_quast = "NA"
}

quast_contig_purged_overview = cbind (hap1_contig_purged_n50_quast_mb, hap1_contig_purged_l50_quast, hap1_contig_purged_n90_quast_mb, hap1_contig_purged_l90_quast, hap1_contig_purged_assembly_length_quast_gb, hap1_contig_purged_largest_contig_quast_mb, hap1_contig_purged_number_quast, hap1_contig_purged_GC_quast, hap2_contig_purged_n50_quast_mb, hap2_contig_purged_l50_quast, hap2_contig_purged_n90_quast_mb, hap2_contig_purged_l90_quast, hap2_contig_purged_assembly_length_quast_gb, hap2_contig_purged_largest_contig_quast_mb, hap2_contig_purged_number_quast, hap2_contig_purged_GC_quast)


#If assembly is scaffolded
if (file.exists(args[35])){
	#quast_table_scaffold = read.table("~/Downloads/Greenland_cockle_R/Greenland_cockle_004_report.tsv", sep = "\t", header=T, comment.char="@")
	quast_table_scaffold = read.table(args[35], sep = "\t", header=T, comment.char="@")
	
	hap1_scaffold_n50_quast_mb = quast_table_scaffold[quast_table_scaffold$Assembly=="N50",2]/1000000
	hap1_scaffold_l50_quast = quast_table_scaffold[quast_table_scaffold$Assembly=="L50",2]
	hap1_scaffold_n90_quast_mb = quast_table_scaffold[quast_table_scaffold$Assembly=="N90",2]/1000000
	hap1_scaffold_l90_quast = quast_table_scaffold[quast_table_scaffold$Assembly=="L90",2]
	hap1_scaffold_assembly_length_quast_gb = quast_table_scaffold[quast_table_scaffold$Assembly =="Total length",2]/1000000000
	hap1_scaffold_largest_contig_quast_mb = quast_table_scaffold[quast_table_scaffold$Assembly =="Largest contig",2]/1000000
	hap1_scaffold_number_quast = quast_table_scaffold[quast_table_scaffold$Assembly =="# contigs",2]
        hap1_scaffold_GC_quast = quast_table_scaffold[quast_table_scaffold$Assembly =="GC (%)", 2]
} else {
	hap1_scaffold_n50_quast_mb = "NA"
	hap1_scaffold_l50_quast = "NA"
	hap1_scaffold_n90_quast_mb = "NA"
	hap1_scaffold_l90_quast = "NA"
	hap1_scaffold_assembly_length_quast_gb = "NA"
	hap1_scaffold_largest_contig_quast_mb = "NA"
	hap1_scaffold_number_quast = "NA"
	hap1_scaffold_GC_quast = "NA"
}

if ((file.exists(args[35])) & assembly_method == "hifiasm") {
	quast_table_scaffold = read.table(args[35], sep = "\t", header=T, comment.char="@")
	hap2_scaffold_n50_quast_mb = quast_table_scaffold[quast_table_scaffold$Assembly=="N50",3]/1000000
	hap2_scaffold_l50_quast = quast_table_scaffold[quast_table_scaffold$Assembly=="L50",3]
	hap2_scaffold_n90_quast_mb = quast_table_scaffold[quast_table_scaffold$Assembly=="N90",3]/1000000
	hap2_scaffold_l90_quast = quast_table_scaffold[quast_table_scaffold$Assembly=="L90",3]
	hap2_scaffold_assembly_length_quast_gb = quast_table_scaffold[quast_table_scaffold$Assembly =="Total length",3]/1000000000
	hap2_scaffold_largest_contig_quast_mb = quast_table_scaffold[quast_table_scaffold$Assembly =="Largest contig",3]/1000000
	hap2_scaffold_number_quast = quast_table_scaffold[quast_table_scaffold$Assembly =="# contigs",3]
        hap2_scaffold_GC_quast = quast_table_scaffold[quast_table_scaffold$Assembly =="GC (%)", 3]
} else {
	hap2_scaffold_n50_quast_mb = "NA"
	hap2_scaffold_l50_quast = "NA"
	hap2_scaffold_n90_quast_mb = "NA"
	hap2_scaffold_l90_quast = "NA"
	hap2_scaffold_assembly_length_quast_gb = "NA"	
	hap2_scaffold_largest_contig_quast_mb = "NA"
	hap2_scaffold_number_quast = "NA"
	hap2_scaffold_GC_quast = "NA"
}

quast_scaffold_overview = cbind (hap1_scaffold_n50_quast_mb, hap1_scaffold_l50_quast, hap1_scaffold_n90_quast_mb, hap1_scaffold_l90_quast, hap1_scaffold_assembly_length_quast_gb, hap1_scaffold_largest_contig_quast_mb, hap1_scaffold_number_quast, hap1_scaffold_GC_quast, hap2_scaffold_n50_quast_mb, hap2_scaffold_l50_quast, hap2_scaffold_n90_quast_mb, hap2_scaffold_l90_quast, hap2_scaffold_assembly_length_quast_gb, hap2_scaffold_largest_contig_quast_mb, hap2_scaffold_number_quast, hap2_scaffold_GC_quast)


#Busco scores
#Busco are generated only on the latest file (scaffold>purged>contig) as it takes a bunch of time
#lineage1
if (file.exists(args[36])) {
#busco_report_lin1= fromJSON(file="~/Downloads/Greenland_cockle_R/short_summary.specific.vertebrata_odb10.Greenland_cockle_004.asm.bp.hap1.p_ctg.purged_scaffolds_final.fa.json")
	busco_report_lin1= fromJSON(file=args[36])
#	lin1 = busco_report_lin1[["lineage_dataset"]][["name"]]
#	busco_complete_single_lin1 = busco_report_lin1[["results"]][["Single copy"]]
#	busco_complete_duplicated_lin1 = busco_report_lin1[["results"]][["Multi copy"]]
#	busco_complete_single_duplicated_lin1 = busco_report_lin1[["results"]][["Complete"]]
#	busco_fragmented_lin1 = busco_report_lin1[["results"]][["Fragmented"]]
#	busco_missing_lin1 = busco_report_lin1[["results"]][["Missing"]]
        lin1_temp = busco_report_lin1[["dataset"]]
        lin1=sub('.*/', '', lin1_temp)
        busco_total_num_lin1 = as.numeric(busco_report_lin1[["dataset_total_buscos"]])
        busco_complete_single_lin1_num = busco_report_lin1[["S"]]
        busco_complete_single_lin1 = busco_complete_single_lin1_num*100/busco_total_num_lin1
        busco_complete_duplicated_lin1_num = busco_report_lin1[["D"]]
        busco_complete_duplicated_lin1 = busco_complete_duplicated_lin1_num*100/busco_total_num_lin1
        busco_complete_single_duplicated_lin1_num = busco_report_lin1[["C"]]
        busco_complete_single_duplicated_lin1 = busco_complete_single_duplicated_lin1_num*100/busco_total_num_lin1
        busco_fragmented_lin1_num = busco_report_lin1[["F"]]
        busco_fragmented_lin1 = busco_fragmented_lin1_num*100/busco_total_num_lin1
        busco_missing_lin1_num = busco_report_lin1[["M"]]
        busco_missing_lin1 = busco_missing_lin1_num*100/busco_total_num_lin1
} else {
	busco_report_lin1 = "NA"
	lin1  = "NA"
	busco_complete_single_lin1 = "NA"
	busco_complete_duplicated_lin1 = "NA"
	busco_complete_single_duplicated_lin1 = "NA"
	busco_fragmented_lin1 = "NA"
	busco_missing_lin1 = "NA"
}
busco_lin1_overview = cbind(lin1, busco_complete_single_lin1, busco_complete_duplicated_lin1, busco_complete_single_duplicated_lin1, busco_fragmented_lin1, busco_missing_lin1)
#colnames(busco_lin1_overview) = c(paste0("busco_complete_single_",lin1), paste0("busco_complete_duplicated_",lin1), paste0("busco_complete_single_duplicated_single_",lin1), paste0("busco_fragmented_",lin1), paste0("busco_missing_",lin1))

#lineage2
#if ( lineage2 != "null") {
if (file.exists(args[37])) {
	#busco_report_lin2= fromJSON(file="~/Downloads/Greenland_cockle_R/short_summary.specific.vertebrata_odb10.Greenland_cockle_004.asm.bp.hap1.p_ctg.purged_scaffolds_final.fa.json")
	busco_report_lin2= fromJSON(file=args[37])
        lin2_temp = busco_report_lin2[["dataset"]]
        lin2=sub('.*/', '', lin2_temp)
        busco_total_num_lin2 = as.numeric(busco_report_lin2[["dataset_total_buscos"]])
        busco_complete_single_lin2_num = busco_report_lin2[["S"]]
        busco_complete_single_lin2 = busco_complete_single_lin2_num*100/busco_total_num_lin2
        busco_complete_duplicated_lin2_num = busco_report_lin2[["D"]]
        busco_complete_duplicated_lin2 = busco_complete_duplicated_lin2_num*100/busco_total_num_lin2
        busco_complete_single_duplicated_lin2_num = busco_report_lin2[["C"]]
        busco_complete_single_duplicated_lin2 = busco_complete_single_duplicated_lin2_num*100/busco_total_num_lin2
        busco_fragmented_lin2_num = busco_report_lin2[["F"]]
        busco_fragmented_lin2 = busco_fragmented_lin2_num*100/busco_total_num_lin2
        busco_missing_lin2_num = busco_report_lin2[["M"]]
        busco_missing_lin2 = busco_missing_lin2_num*100/busco_total_num_lin2
#	lin2 = busco_report_lin2[["lineage_dataset"]][["name"]]
#	busco_complete_single_lin2 = busco_report_lin2[["results"]][["Single copy"]]
#	busco_complete_duplicated_lin2 = busco_report_lin2[["results"]][["Multi copy"]]
#	busco_complete_single_duplicated_lin2 = busco_report_lin2[["results"]][["Complete"]]
#	busco_fragmented_lin2 = busco_report_lin2[["results"]][["Fragmented"]]
#	busco_missing_lin2 = busco_report_lin2[["results"]][["Missing"]]
} else {
	lin2 = "null"
	busco_complete_single_lin2 = "NA"
	busco_complete_duplicated_lin2 = "NA"
	busco_complete_single_duplicated_lin2 = "NA"
	busco_fragmented_lin2 = "NA"
	busco_missing_lin2 = "NA"
}
	
busco_lin2_overview = cbind(lin2, busco_complete_single_lin2, busco_complete_duplicated_lin2, busco_complete_single_duplicated_lin2, busco_fragmented_lin2, busco_missing_lin2)
#colnames(busco_lin2_overview) = c(paste0("busco_complete_single_",lin2), paste0("busco_complete_duplicated_",lin2), paste0("busco_complete_single_duplicated_single_",lin2), paste0("busco_fragmented_",lin2), paste0("busco_missing_",lin2))

#lineage3
#if ( lineage3 != "null") {
if (file.exists(args[38])) {
	#busco_report_lin3= fromJSON(file="~/Downloads/Greenland_cockle_R/short_summary.specific.vertebrata_odb10.Greenland_cockle_004.asm.bp.hap1.p_ctg.purged_scaffolds_final.fa.json")
	busco_report_lin3= fromJSON(file=args[38])
        lin3_temp = busco_report_lin3[["dataset"]]
        lin3=sub('.*/', '', lin3_temp)
        busco_total_num_lin3 = as.numeric(busco_report_lin3[["dataset_total_buscos"]])
        busco_complete_single_lin3_num = busco_report_lin3[["S"]]
        busco_complete_single_lin3 = busco_complete_single_lin3_num*100/busco_total_num_lin3
        busco_complete_duplicated_lin3_num = busco_report_lin3[["D"]]
        busco_complete_duplicated_lin3 = busco_complete_duplicated_lin3_num*100/busco_total_num_lin3
        busco_complete_single_duplicated_lin3_num = busco_report_lin3[["C"]]
        busco_complete_single_duplicated_lin3 = busco_complete_single_duplicated_lin3_num*100/busco_total_num_lin3
        busco_fragmented_lin3_num = busco_report_lin3[["F"]]
        busco_fragmented_lin3 = busco_fragmented_lin3_num*100/busco_total_num_lin3
        busco_missing_lin3_num = busco_report_lin3[["M"]]
        busco_missing_lin3 = busco_missing_lin3_num*100/busco_total_num_lin3
#	lin3 = busco_report_lin3[["lineage_dataset"]][["name"]]
#	busco_complete_single_lin3 = busco_report_lin3[["results"]][["Single copy"]]
#	busco_complete_duplicated_lin3 = busco_report_lin3[["results"]][["Multi copy"]]
#	busco_complete_single_duplicated_lin3 = busco_report_lin3[["results"]][["Complete"]]
#	busco_fragmented_lin3 = busco_report_lin3[["results"]][["Fragmented"]]
#	busco_missing_lin3 = busco_report_lin3[["results"]][["Missing"]]
} else {
	lin3 = "null"
	busco_complete_single_lin3 = "NA"
	busco_complete_duplicated_lin3 = "NA"
	busco_complete_single_duplicated_lin3 = "NA"
	busco_fragmented_lin3 = "NA"
	busco_missing_lin3 = "NA"
}

busco_lin3_overview = cbind(lin3, busco_complete_single_lin3, busco_complete_duplicated_lin3, busco_complete_single_duplicated_lin3, busco_fragmented_lin3, busco_missing_lin3)
#colnames(busco_lin3_overview) = c(paste0("busco_complete_single_",lin3), paste0("busco_complete_duplicated_",lin3), paste0("busco_complete_single_duplicated_single_",lin3), paste0("busco_fragmented_",lin3), paste0("busco_missing_",lin3))

#lineage4
#if ( lineage4 != "null") {
if (file.exists(args[39])) {
	#busco_report_lin4= fromJSON(file="~/Downloads/Greenland_cockle_R/short_summary.specific.vertebrata_odb10.Greenland_cockle_004.asm.bp.hap1.p_ctg.purged_scaffolds_final.fa.json")
	busco_report_lin4= fromJSON(file=args[39])
        lin4_temp = busco_report_lin4[["dataset"]]
        lin4=sub('.*/', '', lin4_temp)
        busco_total_num_lin4 = as.numeric(busco_report_lin4[["dataset_total_buscos"]])
        busco_complete_single_lin4_num = busco_report_lin4[["S"]]
        busco_complete_single_lin4 = busco_complete_single_lin4_num*100/busco_total_num_lin4
        busco_complete_duplicated_lin4_num = busco_report_lin4[["D"]]
        busco_complete_duplicated_lin4 = busco_complete_duplicated_lin4_num*100/busco_total_num_lin4
        busco_complete_single_duplicated_lin4_num = busco_report_lin4[["C"]]
        busco_complete_single_duplicated_lin4 = busco_complete_single_duplicated_lin4_num*100/busco_total_num_lin4
        busco_fragmented_lin4_num = busco_report_lin4[["F"]]
        busco_fragmented_lin4 = busco_fragmented_lin4_num*100/busco_total_num_lin4
        busco_missing_lin4_num = busco_report_lin4[["M"]]
        busco_missing_lin4 = busco_missing_lin4_num*100/busco_total_num_lin4
#	lin4 = busco_report_lin4[["lineage_dataset"]][["name"]]
#	busco_complete_single_lin4 = busco_report_lin4[["results"]][["Single copy"]]
#	busco_complete_duplicated_lin4 = busco_report_lin4[["results"]][["Multi copy"]]
#	busco_complete_single_duplicated_lin4 = busco_report_lin4[["results"]][["Complete"]]
#	busco_fragmented_lin4 = busco_report_lin4[["results"]][["Fragmented"]]
#	busco_missing_lin4 = busco_report_lin4[["results"]][["Missing"]]
} else {
	lin4 = "null"
	busco_complete_single_lin4 = "NA"
	busco_complete_duplicated_lin4 = "NA"
	busco_complete_single_duplicated_lin4 = "NA"
	busco_fragmented_lin4 = "NA"
	busco_missing_lin4 = "NA"
}

busco_lin4_overview = cbind(lin4, busco_complete_single_lin4, busco_complete_duplicated_lin4, busco_complete_single_duplicated_lin4, busco_fragmented_lin4, busco_missing_lin4)
#colnames(busco_lin4_overview) = c(paste0("busco_complete_single_",lin4), paste0("busco_complete_duplicated_",lin4), paste0("busco_complete_single_duplicated_single_",lin4), paste0("busco_fragmented_",lin4), paste0("busco_missing_",lin4))

#Percentage of the genome included in the number of expected chromsome
#Load $id_scaffolds_final.chrom.sizes
if (file.exists(args[40])) {
  chrom_size_table = read.table(args[40], sep = " ")
  colnames(chrom_size_table) = c("chromsome", "size")
  ass_length_adding_scaffold_bp = sum(as.numeric(chrom_size_table$size))
  ass_length_adding_scaffold_gb = ass_length_adding_scaffold_bp / 1000000000
  #Extract table with the expected number of chromosome
  #chrom_num = 19
  subset_chrom_size_table = chrom_size_table [1:chrom_num,]
  #cumulate size
  cum_length_n_chrom_bp = sum(subset_chrom_size_table$size)
  #hap_gen_size_Gb = 1.38
  perc_ass_assigned_to_expected_n_of_chr = cum_length_n_chrom_bp*100/ass_length_adding_scaffold_bp
} else {
  ass_length_adding_scaffold_gb = "NA"
  perc_ass_assigned_to_expected_n_of_chr = "NA"

}


##Write file

overview_table = cbind(xxid, taxon_name, taxon_taxid, lineage, lineage2, lineage3, lineage4, hap_gen_size_Gb, ploidy, chrom_num, pacbio_n_files, pacbio_input_type, pacbio_rq, hic_n_paired_files, pe150_n_paired_files, ont_n_files, pipeline_version, assembly_method, assembly_secondary_mode, polishing_method, purging_method, scaffolding_method, manual_curation, mitohifi, lonqc_overview, kraken_pacbio_overview, kraken_hic_overview, quast_contig_overview, quast_contig_purged_overview, quast_scaffold_overview, busco_lin1_overview, busco_lin2_overview, busco_lin3_overview, busco_lin4_overview, ass_length_adding_scaffold_gb, perc_ass_assigned_to_expected_n_of_chr)

write.table(overview_table, file="overview_results.tsv", quote=FALSE, row.names = FALSE, sep="\t")



