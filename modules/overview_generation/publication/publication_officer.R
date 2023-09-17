#officer_url = "https://cran.r-project.org/src/contrib/officer_0.6.2.tar.gz"
#install.packages(officer_url, repos=NULL, type="source", lib="/projects/cbp/scratch/singularity/")

#library("officer")
library("officer", lib.loc="/projects/cbp/scratch/singularity/")
#library(magrittr)

args <- commandArgs(trailingOnly = TRUE)

outdir = "/projects/cbp/scratch/Arctic_surfclam_005/V1/"

#Extract as much info from the overview table
sample_overview_file = read.table(args[1], sep = "\t", header=T)
#xxid, taxon_name, taxon_taxid, lineage, lineage2, lineage3, lineage4, hap_gen_size_Gb, ploidy, chrom_num, pacbio_n_files, pacbio_input_type, pacbio_rq, hic_n_paired_files, pe150_n_paired_files, ont_n_files, pipeline_version, assembly_method, assembly_secondary_mode, polishing_method, purging_method, scaffolding_method, manual_curation, mitohifi, lonqc_overview, kraken_pacbio_overview, kraken_hic_overview, quast_contig_overview, quast_contig_purged_overview, quast_scaffold_overview, busco_lin1_overview, busco_lin2_overview, busco_lin3_overview, busco_lin4_overview, ass_length_adding_scaffold_gb, perc_ass_assigned_to_expected_n_of_chr

#lonqc_overview = cbind(pacbio_n_reads_longqc, pacbio_longest_read_longqc_kb, pacbio_mean_read_length_longqc_kb, pacbio_coverage_x)
#kraken_pacbio_overview = cbind(kraken_pacbio_actinopteri_perc, kraken_pacbio_amphibia_perc, kraken_pacbio_aves_perc, kraken_pacbio_bivalvia_perc, kraken_pacbio_insecta_perc, kraken_pacbio_magnoliopsida_perc, kraken_pacbio_mammallia_perc, kraken_pacbio_unclassified_perc, kraken_pacbio_other_perc)
#kraken_hic_overview = cbind(kraken_hic_actinopteri_perc, kraken_hic_amphibia_perc, kraken_hic_aves_perc, kraken_hic_bivalvia_perc, kraken_hic_insecta_perc, kraken_hic_magnoliopsida_perc, kraken_hic_mammallia_perc, kraken_hic_unclassified_perc, kraken_hic_other_perc, kraken_hic_root_numer_reads, kraken_hic_unclassified_numer_reads, kraken_hic_root_unclassified_numer_reads)
#quast_contig_overview = cbind (hap1_contig_n50_quast_mb, hap1_contig_l50_quast, hap1_contig_n90_quast_mb, hap1_contig_l90_quast, hap1_contig_assembly_length_quast_gb, hap1_contig_largest_contig_quast_mb, hap1_contig_number_quast, hap2_contig_n50_quast_mb, hap2_contig_l50_quast, hap2_contig_n90_quast_mb, hap2_contig_l90_quast, hap2_contig_assembly_length_quast_gb, hap2_contig_largest_contig_quast_mb, hap2_contig_number_quast)
#quast_contig_purged_overview = cbind (hap1_contig_purged_n50_quast_mb, hap1_contig_purged_l50_quast, hap1_contig_purged_n90_quast_mb, hap1_contig_purged_l90_quast, hap1_contig_purged_assembly_length_quast_gb, hap1_contig_purged_largest_contig_quast_mb, hap1_contig_purged_number_quast, hap2_contig_purged_n50_quast_mb, hap2_contig_purged_l50_quast, hap2_contig_purged_n90_quast_mb, hap2_contig_purged_l90_quast, hap2_contig_purged_assembly_length_quast_gb, hap2_contig_purged_largest_contig_quast_mb, hap2_contig_purged_number_quast)
#quast_scaffold_overview = cbind (hap1_scaffold_n50_quast_mb, hap1_scaffold_l50_quast, hap1_scaffold_n90_quast_mb, hap1_scaffold_l90_quast, hap1_scaffold_assembly_length_quast_gb, hap1_scaffold_largest_contig_quast_mb, hap1_scaffold_number_quast, hap2_scaffold_n50_quast_mb, hap2_scaffold_l50_quast, hap2_scaffold_n90_quast_mb, hap2_scaffold_l90_quast, hap2_scaffold_assembly_length_quast_gb, hap2_scaffold_largest_contig_quast_mb, hap2_scaffold_number_quast)
#busco_lin1_overview = cbind(lin1, busco_complete_single_lin1, busco_complete_duplicated_lin1, busco_complete_single_duplicated_lin1, busco_fragmented_lin1, busco_missing_lin1)


#Common name is xxid without with _ replace by spaces and removed the last argument
id = sample_overview_file$xxid
common_name_temp = sub("_", " ", sample_overview_file$xxid) 
common_name = substr(common_name_temp,1,nchar(common_name_temp)-4)

scientific_name = sub("_", " ", sample_overview_file$taxon_name)

genome_size = sample_overview_file$hap_gen_size_Gb
chrom_number = sample_overview_file$chrom_num

pacbio_minrq = sample_overview_file$pacbio_rq
pipeline_version = sample_overview_file$pipeline_version

assembly_method = sample_overview_file$assembly_method
assembly_secondary_mode = sample_overview_file$assembly_secondary_mode
polishing_method = sample_overview_file$polishing_method
purging_method = sample_overview_file$purging_method
scaffolding_method = sample_overview_file$scaffolding_method
manual_curation = sample_overview_file$manual_curation
mitohifi = sample_overview_file$mitohifi

pacbio_coverage = sample_overview_file$pacbio_coverage_x

contig_N50 = sample_overview_file$hap1_contig_purged_n50_quast_mb
contig_number = sample_overview_file$hap1_contig_purged_number_quast

scaffold_N50 = sample_overview_file$hap1_scaffold_n50_quast_mb
scaffold_N90 = sample_overview_file$hap1_scaffold_n90_quast_mb
scaffold_number = sample_overview_file$hap1_scaffold_number_quast
assembly_length = sample_overview_file$hap1_scaffold_assembly_length_quast_gb
length_longest_scaffold = sample_overview_file$hap1_scaffold_largest_contig_quast_mb

lin1 = sample_overview_file$lin1
Busco_lin1 = sample_overview_file$busco_complete_single_duplicated_lin1
Busco_complete_lin1 = sample_overview_file$busco_complete_single_lin1
Busco_frag_lin1 = sample_overview_file$busco_fragmented_lin1
Busco_dup_lin1 = sample_overview_file$busco_complete_duplicated_lin1
Busco_missing_lin1 = sample_overview_file$busco_missing_lin1

perc_ass_chr = sample_overview_file$perc_ass_assigned_to_expected_n_of_chr

assembly_name = "<To complete: assembly name>"
assembly_number = "<To complete: assembly number>"


doc_1 <-  read_docx()
doc_1 <- body_add_par(doc_1, "First page", style = "heading 1")

# Add a title to the document
doc_1 <- body_add_par(doc_1, "Title", style = "heading 2")
doc_1 <- body_add_par(doc_1, paste0("The genome sequence of the ",common_name, ", ", scientific_name), style = "Normal")

# Add a sub title (author list)
doc_1 <- body_add_par(doc_1, "Author list and affiliation", style = "heading 2")
doc_1 <- body_add_par(doc_1, "<To complete>", style = "Normal")

# Add a sub title (abstract)
doc_1 <- body_add_par(doc_1, "Abstract", style = "heading 2")
doc_1 <- body_add_par(doc_1, paste0("We present a genome assembly of ", scientific_name, " (the ", common_name, "; <TO COMPLETE : lineages>. The genome sequence is ", genome_size, " gigabases in size. The assembly is composed of ", scaffold_number, ", with a N50 of ", scaffold_N50, " and a BUSCO score of ", Busco_lin1, " for lineage ", lin1, "."), style = "Normal")

# Add a sub title (Keywords)
doc_1 <- body_add_par(doc_1, "Keywords", style = "heading 2")
doc_1 <- body_add_par(doc_1, paste0(scientific_name, ", ", common_name, ", genome sequence"), style = "Normal")

# Add a sub title (Taxonomy)
doc_1 <- body_add_par(doc_1, "Species taxonomy", style = "heading 2")
doc_1 <- body_add_par(doc_1, "<To complete>", style = "Normal")

#Article
doc_1 <- body_add_par(doc_1, "Article", style = "heading 1")

# Add a sub title (Intro)
doc_1 <- body_add_par(doc_1, "Introduction", style = "heading 2")
doc_1 <- body_add_par(doc_1, paste0("The ", common_name, ", ", scientific_name, "<TO COMPLETE: Description of the specie, geographical distribution, relevance to biosystems, endenegered status (if applicable), etc>. The genome of ", scientific_name, " was sequenced as part of the Canadian BioGenome Project (CBP). The ", scientific_name, " genome will provide insights into genomic diversity and architecture, and inform conservation genomics applications."), style = "Normal") 

# Add a sub title (Method)
doc_1 <- body_add_par(doc_1, "Methods", style = "heading 2")
doc_1 <- body_add_par(doc_1, "Sample collection", style = "heading 3")

doc_1 <- body_add_par(doc_1, "<TO COMPLETE : Include Number of individuals, geographical location, sex, development stage, tissue, conservation method, ethic and permit numbers (if applicable), etc.>", style = "Normal")


doc_1 <- body_add_par(doc_1, "Sample extraction, library construction and sequencing", style = "heading 3")
doc_1 <- body_add_par(doc_1, "<TO COMPLETE (Possible text) : High-molecular weight (HMW) DNA was extracted from <tissue type> using the <extraction kit information> at <extraction center>. PacBio genome libraries were constructed using <library kit info> and sequenced on <PacBio instrument> at <location PacBio sequencing>. A Hi-C library was constructed using the <Hi-C kit info> at <location> and subjected to PE150 sequencing on an <Hi-C sequencer> instrument at <location>. If short read, RNA or other data is generated for this specie, add the information>", style = "Normal")

doc_1 <- body_add_par(doc_1, "Genome assembly", style = "heading 3")


if ((assembly_method == "hifiasm") && (purging_method == "purge_dups") && (scaffolding_method == "yahs") && (manual_curation == "none")) {
	# Hifiasm + purge_dups + yahs (no manual curation)
	doc_1 <- body_add_par(doc_1, paste0("Assembly was carried out using hifiasm with the ", assembly_secondary_mode, " mode (Cheng et al., 2021). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."), style = "Normal")
} else if ((assembly_method == "hifiasm") && (purging_method == "purge_dups") && (scaffolding_method == "yahs") && (manual_curation == "yes")) {
        # Hifiasm + purge_dups + yahs + manual curation using JUICER
	doc_1 <- body_add_par(doc_1, paste0("Assembly was carried out using hifiasm with the ", assembly_secondary_mode, " (Cheng et al., 2021). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). Manual curation was performed using Juicer (Durand et al., 2016). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."), style = "Normal")
} else if ((assembly_method == "flye") && (purging_method == "purge_dups") && (scaffolding_method == "yahs") && (manual_curation == "none")) {
	# Flye + purge_dups + yahs  (no manual curation)
        doc_1 <- body_add_par(doc_1, paste0("Assembly was carried out using flye (Kolmogorov et al., 2019). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."), style = "Normal")
} else if ((assembly_method == "flye") && (purging_method == "purge_dups") && (scaffolding_method == "yahs") && (manual_curation == "yes")) {
        # Flye + purge_dups + yahs + manual curation using JUICER
	doc_1 <- body_add_par(doc_1, paste0("Assembly was carried out using flye (Kolmogorov et al., 2019). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). Manual curation was performed using Juicer (Durand et al., 2016). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."), style = "Normal")
} else {
	doc_1 <- body_add_par(doc_1, paste0("<TO COMPLETE, the method appears different from the mainstream pipeline>"), style = "Normal")
}

if (mitohifi == "yes") {
	doc_1 <- body_add_par(doc_1, paste0("The mitochondrial genome was assembled using MitoHiFi (Uliano-Silva et al., 2021)."), style = "Normal")
}

##MAY NEED TO INCLUDE MERQURY IF FIGURE OR QV USED
#To evaluate the assembly, MerquryFK was used to estimate consensus quality (QV) scores and k-mer completeness (Rhie et al., 2020).

#May need to include barcoding

# Add a sub title (Results)
doc_1 <- body_add_par(doc_1, "Results", style = "heading 2")
doc_1 <- body_add_par(doc_1, "Genome sequence report", style = "heading 3")

if (polishing_method == "none") {
	doc_1 <- body_add_par(doc_1, paste0("The genome of <TO COMPLETE : number of individual(s)> ", common_name, " collected from <TO COMPLETE : location> was sequenced. A total of ", pacbio_coverage, "-fold coverage in PacBio long reads (read quality > ", pacbio_minrq, ") were generated. Contigs were then scaffolded with Hi-C data <TO COMPLETE : < generated from the same individual> or <generated from a different individual collected from <location>>. The final assembly has a total length of ", assembly_length, "Gb organized in ", scaffold_number, "sequence scaffolds with a scaffold N50 of ", scaffold_N50, " Mb (Table 1). ", perc_ass_chr, " of the assembly sequence was assigned to the ", chrom_number, "longest scaffolds representing the species’ known ", chrom_number, " autosomes (<TO COMPLETE : reference>) (numbered by sequence length; Figure 1–Figure 4; Table 2). Determining gene coverage using BUSCO, we estimated ", Busco_lin1, "% gene completeness using the ", lin1, " reference set (Manni et al., 2021)."), style = "Normal")
} else {
	doc_1 <- body_add_par(doc_1, paste0("<TO COMPLETE, the method appears different from the mainstream pipeline>"), style = "Normal")
}

doc_1 <- body_add_par(doc_1, "Genome annotation", style = "heading 3")
doc_1 <- body_add_par(doc_1, paste0("Annotation for the ", common_name, " genome assembly (", assembly_name, " (", assembly_number, ")) was generated by the Ensembl Rapid Release gene annotation pipeline (Aken et al., 2016). The resulting Ensembl annotation includes <TO COMPLETE> transcripts assigned to <TO COMPLETE> coding and <TO COMPLETE> non-coding genes (", scientific_name, " - Ensembl Rapid Release). The ", common_name, " assembly was also annotated for <TO COMPLETE> protein sequences using RefSeq (<TO COMPLETE>)."), style = "Normal")


# Add a sub title (Data availability)
doc_1 <- body_add_par(doc_1, "Data availability", style = "heading 2")
doc_1 <- body_add_par(doc_1, "Underlying data", style = "heading 3")
doc_1 <- body_add_par(doc_1, paste0("National Centre for Biotechnology Information BioProject: ", common_name, " (", scientific_name, ") genome sequencing and assembly, ", assembly_name, ". Accession number: <TO COMPLETE>."), style = "Normal")
doc_1 <- body_add_par(doc_1, paste0("The genome sequence is released openly for reuse. The ", scientific_name, " genome sequencing initiative is part of the Canadian BioGenome Project. All raw sequence data and the assembly have been deposited in INSDC databases. The genome is annotated through the Reference Sequence (RefSeq) database in BioProject accession number <TO COMPLETE BIOPROJECT NUMBER>. Raw data and assembly accession identifiers are reported in Table 1."), style = "Normal") 


# Add a sub title (References)
doc_1 <- body_add_par(doc_1, "References", style = "heading 2")
doc_1 <- body_add_par(doc_1, "<TO DO>", style = "Normal")

#Figures
doc_1 <- body_add_par(doc_1, "Figures", style = "heading 1")

#FIG 1 SNAIL PLOT
doc_1 <- body_add_par(doc_1, paste0("Figure 1. Genome assembly of ", scientific_name, ", ", assembly_name, ": metrics."), style = "heading 2")
#fig1 <- file.path("/projects/cbp/scratch/Arctic_surfclam_005/V1/QC/overview/Fig1.png")
if (file.exists(file.path(paste0(outdir,"/blobtools/",id,".snail.png")))){
	fig1 <- file.path(paste0(outdir,"/blobtools/",id,".snail.png"))
	doc_1 <- body_add_img(doc_1, src = fig1, height = 4, width = 4, style = "centered")
} else {
        doc_1 <- body_add_par(doc_1, " ")
	doc_1 <- body_add_par(doc_1, "SNAIL PLOT NOT GENERATED")
        doc_1 <- body_add_par(doc_1, " ")
}
doc_1 <- body_add_par(doc_1, paste0("Snail plot showing N50 metrics, base pair composition and BUSCO gene completeness for ", scientific_name, " (", assembly_name, ") generated from Blobtoolkit v.2.6.4 (Challis et al., 2020). The plot is divided into 1,000 size-ordered bins around the circumference with each bin representing 0.1% of the ", assembly_length, " bp assembly. The distribution of chromosome lengths is shown in dark grey with the plot radius scaled to the longest chromosome present in the assembly (", length_longest_scaffold, " bp) shown in red. Orange and pale-orange arcs show the N50 and N90 chromosome lengths (", scaffold_N50, " bp and ", scaffold_N90, " bp, respectively). The pale grey spiral shows the cumulative chromosome count on a log scale with white scale lines showing successive orders of magnitude. The blue and pale-blue area around the outside of the plot displays the distribution of GC (blue), AT (pale blue) and N (white) percentages using the same bins as the inner plot. A summary of complete (", Busco_complete_lin1, "%), fragmented (", Busco_frag_lin1, "%), duplicated (", Busco_dup_lin1, "%), and missing (", Busco_missing_lin1, "%) BUSCO genes in the ", lin1, " set is show in the top right."), style = "Normal")


#Fig 2 : GC content
doc_1 <- body_add_par(doc_1, paste0("Figure 2. Genome assembly of ", scientific_name, ", ", assembly_name, ": GC-content."), style = "heading 2")
doc_1 <- body_add_par(doc_1, " ")

#fig2 <- file.path("/projects/cbp/scratch/Arctic_surfclam_005/V1/QC/overview/Fig1.png")
if (file.exists(file.path(paste0(outdir,"/blobtools/",id,".blob.circle.png")))){
        fig2 <- file.path(paste0(outdir,"/blobtools/",id,".blob.circle.png"))
        doc_1 <- body_add_img(doc_1, src = fig2, height = 4, width = 4, style = "centered")
} else {
        doc_1 <- body_add_par(doc_1, "GC COVERAGE PLOT NOT GENERATED")
}
doc_1 <- body_add_par(doc_1, " ")
doc_1 <- body_add_par(doc_1, paste0("GC-coverage plot of ", scientific_name, " (", assembly_name, ") generated from Blobtoolkit v.2.6.4 (Challis et al., 2020). Scaffolds are coloured by phylum with <TO COMPLETE> represented by blue and no-hit represented by pale blue. Circles are sized in proportion to scaffold length. Histograms show the distribution of scaffold length sum along each axis."), style = "Normal")


#Fig 3 cumulative plot
doc_1 <- body_add_par(doc_1, paste0("Figure 3. Genome assembly of ", scientific_name, ", ", assembly_name, ": cumulative sequence length."), style = "heading 2")
doc_1 <- body_add_par(doc_1, " ")
#fig3 <- file.path("/projects/cbp/scratch/Arctic_surfclam_005/V1/QC/overview/Fig1.png")
if (file.exists(file.path(paste0(outdir,"/blobtools/",id,".cumulative.png")))){
        fig3 <- file.path(paste0(outdir,"/blobtools/",id,".cumulative.png"))
        doc_1 <- body_add_img(doc_1, src = fig3, height = 4, width = 4, style = "centered")
} else {
        doc_1 <- body_add_par(doc_1, "CUMULATIVE PLOT NOT GENERATED")
}
doc_1 <- body_add_par(doc_1, " ")
doc_1 <- body_add_par(doc_1, paste0("Cumulative sequence length of ", scientific_name, " (", assembly_name, ") generated from Blobtoolkit v.2.6.4 (Challis et al., 2020). The grey line shows the cumulative length for all scaffolds. Coloured lines show cumulative lengths of scaffolds assigned to each phylum using the BUSCO genes tax rule, with <TO COMPLETE> represented by blue and no-hit represented by pale blue."), style = "Normal")


## Fig 4 Hi-C Map
doc_1 <- body_add_par(doc_1, paste0("Figure 4. Genome assembly of ", scientific_name, ", ", assembly_name, ": Hi-C contact map."), style = "heading 2")
doc_1 <- body_add_par(doc_1, " ")
if (file.exists(file.path(paste0(outdir,"/QC/pretext/",id,"FullMap.png")))){
        fig4 <- file.path(paste0(outdir,"/QC/pretext/",id,"FullMap.png"))
        doc_1 <- body_add_img(doc_1, src = fig4, height = 4, width = 4, style = "centered")
} else { 
        doc_1 <- body_add_par(doc_1, "PRETEXT HIC MAP NOT GENERATED")
}
doc_1 <- body_add_par(doc_1, " ")
doc_1 <- body_add_par(doc_1, paste0("HiC contact map of ", assembly_name, " assembly visualized using PRETEXT <REFERENCE>. Scaffolds are shown in order of size from left to right and top to bottom."), style = "Normal")

#Tables
doc_1 <- body_add_par(doc_1, "Tables", style = "heading 1")
doc_1 <- body_add_par(doc_1, paste0("Table 1. Genome data for ", scientific_name, ", ", assembly_name, "."), style = "heading 2")
table_1_1_first_col = c("Project accession data", "Assembly identifier", "Species", "Specimen", "NCBI Taxonomy ID", "BioProject", "BioSample ID", "Isolate Information")
table1_1_second_col = c(" ",  assembly_name, scientific_name, "<TO COMPLETE>", "<TO COMPLETE>", "<TO COMPLETE>", "<TO COMPLETE>", "<TO COMPLETE>")
data_table1_1 = data.frame(table_1_1_first_col, table1_1_second_col)
table_1_2_first_col = c("Raw data accessions", "PacBio data", "Hi-C Illumina")
table1_2_second_col = c(" ", "<TO COMPLETE>", "<TO COMPLETE>")
data_table1_2 = data.frame(table_1_2_first_col, table1_2_second_col)
table_1_3_first_col = c("Genome assembly", "Assembly accession", "Assembly name", "Span (Mb)", "Number of contigs", "Contig N50 length (Mb)", "Number of scaffolds", "Scaffold N50 length (Mb)", "Longest scaffold (Mb)", "BUSCO* genome score")
table1_3_second_col = c(" ", "<To complete>", assembly_name, assembly_length, contig_number, contig_N50, scaffold_number, scaffold_N50, length_longest_scaffold, paste0("C:",Busco_lin1,"%[S:",Busco_complete_lin1,"%,D:",Busco_dup_lin1,"%],F:",Busco_frag_lin1,"%,M:",Busco_missing_lin1,"%,n=<To complete>"))
data_table1_3 = data.frame(table_1_3_first_col, table1_3_second_col)

doc_1 <- body_add_par(doc_1, " ")
doc_1 <- body_add_table(doc_1, value = data_table1_1, style = "table_template", header = FALSE, first_column=TRUE, first_row = TRUE, alignment="l")
doc_1 <- body_add_par(doc_1, paste0(" "), style = "Normal")
doc_1 <- body_add_table(doc_1, value = data_table1_2, style = "table_template", header = FALSE, first_column=TRUE, first_row = TRUE, alignment="l")
doc_1 <- body_add_par(doc_1, paste0(" "), style = "Normal")
doc_1 <- body_add_table(doc_1, value = data_table1_3, style = "table_template", header = FALSE, first_column=TRUE, first_row = TRUE, alignment="l")
doc_1 <- body_add_par(doc_1, paste0("* BUSCO scores based on the ", lin1, " BUSCO set using v5.0.0. C= complete [S= single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison."), style = "Normal")

doc_1 <- body_add_par(doc_1, paste0("Table 2. Chromosomal pseudomolecules in the genome assembly of ", scientific_name, ", ", assembly_name, "."), style = "heading 2")
doc_1 <- body_add_par(doc_1, "TO DO")

doc_1 <- body_add_par(doc_1, "3. Software tools used.", style = "heading 2")
doc_1 <- body_add_par(doc_1, "TO DO")



print(doc_1, target=paste0(common_name, "_publication_template.docx"))











# Add a hyperlink
#list of hyperlinks : Bioproject, Ensembl release, CBP website? Link to pipeine? References? Tables? Figures?
#my_link <- pot('Click here to visit STHDA web site!', 
#          hyperlink = 'http://www.sthda.com/english',
#          format=textBoldItalic(color = 'blue', underline = TRUE ))
#addParagraph(doc, my_link)






#http://www.sthda.com/english/wiki/create-and-format-word-documents-using-r-software-and-reporters-package

#https://ardata-fr.github.io/officeverse/officer-for-word.html
