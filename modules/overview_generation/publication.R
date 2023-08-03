install.packages('ReporteRs')
library('ReporteRs')

# Create a word document to contain R outputs
doc <- docx()
# Add a title to the document
doc <- addTitle(doc, paste0("The genome sequence of the ",common_name, ", ", scientific_name) level=1)

# Add a sub title (author list)
doc <- addTitle(doc, "Author list and affiliation", level = 2)
doc <- addParagraph(doc, "<TO COMPLETE>")

# Add a sub title (abstract)
doc <- addTitle(doc, "Abstract", level = 2)
doc <- addParagraph(doc, paste0("We present a genome assembly of ", scientific_name, " (the ", common_name, "; <TO COMPLETE : lineages>. The genome sequence is ", genome_size, " gigabases in size. The assembly is composed of ", scaffold_number, ", with a N50 of ", scaffold_N50, " and a BUSCO score of ", Busco_lin1, " for lineage ", lin1, "."))

# Add a sub title (Keywords)
doc <- addTitle(doc, "Keywords", level = 2)
doc <- addParagraph(doc, paste0(scientific_name, ", ", common_name, ", genome sequence"))

# Add a sub title (Taxonomy)
doc <- addTitle(doc, "Species taxonomy", level = 2)
doc <- addParagraph(doc, "<TO DO>")

# Add a sub title (Intro)
doc <- addTitle(doc, "Introduction", level = 2)
doc <- addParagraph(doc, paste0("The ", common_name, ", ", scientific_name, "<TO COMPLETE: Description of the specie, geographical distribution, relevance to biosystems, endenegered status (if applicable), etc>. The genome of ", scientific_name, " was sequenced as part of the Canadian BioGenome Project (CBP). The ", scientific_name, " genome will provide insights into genomic diversity and architecture, and inform conservation genomics applications."))

# Add a sub title (Method)
doc <- addTitle(doc, "Methods", level = 2)
doc <- addTitle(doc, "Sample collection", level = 3)
doc <- addParagraph(doc, "<TO COMPLETE : Include Number of individuals, geographical location, sex, development stage, tissue, conservation method, ethic and permit numbers (if applicable), etc.>")

doc <- addTitle(doc, "Sample extraction, library construction and sequencing", level = 3)
doc <- addParagraph(doc, "<TO COMPLETE (Possible text) : High-molecular weight (HMW) DNA was extracted from <tissue type> using the <extraction kit information> at <extraction center>. PacBio genome libraries were constructed using <library kit info> and sequenced on <PacBio instrument> at <location PacBio sequencing>. A Hi-C library was constructed using the <Hi-C kit info> at <location> and subjected to PE150 sequencing on an <Hi-C sequencer> instrument at <location>. If short read, RNA or other data is generated for this specie, add the information>")

doc <- addTitle(doc, "Genome assembly", level = 3)
if ((assembly_method == "hifiasm") & (purging_method == "purge_dups") & (scaffolding_method == "yahs")  & (manual_curation == "none")) {
	# Hifiasm + purge_dups + yahs (no manual curation)
	doc <- addParagraph(doc, paste0("Assembly was carried out using hifiasm with the ", assembly_secondary_mode, " mode (Cheng et al., 2021). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."))
} else if ((assembly_method == "hifiasm") & (purging_method == "purge_dups") & (scaffolding_method == "yahs") & (manual_curation == "yes")) {
        # Hifiasm + purge_dups + yahs + manual curation using JUICER
        doc <- addParagraph(doc, paste0("Assembly was carried out using hifiasm with the ", assembly_secondary_mode, " (Cheng et al., 2021). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). Manual curation was performed using Juicer (Durand et al., 2016). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."))
} else if ((assembly_method == "flye") & (purging_method == "purge_dups") & (scaffolding_method == "yahs") & (manual_curation == "none")) {
        # Flye + purge_dups + yahs  (no manual curation)
        doc <- addParagraph(doc, paste0("Assembly was carried out using flye (Kolmogorov et al., 2019). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."))
} else if ((assembly_method == "flye") & (purging_method == "purge_dups") & (scaffolding_method == "yahs") & (manual_curation == "yes")) {
        # Flye + purge_dups + yahs + manual curation using JUICER
        doc <- addParagraph(doc, paste0("Assembly was carried out using flye (Kolmogorov et al., 2019). Purging was done using purge_dups (Guan et al., 2020). Scaffolding with Hi-C data was carried out using YaHS (Zhou et al., 2023). The Hi-C contact maps were generated using Pretext (Harry, 2022). Manual curation was performed using Juicer (Durand et al., 2016). The final sequence was analyzed using BlobToolKit (Challis et al., 2020) and BUSCO scores were generated (Manni et al., 2021; Simão et al., 2015). The steps listed before as well as the generation of the manuscript template are organized in a nextflow pipeline available at : https://github.com/bcgsc/Canadian_Biogenome_Project. Software tools and versions are listed in Table 3."))
} else {
        doc <- addParagraph(doc, paste0("<TO COMPLETE, the method appears different from the mainstream pipeline>"))
}

if (mitohifi == "yes") {
        doc <- addParagraph(doc, paste0("The mitochondrial genome was assembled using MitoHiFi (Uliano-Silva et al., 2021)."))
}

##MAY NEED TO INCLUDE MERQURY IF FIGURE OR QV USED
#To evaluate the assembly, MerquryFK was used to estimate consensus quality (QV) scores and k-mer completeness (Rhie et al., 2020).

#May need to include barcoding

# Add a sub title (Results)
doc <- addTitle(doc, "Results", level = 2)
doc <- addTitle(doc, "Genome sequence report", level = 3)
if (polishing_method == "none") {
	doc <- addParagraph(doc, paste0("The genome of <TO COMPLETE : number of individual(s)> ", common_name, " collected from <TO COMPLETE : location> was sequenced. A total of ", pacbio_coverage, "-fold coverage in PacBio long reads (read quality > ", pacbio_minrq, ") were generated. Contigs were then scaffolded with Hi-C data <TO COMPLETE : < generated from the same individual> or <generated from a different individual collected from <location>>. The final assembly has a total length of ", assembly_length, "Gb organized in ", scaffold_number, "sequence scaffolds with a scaffold N50 of ", scaffold_N50, " Mb (Table 1). ", perc_ass_chr, " of the assembly sequence was assigned to the ", chrom_number, "longest scaffolds representing the species’ known ", chrom_number, " autosomes (<TO COMPLETE : reference>) (numbered by sequence length; Figure 1–Figure 4; Table 2). Determining gene coverage using BUSCO, we estimated ", Busco_lin1, "% gene completeness using the ", lin1, " reference set (Manni et al., 2021)."))
} else {
        doc <- addParagraph(doc, paste0("<TO COMPLETE, the method appears different from the mainstream pipeline>"))
}


doc <- addTitle(doc, "Genome annotation", level = 3)
doc <- addParagraph(doc, paste0("Annotation for the ", common_name, " genome assembly (", assembly_name, " (", assembly_number, ")) was generated by the Ensembl Rapid Release gene annotation pipeline (Aken et al., 2016). The resulting Ensembl annotation includes <TO COMPLETE> transcripts assigned to <TO COMPLETE> coding and <TO COMPLETE> non-coding genes (", scientific name, " - Ensembl Rapid Release). The ", common_name, " assembly was also annotated for <TO COMPLETE> protein sequences using RefSeq (<TO COMPLETE>)."))


# Add a sub title (Data availability)
doc <- addTitle(doc, "Data availability", level = 2)
doc <- addTitle(doc, "Underlying data", level = 3)
doc <- addParagraph(doc, paste0("National Centre for Biotechnology Information BioProject: ", common_name, " (", scientific_name, ") genome sequencing and assembly, ", assembly name, ". Accession number: <TO COMPLETE>.")
doc <- addParagraph(doc, paste0("The genome sequence is released openly for reuse. The ", scientific_name, " genome sequencing initiative is part of the Canadian BioGenome Project. All raw sequence data and the assembly have been deposited in INSDC databases. The genome is annotated through the Reference Sequence (RefSeq) database in BioProject accession number <TO COMPLETE BIOPROJECT NUMBER>. Raw data and assembly accession identifiers are reported in Table 1."))


# Add a sub title (References)
doc <- addTitle(doc, "References", level = 2)
doc <- addParagraph(doc, "<TO DO>")

#Figures
doc <- addTitle(doc, "Figures", level = 1)
doc <- addTitle(doc, paste0("Figure 1. Genome assembly of ", scientific_name, ", ", assembly_name, ": metrics."), level = 2)
doc <- addImage(doc, "Fig1.png")
doc <- addParagraph(doc, paste0("Snail plot showing N50 metrics, base pair composition and BUSCO gene completeness for ", scientific_name, " (", assembly_name, ") generated from Blobtoolkit v.2.6.4 (Challis et al., 2020). The plot is divided into 1,000 size-ordered bins around the circumference with each bin representing 0.1% of the ", assembly_length, " bp assembly. The distribution of chromosome lengths is shown in dark grey with the plot radius scaled to the longest chromosome present in the assembly (", length_longest_scaffold, " bp) shown in red. Orange and pale-orange arcs show the N50 and N90 chromosome lengths (", scaffold_N50, " bp and ", scaffold_N90, " bp, respectively). The pale grey spiral shows the cumulative chromosome count on a log scale with white scale lines showing successive orders of magnitude. The blue and pale-blue area around the outside of the plot displays the distribution of GC (blue), AT (pale blue) and N (white) percentages using the same bins as the inner plot. A summary of complete (", Busco_complete_lin1, "%), fragmented (", Busco_frag_lin1, "%), duplicated (", Busci_dup_lin1, "%), and missing (", Busco_missing_lin1, "%) BUSCO genes in the ", lin1, " set is show in the top right."))

doc <- addTitle(doc, paste0("Figure 2. Genome assembly of ", scientific_name, ", ", assembly_name, ": GC-content."), level = 2)
doc <- addImage(doc, "Fig2.png")
doc <- addParagraph(doc, paste0("GC-coverage plot of ", scientific_name, " (", assembly_name, ") generated from Blobtoolkit v.2.6.4 (Challis et al., 2020). Scaffolds are coloured by phylum with <TO COMPLETE> represented by blue and no-hit represented by pale blue. Circles are sized in proportion to scaffold length. Histograms show the distribution of scaffold length sum along each axis."))

doc <- addTitle(doc, paste0("Figure 3. Genome assembly of ", scientific_name, ", ", assembly_name, ": cumulative sequence length."), level = 2)
doc <- addImage(doc, "Fig3.png")
doc <- addParagraph(doc, paste0("Cumulative sequence length of ", scientific_name, " (", assembly_name, ") generated from Blobtoolkit v.2.6.4 (Challis et al., 2020). The grey line shows the cumulative length for all scaffolds. Coloured lines show cumulative lengths of scaffolds assigned to each phylum using the BUSCO genes tax rule, with <TO COMPLETE> represented by blue and no-hit represented by pale blue."))

doc <- addTitle(doc, paste0("Figure 4. Genome assembly of ", scientific_name, ", ", assembly_name, ": Hi-C contact map."), level = 2)
doc <- addImage(doc, "Fig4.png")
doc <- addParagraph(doc, paste0("HiC contact map of ", assembly_name, " assembly visualized using PRETEXT <REFERENCE>. Scaffolds are shown in order of size from left to right and top to bottom."))

#Tables
doc <- addTitle(doc, "Tables", level = 1)
doc <- addTitle(doc, paste0("Table 1. Genome data for ", scientific_name, ", ", assembly_name, "."), level = 2)
data_table1 <- matrix(c("Project accession data", " ", "Assembly identifier", assembly_name, "Species", scientific_name, "Specimen", "<TO COMPLETE>", "NCBI Taxonomy ID", "<TO COMPLETE>", "BioProject", "<TO COMPLETE>", "BioSample ID", "<TO COMPLETE>", "Isolate Information", "<TO COMPLETE>"), ncol=2, byrow=TRUE)
Table1 <- vanilla.table(data_table1)
Table1 <- setZebraStyle(Table1, odd = '#eeeeee', even = 'white')
doc <- addFlexTable( doc, Table1)
doc <- addParagraph(doc, paste0("* BUSCO scores based on the ", lin1, " BUSCO set using v5.0.0. C= complete [S= single copy, D=duplicated], F=fragmented, M=missing, n=number of orthologues in comparison."))











# Add a hyperlink
#list of hyperlinks : Bioproject, Ensembl release, CBP website? Link to pipeine? References? Tables? Figures?
#my_link <- pot('Click here to visit STHDA web site!', 
#          hyperlink = 'http://www.sthda.com/english',
#          format=textBoldItalic(color = 'blue', underline = TRUE ))
#doc <- addParagraph(doc, my_link)






#http://www.sthda.com/english/wiki/create-and-format-word-documents-using-r-software-and-reporters-package

# Write the Word document to a file 
writeDoc(doc, file = paste0(common_name, "_publication_template.docx"))
