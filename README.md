# CBP Genome Assembly pipeline
This pipeline is used by the [Canadian Biogenome project](http://earthbiogenome.ca) to generate genome assemblies from a variety of species.

The pipeline is built using [nextflow](https://www.nextflow.io/).

In short, each step of the pipeline is included in a module. Most of the modules uses one container which makes it much easier to maintain and update software dependencies. **Some modules rely on locally installed tools**. Future updates of the pipeline may include better portability.

A lot of the modules available in this pipeline were developed by members of the nf-core/genomeassembler group, if you want to participate, feel free to join the community.

## **Table of Contents**
* **[Input data](#input-data)**
* **[Output data](#output-files)**
* **[Process](#process)**
  * [Running the pipeline with test data](#running-the-pipeline-with-test-data-(will-work-once-the-repo-is-public))
  * [Running the pipeline with your own data](#running-the-pipeline-with-your-own-data)
* **[Credits](#credits)**
* **[Details on the test dataset](#details-on-the-test-dataset)**



## **Input data**
The pipeline was developped to take as input PacBio ccs files (bam) and Hi-C files (fastq.gz). The pipeline also support the inclusion of nanopore data and short-reads for polishing.

The pipeline also require the specie [Taxonomy ID](https://en.wikipedia.org/wiki/Help:Taxon_identifiers#:~:text=Taxon%20identifiers%20are%20identifiers%20assigned,reference%20for%20each%20catalogued%20taxon.) in order to query [GoaT database](https://goat.genomehubs.org) and retrieve additional information on the specie such as the scientific name, the ploidy or the genome size.

You can find the taxonomic id of your specie of interest on NCBI.
For example, the taxonoic ID for the Steller Sea Lion is 34886 [NCBI link](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=34886)



## **Output files**
The pipeline generates many files and intermediate files, most are self explanatory. 



## **Process**
An overview of the pipeline is visible on the following subway map. Some parts of the pipeline may have been commented out in this version as they relied on locaaly installed software. The code is still available in case you also want to locally install the software and try it out.

By default, the pipeline will use hifiasm with PacBio data for the assembly, and if Hi-C data is available, YAHS is used for the scaffolding.
Other assembler and scaffolder are available within the pipeline, to change, you need to edit the nextflow.config file.

Software used that would require local installation:

- [LongQC](https://github.com/yfukasawa/LongQC)
- [MitoHifi](https://github.com/marcelauliano/MitoHiFi)
- [Juicer](https://github.com/aidenlab/juicer)


Software that relies on locally downloaded files / databases :

- [Kraken](http://ccb.jhu.edu/software/kraken/)

<p align="center">
    <img title="The Canadian Biogenome Project Workflow" src="res/CBP_workflow.png" width=50%>
</p>
<p align="center">
Figure : Overview of the Canadian Biogenome project assembly pipeline
</p>


## Running the pipeline with test data
To run this pipeline, you need nextflow and conda installed on your system.

A set of test data are available in this repo to allow you to test the pipeline with just one command line:

```
nextflow run bcgsc/Canadian_Biogenome_Project -r V2 -latest -profile conda
```

The outputs are organized in several subfolder that are self-explenatory.



## Running the pipeline with your own data
Clone the repository in your local environment:

```
git clone https://github.com/bcgsc/Canadian_Biogenome_Project.git
cd Canadian_Biogenome_Project
```

Modify the nextflow.config file:
IMPORTANT : In the nextflow config file, to comment out a line, the type is : // (instead of # in bash scripts)

- Indicate the specie ID (ex : "Steller_sea_lion_001) and the taxon ID (ex : 34886)
```
//Specie parameters
        id                      = "Steller_sea_lion_001"
        taxon_taxid             = "34886"
//	related_genome		= "GCA_009762305.2"
	string_telomere		= "TTAGGG"
        pipeline_version        = "V2"
```


- Indicate the location and filenames of the input files (PacBio data, Hi-C data, etc)
For PacBio data, unlaigned bam files are expected.
Reads from the latest PacBio technology (CCS) are expected, if other reads are used, indicate it and use another aligner to generate the assembly (hifiasm, the default assembler in this pipeline, is designed for HiFi reads, PacBio latest technology).
If you have several PacBio SMRT cells for one specie, you can indicate them all (up to 4). No need to merge the data prior to launching the pipeline.
```
//PacBio input
	pacbio_input_type	= "ccs" //'ccs' or 'clr' or 'subreads'
	bam_cell1		= "$baseDir/example_input/subset_covid_hifi.bam"
//      bam_cell2               = "${raw_data_path}/pacbio/"

//HiC Illumina input
        hic_read1               = "$baseDir/example_input/test_1.fastq.gz"
        hic_read2               = "$baseDir/example_input/test_2.fastq.gz"
        Illumina_prefix         = "test"
```

(optionnaly) - As the pipeline relies on hifiasm by default, only reads with rq>0.99 (≥Q20 (HiFi reads only, Probability of incorrect base call : 1 in 100), equivalent of using extracthifi software) are used.
If you want to filter mor or less read, you can change the read quality threashold in the nextf.config file:
```
        pacbio_rq               = "0.99"
```

(optionnaly) - Indicate the BUSCO lineage or lineages that BUSCO should use to assess the completness of the assembly.
If no lineage are indicated (if the lines are comented out as in the below example), BUSCO will be set to [auto-lineage](https://busco.ezlab.org/busco_userguide.html#automated-lineage-selection)
If you indiciate a specific lineage, this lineage needs to be [downloaded](https://busco.ezlab.org/busco_userguide.html#offline) and located in the path section in nextflow.config
```
//Optional (if not indicated, autolineage for busco)
//        lineage                 = ""
//        lineage2                = "vertebrata_odb10"
//        lineage3                = "metazoa_odb10"
//        lineage4                = "eukaryota_odb10"
```

(optionnaly) - Include additional dataset (nanopore data such as ultra-long reads or long reads; Illumina short read data for polishing)
By default, this pipeline relies on PacBio HiFi reads and Hi-C data only. In some cases, nanopore data may be generated, in which case, polishing using Illumina short-read may be preferable. If such approach is done, you can indicate the nanopore and Illumina short-read data in this section
```
//ONT input
//	ont_fastq_1		= "${raw_data_path}/nanopore/"

//Illumina short reads input
//	illumina_SR_read1	= "${raw_data_path}/SR/"
//	illumina_SR_read2	= "${raw_data_path}/SR/"
```


(optionnaly) - To use another assembler than hifiasm, modify the method section in the nextflow.config file (Details indiacted in the nextflow.config file)
```
//Method
        assembly_method         = "hifiasm"	// 'hifiasm' or 'canu' of 'flye' or 'verkko'
	assembly_secondary_mode	= "pacbio"	// Depends on the assembly method selected, details in the following lines :
// With hifiasm : 'pacbio' (uses pacbio data only), 'pacbio+hic' (--h1 //--h2 : include Hi-C integration, requires Hi-C reads, VGP says that the output requires additional manual curation), 'pacbio+ont' (--ul : Ultra-long ONT integration), 'pacbio+ont+hic'
// With canu : 'hicanu' (-pacbio-hifi : uses HiFi data only), 'ont' (-nanopore : uses nanopore data only), 'clr' (-pacbio : for clr reads (lower quality than hifi))
// With flye : 'hifi' (--pacbio-hifi mode), 'ont' (--nano-raw mode), 'pacbio+ont', 'clr' (--pacbio-raw)
// With verkko : 'pacbio' (--hifi: uses HiFi data only), 'ont' (--nano : uses nanopore data only), 'pacbio+ont' (--hifi --nano)
	polishing_method	= "none"	// 'pilon' or 'none'
	purging_method		= "purge_dups"	//DO NOT MODIFY
	scaffolding_method	= "yahs"	// 'yahs' or 'salsa'
```

(optionnaly) -  If you want to run additional steps of the pipeline. 
Most of them have been set to 'no' by default as they require local installation of tools or databases.

generate a [Jupiter plot](https://github.com/JustinChu/JupiterPlot) to compare the obtained assembly to the assembly of a closely related specie using [Circos](http://circos.ca), specify the GenBank assembly number of the closely related specie.
For example, if we wanted to compare to the [California sea lion](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009762305.2/), we would indicate "GCA_009762305.2"

- To run some additional steps (Kraken,  

Launch the pipeline

```
nextflow run main.nf -profile singularity
```



## **Credits**

The pipeline was originnally written by [@scorreard](https://github.com/scorreard) with the help and input from :

- Members of the Jones lab (Canada's Michael Smith Genome Sciences Centre, Vancouver, Canada).

	- Special thanks to [@Glenn Chang](https://github.com/Glenn032787) for reviewing that repo.

- Members of the Earth Biogenome Project and other affiliated projects.
- Members of the nf-core / nextflow community.



## **Details on the test dataset**

The PacBio data is a subset of covid sequences obtained with this command lines :

```
wget https://downloads.pacbcloud.com/public/dataset/HiFiViral/Jan_2022/m64187e_211217_130958.hifi_reads.bam
samtools view -b m64187e_211217_130958.hifi_reads.bam -s 123.001 > subset_covid_hifi.bam
```

The Hi-C data was downloaded from one of the nf-core test dataset

```
wget https://github.com/nf-core/test-datasets/blob/modules/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz?raw=true
wget https://github.com/nf-core/test-datasets/blob/modules/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz?raw=true
```
