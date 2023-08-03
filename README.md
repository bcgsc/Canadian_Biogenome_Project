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

The pipeline also require information related to the specie of interest such as genome size or ploidy. This information can be found on [GoaT](https://goat.genomehubs.org).



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

- [Busco](https://busco.ezlab.org/busco_userguide.html#download-and-automated-update)
- [Kraken](http://ccb.jhu.edu/software/kraken/)

<p align="center">
    <img title="The Canadian Biogenome Project Workflow" src="res/CBP_workflow.png" width=50%>
</p>
<p align="center">
Figure : Overview of the Canadian Biogenome project assembly pipeline
</p>


## Running the pipeline with test data (will work once the repo is public)
To run this pipeline, you need nextflow, conda and singularity installed on your system.

A set of test data are available in this repo to allow you to test the pipeline with just one command line:

```
nextflow run bcgsc/Canadian_Biogenome_Project -latest -r V2
```

The outputs are organized in several subfolder that are self-explenatory.



## Running the pipeline with your own data
Clone the repository in your local environment:

```
git clone https://github.com/bcgsc/Canadian_Biogenome_Project.git
cd Canadian_Biogenome_Project
```

Modify the nextflow.config file:

- Indicate the location of the input file
- Indicate the required information
	

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
