# COVID19 NGS pipeline 

## Introduction

Processing pipeline for WGS of COVID19 genomes sequenced. 

## Installation

Clone this repository using

```
    git clone git@github.com:theboocock/COVID19-NGS.git
```

Install dependenncies via Conda by running: 

```
    conda env create -f envs/ngs_covid19.yaml
```

Once installed, the environment needs to be activated by running:

```
    conda activate ngs_covid19
```

To download external dependencies, and reference databases also run `./check_install.sh` after activating the conda environnmennt
from the base directory in which the pipeline is installed. 

## Pipelines

This package consists of different pipelines for processing COVID19 WGS data.  

To date there are two main pipelines 


1. Processing NGS genomes, performing QC, and integrating this data with the publicly avaliable genomes.  (Snakefile=Snakemake_ngs_covid19_dev)
2. Proprocess public genomes from gissaid into useful file-formats for genomics such as VCF and phylogenetic trees. (Snakefile=Snakefile_proess_gisaid)


### Runinng Snakefile_ngs_covid19

The pipeline expects only a single space-delimited file as input. The file has 3 columns with no header line. The first is the sample nname, the second is the path to read1, and the third is the path to read2. If your only performed single ended sequecing make sure the 3rd column contains the string "NA".


Example file input 

```
UID sample  fastq_one fastq_two lib_type uid_sample_type
U0001   1   1_R1_fastq.gz  1_R2_fastq.gz    neb U0001_neb
U0002   2   2_R1_fastq.gz   2_R2_fastq.gz   amp U0002_amp
``` 

To run the pipeline  

```
snakemake -j 8 --use-conda --config merge_rule="column_uid_sample_type" sample_list=<SAMPLE LIST IN> --snakefile <PATH TO Snakefile_ngs_covid19_dev> 
```


### Snakefile_ngs_covid19 output files. 

### Description of the output files coming soon. 
 

### Running Snakefile_process_gisaid

First download the avaliable sequences from gisaid using the download button. 

```
snakemake -j 8 --config  gisaid_in=<GIS FASTA>  --snakefile Snakefile_process_gisaid

``` 

### Description of the output files coming soon. 
