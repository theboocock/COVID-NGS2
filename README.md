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


1. Processing NGS genomes, performing QC, and integrating this data with the publicly avaliable genomes.  (Snakefile=Snakemake_ngs_covid19)
2. Proprocess public genomes from gissaid into useful file-formats for genomics such as VCF and phylogenetic trees. (Snakefile=Snakefile_proess_gisaid)


### Runinng Snakefile_ngs_covid19

The pipeline expects only a single space-delimited file as input. The file has 3 columns with no header line. The first is the sample nname, the second is the path to read1, and the third is the path to read2. If your only performed single ended sequecing make sure the 3rd column contains the string "NA".


Example file input 

```
SRR10902284 SRR10902284_1.fastq.gz NA
SRR10903401 SRR10903401_1.fastq.gz SRR10903401_2.fastq.gz
``` 

To run the pipeline  

```
snakemake -j 8 --config sample_list=<SAMPLE LIST IN> --snakefile <PATH TO Snakefile_ngs_covid19> 
```


### Snakefile_ngs_covid19 output files. 

### Description of the output files coming soon. 
 

### Running Snakefile_process_gisaid

First download the avaliable sequences from gisaid using the download button. 

```
snakemake -j 8 --config  gisaid_in=<GIS FASTA>  --snakefile Snakefile_process_gisaid

``` 

### Description of the output files coming soon. 
