#!/bin/bash
#$ -l h_data=32G,highp, -pe shared 1 
#$ -l time=124:00:00
#$ -cwd 
#$ -M theboocock@ucla.edu
#$ -m a
#$ -e logs/
#$ -o logs/


INPUT=$1
OUTPUT_VCF_NO_SNPS=$2
OUTPUT_VCF_INDELS=$3
LOG=$4
MIN_QUAL=$5
MIN_DEPTH=$6

#### Biallelic sites only 
bcftools view -v snps -m1 -M2 -e " QUAL < ${MIN_QUAL} "  $INPUT   | bgzip -c > $OUTPUT_VCF_NO_SNPS
bcftools view -e " QUAL < ${MIN_QUAL} ||  MIN(FORMAT/DP) < ${MIN_DEPTH} "  $INPUT   | bgzip -c > $OUTPUT_VCF_INDELS
tabix -p vcf $OUTPUT_VCF_NO_SNPS
tabix -p vcf $OUTPUT_VCF_INDELS
