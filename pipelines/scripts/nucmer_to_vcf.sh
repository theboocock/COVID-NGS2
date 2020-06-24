#!/bin/bash


NUCMER_TO_VCF=$1
#"{SCRIPTS_DIR}/nucmer_to_vcf {NUCMER_TO_VCF} {input} {output} {REFERENCE} {log}"
INPUT_SNP=$2
OUTPUT_VCF=$3
REFERENCE=$4
LOG=$5
#./my-mummer-2-vcf.py -g ../../ref/sars-cov-2.fasta --input-header --output-header -s ../../../ngs/test/align_public/outputs/snps/USA_NH_0008_2020.snp
echo $NUCMER_TO_VCF -g $REFERENCE --input-header --output-header -n -s $INPUT_SNP | bgzip -c > $OUTPUT_VCF 2> $LOG
$NUCMER_TO_VCF -g $REFERENCE --input-header --output-header -n -s $INPUT_SNP | bgzip -c > $OUTPUT_VCF 2> $LOG
tabix -p vcf $OUTPUT_VCF   
