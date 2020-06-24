#!/bin/bash

OUTPUT_VCF=$1
LOG=$2
INPUT_FILES=$3

#echo $1
#echo "bcftools merge --missing-to-ref -l $INPUT_FILES | bgzip -c >  $OUTPUT_VCF"  > command.log
bcftools merge --missing-to-ref -l $INPUT_FILES | bgzip -c >  $OUTPUT_VCF 2> $LOG
tabix -p vcf $OUTPUT_VCF
