#!/bin/bash


GATK_PATH=$1
INPUT=$2
INPUT=$(echo $INPUT| sed "s/_complete.txt//g")
echo $INPUT
OUTPUT=$3
LOG=$4
REF=$5

$GATK_PATH GenotypeGVCFs \
    -R $REF \
    -V gendb://${INPUT} \
    -O $OUTPUT 2> $LOG

