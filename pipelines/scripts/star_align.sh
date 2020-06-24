#!/bin/bash


THREADS=$1
READ_ONE=$2
READ_TWO=$3
GENOME_DIR=$4
READ_GROUP=$5
OUTPUT=$6
LOG=$7

STAR --runThreadN $THREADS --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn $READ_ONE $READ_TWO --outSAMtype BAM Unsorted --outSAMattrRGline ${READ_GROUP} > $LOG
mv Aligned.out.bam $OUTPUT 
