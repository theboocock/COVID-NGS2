#!/bin/bash
#$ -N show-snps
#$ -cwd
#$ -l h_data=2G
#$ -l h_rt=01:00:00


INPUT_FILE=$1
INPUT_FULL=$(readlink -e $INPUT_FILE)
OUTPUT_FILE=$2
REFERENCE=$3
LOG=$4
CWD=$5
BASE_INPUT=$(basename $(head -1 $INPUT_FILE | awk '{print $2}'))
BASE_INPUT=$(echo $BASE_INPUT | awk '{split($1,a,"\.nucmer"); print a[1]}')
echo $BASE_INPUT
mkdir -p outputs/public_fasta/
OUTPUT_FOLDER=outputs/public_fasta/$BASE_INPUT

cp $CWD/outputs/public_fasta/$BASE_INPUT $OUTPUT_FOLDER

touch $LOG

cat $INPUT_FILE | sed "1s|^.*$|${REFERENCE} ${OUTPUT_FOLDER}|" > tmp.file
show-snps -T -C tmp.file > $OUTPUT_FILE 
