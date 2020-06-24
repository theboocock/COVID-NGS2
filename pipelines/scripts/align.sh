#!/bin/bash
#$ -N run_minimap2
#$ -cwd
#$ -l h_data=16G
#$ -l h_rt=01:00:00

RG=$1
REF=$2
READ_ONE=$3
READ_TWO=$4
OUTPUT=$5
LOG=$6
THREADS=$7


if [[ "$(basename ${READ_TWO})" == "NA" ]]
then
    bwa mem -M -t ${THREADS} -R "${RG}" $REF $READ_ONE 2> $LOG | samtools view -Sb - > $OUTPUT 
else
    bwa mem -M -t ${THREADS} -R "${RG}" $REF $READ_ONE $READ_TWO 2> $LOG | samtools view -Sb - > $OUTPUT 
fi
