#!/bin/bash

RG=$1
REF=$2
READ_ONE=$3
OUTPUT=$4

bwa mem -t 1 -R "${RG}" $REF $READ_ONE | samtools sort > $OUTPUT
samtools index $OUTPUT
