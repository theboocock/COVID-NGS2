#!/bin/bash
#$ -N run_minimap2
#$ -cwd
#$ -l h_data=2G
#$ -l h_rt=01:00:00

PUBLIC_FASTA=$1
EXTRACT_FASTA=$2
OUTPUT=$3
cat $PUBLIC_FASTA | bioawk -v seqname=$EXTRACT_FASTA -c 'fastx' '{if($name == seqname){print ">"$name"\n"$seq}}' > $OUTPUT
