#!/bin/bash
#$ -l h_data=32G,highp, -pe shared 1 
#$ -l time=124:00:00
#$ -cwd 
#$ -M theboocock@ucla.edu
#$ -m a
#$ -e logs/
#$ -o logs/


INPUT_BAM=$1
OUTPUT_BAM=$2
CONTIG=$3
MAPQ_FILT=$4
LOG=$5
PICARD=$6
TMP=$7
samtools view -H $INPUT_BAM > header.tmp
cat header.tmp | grep "@HD" > top.tmp
cat header.tmp | grep -A 20 "${CONTIG}" > bot.tmp
cat top.tmp bot.tmp > header.tmp 
echo "Done!" > $LOG
samtools view  -q $MAPQ_FILT $INPUT_BAM $CONTIG > tmp.sam 
cat header.tmp tmp.sam  | samtools view -hbS > tmp2.bam

samtools view -h tmp2.bam | awk 'BEGIN{OFS="\t"} {if($7 != "*"){print $0}}' |  samtools view -hb > tmp3.bam
java -Xmx4g $TMP -jar $PICARD FixMateInformation I=tmp3.bam O=$OUTPUT_BAM
samtools index $OUTPUT_BAM
