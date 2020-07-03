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
SAMPLE_TYPE=$8
SUMMARISE_MAPPING=$9
echo $SAMPLE_TYPE
touch $SUMMARISE_MAPPING
touch $LOG
samtools view -H $INPUT_BAM > header.tmp
cat header.tmp | grep "@HD" > top.tmp
cat header.tmp | grep -A 20 "${CONTIG}" > bot.tmp
cat top.tmp bot.tmp > header.tmp 
samtools view  -q $MAPQ_FILT $INPUT_BAM $CONTIG > tmp.sam  2> /dev/null
cat header.tmp tmp.sam  | samtools view -hbS > tmp2.bam  2> /dev/null

samtools view -h tmp2.bam | awk 'BEGIN{OFS="\t"} {if(($7 == "=" || $0 ~ "@") && $7 != "*" ){print $0}}' |  samtools view -hb > tmp4.bam 2> /dev/null

DIR=$(dirname $0)
if [[ "${SAMPLE_TYPE}" =~ amp ]]; then
    $DIR/amplicon/check_read_one_position_filter.py -i tmp4.bam -o $OUTPUT_BAM -s $SUMMARISE_MAPPING
else
    cp tmp4.bam $OUTPUT_BAM    
#echo ${SAMPLE_TYPE}
    #echo "HEREHREH"
fi


samtools index $OUTPUT_BAM
