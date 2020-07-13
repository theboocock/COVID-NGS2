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
ARCTIC_LOCATIONS=${10}
READ_GROUP=${11}
SARS_REFERENCE=${12}
BAM_AMPLICON=${13}
BAM_UMI=${14}
echo $SAMPLE_TYPE
touch $SUMMARISE_MAPPING
touch $LOG
samtools view -H $INPUT_BAM > header.tmp
cat header.tmp | grep "@HD" > top.tmp
cat header.tmp | grep -A 20 "${CONTIG}" > bot.tmp
cat top.tmp bot.tmp > header.tmp 

# Remove duplicates, supplemental alignments, and filter on mapping quality
# https://broadinstitute.github.io/picard/explain-flags.html
samtools view  -F 3328 -q $MAPQ_FILT $INPUT_BAM $CONTIG > tmp.sam  2> /dev/null
cat header.tmp tmp.sam  | samtools view -hbS > tmp2.bam  2> /dev/null

### both reads must map to the same genome ## 

### So if read2 maps to human sars2 must map to non-humann ### 

samtools view -h tmp2.bam | awk 'BEGIN{OFS="\t"} {if(($7 == "=" || $0 ~ "@") && $7 != "*" ){print $0}}' |  samtools view -hb > tmp4.bam 2> /dev/null

DIR=$(dirname $0)
if [[ "${SAMPLE_TYPE}" =~ amp ]]; then
    $DIR/amplicon/check_read_one_position_filter.py -i tmp4.bam -o $BAM_AMPLICON -s $SUMMARISE_MAPPING -a $ARCTIC_LOCATIONS
    $DIR/amplicon/umi_fastaq.py -i $BAM_AMPLICON -o $OUTPUT_BAM --fastq-one-out out.fastq --fastq-two-out out2.fastq --read-group "${READ_GROUP}"  --reference $SARS_REFERENCE --output-summarise-molecular-info $BAM_UMI

elif [[ "${SAMPLE_TYPE}" =~ arctic ]]; then
    $DIR/amplicon/check_read_one_position_filter.py -i tmp4.bam -o $BAM_AMPLICON -s $SUMMARISE_MAPPING -a $ARCTIC_LOCATIONS
    $DIR/amplicon/umi_fastaq.py -i $BAM_AMPLICON -o $OUTPUT_BAM --fastq-one-out out.fastq --fastq-two-out out2.fastq --read-group "${READ_GROUP}"  --reference $SARS_REFERENCE --output-summarise-molecular-info $BAM_UMI
else
    cp tmp4.bam $OUTPUT_BAM
    samtools view  -q $MAPQ_FILT $INPUT_BAM $CONTIG > tmp.sam  2> /dev/null 
    cat header.tmp tmp.sam  | samtools view -hbS > tmp2.bam  2> /dev/null
    samtools view -h tmp2.bam | awk 'BEGIN{OFS="\t"} {if(($7 == "=" || $0 ~ "@") && $7 != "*" ){print $0}}' |  samtools view -hb > tmp4.bam 2> /dev/null
    cp tmp4.bam $BAM_AMPLICON
    touch $BAM_UMI
fi

samtools index $BAM_AMPLICON
samtools index $OUTPUT_BAM
