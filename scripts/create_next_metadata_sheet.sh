#!/bin/bash
#
# Create samplei list. 
### Script for creatingg the metadata table given the current setup, probably good enough for now.

START_NUMBER=$1
ls *R1*fastq.gz | awk -v pwd="$PWD" -v startnumber="$START_NUMBER" '{split($1,a,"\.fastq.gz");c=$1;b=gsub("R1","R2",c); printf "%s %s %s %03i\n",a[1],pwd"/"$1,pwd"/"c,NR + startnumber-1}'   > samples1.txt
folder=$(dirname $0)

### extract ### 
SAMPLE_SHEETS=$2
OUTPUT=$3
SAMPLE_TYPE=$4
while read line
do
    fastq_one=$(echo $line | awk '{print $2}')
    in_fastq=$(zcat $fastq_one | head -1 )
    run=$(echo $in_fastq | awk -F ":" '{print $2}')
    echo $run >> run_names.txt
done < samples1.txt
paste -d ' ' samples1.txt run_names.txt > samples2.txt
cat $SAMPLE_SHEETS |  grep "\[Data\]"  -A 1000000 | grep -v "\[Data\]" > ss.csv
echo Rscript $folder/create_next_metadata_sample_sheet.R samples2.txt $OUTPUT ss.csv $SAMPLE_TYPE
Rscript $folder/create_next_metadata_sample_sheet.R samples2.txt $OUTPUT ss.csv $SAMPLE_TYPE
rm ss.csv
rm samples1.txt
rm run_names.txt
rm samples2.txt
