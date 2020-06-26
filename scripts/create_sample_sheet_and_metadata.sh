#!/bin/bash

SCRIPTS_DIR=$(dirname $0)
echo $SCRIPTS_DIR
SAMPLE_SHEET=$1
DATE=$2
OUTPUT_SAMPLE_SHEET=$3
START_INDEX=$4
Rscript /u/scratch/s/smilefre/covidngs/COVID19-NGS/scripts/extract_conditions_from_new_sheet.R $SAMPLE_SHEET sample_sheet.csv 
python $SCRIPTS_DIR/create_sample_sheet.py -s sample_sheet.csv -d $DATE -e 48_neb --output-file-neb  sample_sheet_neb.csv --output-file-amp  sample_sheet_amp.csv
cd amp 
bash $SCRATCH/covidngs/COVID19-NGS/scripts/create_sample_list.sh $START_INDEX /u/scratch/s/smilefre/covidngs/COVID19-NGS/ref/clin/clin_05_31.csv metadata.tsv amp_complicated  ../sample_sheet.csv 


