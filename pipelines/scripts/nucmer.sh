#!/bin/bash
#$ -N run_minimap2
#$ -cwd
#$ -l h_data=2G
#$ -l h_rt=01:00:00


REFERENCE=$1
INPUT_FILE=$2
OUTPUT_FILE=$3
LOG=$4


echo "nucmer -maxmatch -c 100 -p nucmer $REFERENCE $INPUT_FILE  > $LOG 2>&1"
nucmer -maxmatch -c 100 -p nucmer $REFERENCE $INPUT_FILE  > $LOG 2>&1
mv nucmer.delta $OUTPUT_FILE

