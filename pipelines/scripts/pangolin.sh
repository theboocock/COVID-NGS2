#!/usr/bin/env bash
#
#

INPUT_FASTA=$1
OUTPUT_LINEAGES=$2
THREADS=$3
PANGOLEARN=$4
if [[ "${PANGOLEARN}" == "TRUE" ]]
then
    pangolin -t  $THREADS $1 
else
    pangolin -t  $THREADS $1 --legacy 
fi
cp lineage_report.csv $2
