#!/usr/bin/env bash
#
#

INPUT_FASTA=$1
OUTPUT_LINEAGES=$2
THREADS=$3
pangolin -t  $THREADS $1
cp lineage_report.csv $2
