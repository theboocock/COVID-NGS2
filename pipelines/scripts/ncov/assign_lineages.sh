#!/bin/bash
#
# @author James Boocock
#

in_fasta=$1
nextstrain_clades_out=$2
gisaid_clades_out=$3
gisaid_clades_ref=$4
nextstrain_clades_ref=$5
genbank_ref_sars2=$6
threads=$7
##TODO: assign clades
SCRIPT_DIR=$(dirname $0)
$SCRIPT_DIR/assign_clades.py --sequences $in_fasta --output $gisaid_clades_out  --clade-definition $gisaid_clades_ref --reference-input $genbank_ref_sars2 --nthreads $threads
$SCRIPT_DIR/assign_clades.py --sequences $in_fasta --output $nextstrain_clades_out --clade-definition $nextstrain_clades_ref --reference-input $genbank_ref_sars2 --nthreads $threads

