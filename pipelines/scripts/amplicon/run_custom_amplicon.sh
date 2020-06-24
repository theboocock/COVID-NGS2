#!/bin/bash


analysis_type=$1
read_one=$2
read_two=$3
read_one_out=$4
read_two_out=$5
oligo_pool=$6
arctic_primers=$7
folder=$(dirname $0)
log=$8

if [ "${analysis_type}" = "neb" ]
then
	ln -s $read_one $read_one_out
	ln -s $read_two $read_two_out 
elif [ "${analysis_type}" = "neb_old" ] 
then
	ln -s $read_one $read_one_out
	ln -s $read_two $read_two_out 
else
	touch $read_one_out 
	touch $read_two_out 
	$folder/custom_amplicon_parser.R $read_one $read_two $read_one_out $read_two_out $oligo_pool $analysis_type $arctic_primers  > $log 2>&1
fi
