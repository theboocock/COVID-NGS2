#!/bin/bash


analysis_type=$1
read_one=$2
read_two=$3
read_one_out=$4
read_two_out=$5
oligo_pool=$6
arctic_primers=$7
all_primers=$8
primer_plates=$9
folder=$(dirname $0)
log=$10

if [ "${analysis_type}" = "neb" ]
then
	cp $read_one $read_one_out
	cp $read_two $read_two_out 
elif [ "${analysis_type}" = "neb_old" ] 
then
	cp $read_one $read_one_out
	cp $read_two $read_two_out 
else
	touch $read_one_out 
	touch $read_two_out 
	$folder/custom_amplicon_parser.R $read_one $read_two $read_one_out $read_two_out $oligo_pool $analysis_type $arctic_primers $all_primers $primer_plates #> $log 2>&1 
    #> $log 2>&1
fi

