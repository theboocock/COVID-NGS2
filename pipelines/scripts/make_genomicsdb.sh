#!/bin/bash
#$ -l h_data=32G,highp, -pe shared 1 
#$ -l time=124:00:00
#$ -cwd 
#$ -M theboocock@ucla.edu
#$ -m a
#$ -e logs/
#$ -o logs/



GATK_PATH=$1
REF=$2
# Strip /complete.txt
OUTPUT=$3
OUTPUT=$(echo $OUTPUT | sed "s/_complete.txt//g")
echo $OUTPUT
INTERVALS=$4

for file in ${@:5} 
do
    in_string="${in_string} -V ${file} "
done
echo $in_string
$GATK_PATH GenomicsDBImport \
$in_string \
--genomicsdb-workspace-path $OUTPUT/${intervals} \
--intervals $INTERVALS

date > ${OUTPUT}_complete.txt

##$GATK_PATH GenotypeGVCFs \
#-R $REF \
#-V gendb://combine_gvcfs_$intervals \
#-O ${intervals}.vcf.gz
#rm -Rf combine_gvcfs_$intervals
###   exit 1


