#!/bin/bash
#$ -l h_data=32G,highp, -pe shared 1 
#$ -l time=124:00:00
#$ -cwd 
#$ -M theboocock@ucla.edu
#$ -m a
#$ -e logs/
#$ -o logs/



intervals=$(sed -n ${SGE_TASK_ID}p /u/scratch/s/smilefre/vcfs/genomes.txt)
intervals=$(echo $intervals | awk '{print $1}')
echo $intervals
GATK_PATH=/u/home/s/smilefre/scripts/vcf/gatk-4.1.4.1/gatk
REF=/u/home/s/smilefre/project-kruglyak/ref/saccer3_with_alt_chin.fasta
for file in /u/scratch/s/smilefre/vcfs/gvcfs_all/*.g.vcf
do
    in_string="${in_string} -V ${file} "
done
$GATK_PATH GenomicsDBImport \
$in_string \
--genomicsdb-workspace-path combine_gvcfs_$intervals \
--intervals $intervals

$GATK_PATH GenotypeGVCFs \
-R $REF \
-V gendb://combine_gvcfs_$intervals \
-O ${intervals}.vcf.gz
rm -Rf combine_gvcfs_$intervals
#   exit 1


