#!/bin/bash

$1 \
GenomicsDBImport \
-R $2 \
-I $3 \
-ERC GVCF \
--native-pair-hmm-threads $5 \
-O $4 
