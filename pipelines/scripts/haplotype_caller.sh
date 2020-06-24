#!/bin/bash



$1 \
HaplotypeCaller \
-R $2 \
-I $3 \
-ERC GVCF \
--native-pair-hmm-threads $5 \
-O $4  2> $5
