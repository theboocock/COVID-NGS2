#!/bin/bash

TMPDIR=$!
PICARD=$2
BAM_IN=$3
OUTPUT=$4
HISTOGRAM=$5

java -Xmx4g $TMPDIR -jar $PICARD CollectInsertSizeMetrics I=$BAM_IN O=$OUTPUT H=$HISTOGRAM M=0.5

