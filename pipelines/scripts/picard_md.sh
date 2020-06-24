#!/bin/bash


java -Xmx2g $1 -jar $2 MarkDuplicates INPUT=$3  OUTPUT=$4  METRICS_FILE=$5 READ_NAME_REGEX=null VALIDATION_STRINGENCY=LENIENT > $6 
