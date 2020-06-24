#!/bin/bash

MEM=$7

java -Xmx${MEM}m $1 -jar $2 SortSam INPUT=$3 OUTPUT=$4 $5 2> $6 
