#!/bin/bash


java -Xmx2g $1 -jar $2 BuildBamIndex INPUT=$3  OUTPUT=$4  2> $5
