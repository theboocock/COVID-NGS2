#!/usr/bin/env Rscript

a = read.delim("metadata_all.tsv",sep="\t")
b = read.delim("../200702_NS500551_0877_AHCFTWAFX2/metadata.tsv", header=T, sep="\t")
library(dplyr)
b$ct = b$ct.x
library(stringr)
a$lib_name_short = sapply(str_split(a$lib_name,"_"), function(x){x[[1]][1]}) 
aa = bind_rows(a,b)
aa$sample = 1:nrow(aa)
write.table(aa,file="metadata_complete_july4.csv",quote=F,row.names=F,sep="\t")
