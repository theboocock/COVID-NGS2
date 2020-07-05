#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)



a= read.delim("metadata.tsv", header=T, sep="\t",stringsAsFactors=F,colClasses=c("character"))
b = read.delim("../200623_NS500551_0874_AHCGH7AFX2/metadata.tsv", header=T, sep="\t",stringsAsFactors=F,colClasses=c("character"))
clin_in = read.csv("../../COVID19-NGS/ref/clin/clin_07_04.csv", header=T, stringsAsFactors=F,colClasses=c("character"))


print(b$ct)

print(dim(b))
aa=  merge(b, clin_in,by="uid")
aa$ct = aa$ct.x
idx_get=(a$date == "")
library(dplyr)
aaa = bind_rows(aa,a)

#aaa$sample_type[aaa$sample_type == "neb_old"] = "neb" 

idx_all = which(is.na(aaa$ct))
aaa$ct[idx_all] = clin_in$ct[(match(aaa$uid[idx_all], clin_in$uid))]
aaa$date[idx_all] = clin_in$date[(match(aaa$uid[idx_all], clin_in$uid))]
aaa$sample = 1:nrow(aaa)
write.table(aaa,row.names=F, file="metadata_all.tsv", quote=F, sep="\t")
#clin_in$uid[

