#!/usr/bin/env Rscript
##
## Script process nextseq runs from Wednesday 24 June to extract the useful information for plotting in R.

args = commandArgs(trailing=T)
in_file = args[1]
output_file = args[2]
in_sheet = read.csv(in_file, header=T)
library(stringr)
uid_split = unlist(lapply(str_split(in_sheet$uid,"_"), function(x){x[1]}))
is_diluted = unlist(lapply(str_split(in_sheet$uid,"_"), function(x){x[2] == "10di"}))
is_diluted[is.na(is_diluted)] = F
in_sheet$diluted = is_diluted
in_sheet$uid = uid_split
in_sheet$ct = unlist(lapply(str_split(in_sheet$ct,"_"), function(x){x[1]}))
# parse experimental design column
split_prep = str_split(in_sheet$prep,"_")
ss_type = unlist(lapply(split_prep, function(x){x[2]}))
in_sheet$nsr=ss_type
ss_type = unlist(lapply(split_prep, function(x){x[1]}))
oligo_type = unlist(lapply(str_split(ss_type,"\\+"), function(x){print(x);x=paste(x,collapse="|",sep=",");return(x)}))
in_sheet$oligos_used = oligo_type
write.csv(in_sheet, file=output_file,quote=F,row.names=F)
