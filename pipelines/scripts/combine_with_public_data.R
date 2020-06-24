#!/usr/bin/env Rscript

library(tidyverse)
args = commandArgs(trailingOnly=T)
metadata_input = read.delim(args[1],header=T,sep="\t")
gisaid_input = read.delim(args[2],header=T,sep="\t")
out_file = args[3] 
metadata_input$genbank_accession= "?"
colnames(metadata_input)[1] = "strain"
out_df = bind_rows(metadata_input,gisaid_input)
write.table(out_df, file=out_file, sep="\t",col.names=T, row.names=F, quote=F)
