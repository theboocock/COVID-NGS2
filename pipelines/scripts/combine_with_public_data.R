#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=T)
metadata_input = read.delim(args[1],header=T,sep="\t",stringsAsFactors=F)
gisaid_input = fread(args[2],header=T,sep="\t")

keep_columns = c("strain","virus","gisaid_epi_isl","genbank_accession","date","region","country","division","location","region_exposure",	
"country_exposure", "division_exposure", "segment", "length" ,"host" ,"age" ,"sex", "pangolin_lineage", "GISAID_clade", "originating_lab", "submitting_lab", "authors" ,"url", "title", "paper_url", "date_submitted")

library(lubridate)

metadata_input$date = as.character(mdy(metadata_input$date))
### TODO: ###

### Remove sequences ####
metadata_input$strain = metadata_input$sample_name_fasta
metadata_input$location = "Los Angeles County"
metadata_input = metadata_input[!duplicated(metadata_input$uid),]
print(head(metadata_input$strain))
###

metadata_input = metadata_input[,colnames(metadata_input) %in% keep_columns]


out_file = args[3] 
out_ncov_la_metadat = args[4]
metadata_input$genbank_accession= "?"
colnames(metadata_input)[1] = "strain"
out_df = bind_rows(metadata_input,gisaid_input)
write.table(out_df, file=out_file, sep="\t",col.names=T, row.names=F, quote=F)
