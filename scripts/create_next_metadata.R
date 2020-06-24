#!/usr/bin/env Rscript 
## columns are
## c("strain","virus","gisaid_epi_gsl", "genbank_accession","date","region","country","division","location","region_exposure","division_exposure","segment","length","host","age","sex","originating_lab","submitting_lab","authors","url",title","date_submitted")


library(stringr)
header_line = c("strain","gisaid_epi_isl","genbank_accession","date",	"region",	"iso","country"	,"division","location"	,"length"	,"host",	"originating_lab",	"submitting_lab"	,"authors","url","update","input_id")

args <- commandArgs(trailingOnly=T)
file_in= args[1]

sample_in = read.table(file_in,header=F, stringsAsFactors=F, colClasses=rep("character",4))
### Final column

###

### TODO: check these arguments to ensure everything matches up
clin_in = read.csv(args[3], header=T)
clin_in$ID= toupper(clin_in$ID)
clin_in$ID. = NULL
uid = str_extract(clin_in$ID,"U[0-9]+")
uid = str_replace_all(uid,"U","")
uids = sprintf("U%04i",as.numeric(uid))
clin_in$ID = uids
sample_type = args[4] 
sample_sheet = args[5]
if(!is.na(sample_sheet)){
    print("Sample sheet found, skipping clinical inforamtion")

}

length_default = 29903
host_default = "Human"
originating_lab = "UCLA"
submitting_lab = "UCLA"
authors = "UCLA"
division = "California"
location = "Los Ageles County"
default_date = "2020-03-20"
region = "North America"

strain_name = sample_in$V4 
fasta_name = sample_in$V4 

df_meta = data.frame(sample=fasta_name, strain=strain_name, sample_type = sample_type, virus="ncov", gisaid_epi_isl="?", gennbak_acession="?",region=region, country="USA",
					division="California", location="Los Angeles", region_exposure=region,country_exposure="USA", division_exposure=division,segment="genome",length=length_default
				  ,host=host_default,age="?",sex="?",originating_lab="UCLA",submitting_lab="UCLA", authors="UCLA",url="NA",title="NA",date_submitted=default_date, lib_name=sample_in$V1, fastq_one =sample_in$V2, fastq_two = sample_in$V3,run_name=sample_in$V5)

sample_id = df_meta$lib_name
uid = str_extract(sample_id,"U[0-9]+")
uid = str_replace_all(uid,"U","")
uids = sprintf("U%04i",as.numeric(uid))
df_meta$uid = uids
print(dim(df_meta))
print(tail(df_meta))
print(tail(clin_in))
df_merged = merge(df_meta,clin_in,by.x="uid",by.y="ID")

write.table(df_merged, file=args[2], quote=F, row.names=F, sep="\t") 
