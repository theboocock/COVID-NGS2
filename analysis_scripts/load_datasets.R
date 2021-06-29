library(lubridate)
library(tidyverse)
library(data.table)
### Load in the QC input file ###
library(Biostrings)
library(stringr)
library(glue)
#source("R/load_datasets.R")
qc_input2 = read.delim("data/qc_report.tsv")
bins_ct = cut(qc_input2$ct,breaks=c(0,15,20,25,40))
qc_input2$bin_ct= bins_ct
qc_input2$library_type = sapply(str_split(qc_input2$id_library_type,"_"),function(x){x[2]})
qc_input2$pass = qc_input2$coverage_3x > 0.8
qc_input2$date_fix = mdy(qc_input2$date)
qc_input2$week = week(qc_input2$date_fix)
qc_input2$lib_type = "V-seq"
qc_input2$lib_type[qc_input2$library_type == "neb"] = "NEB"
qc_input2 = qc_input2 %>% mutate(qc_lineage=ifelse(is.na(status),"failed","passed_qc"))
qc_input2 %>% pivot_wider(uid,names_from=library_type,names_sep =".",values_from=c(lineage,qc_lineage,probability)) %>% filter(!is.na(lineage.amp) &  !is.na(lineage.neb))
#### 100% concordance in SNPs so far. 
qc_input2 %>% ggplot(aes(x=lineage)) + geom_bar()
qc_lineages = qc_input2 %>% filter(pass) %>% filter(!duplicated(uid)) 
library(lubridate)
qc_lineages$date_fix = mdy(qc_lineages$date)
qc_lineages$week =  week(qc_lineages$date_fix)
qc_lineages$month = month(qc_lineages$date_fix)

qc_lineages$strain = qc_lineages$sample_name_fasta
qc_lineages$pangolin_lineage = qc_lineages$lineage
phased_genomes = readDNAStringSet("data/all_3_none.fasta")
phased_genomes_curate = readDNAStringSet("data/phased.fasta")
final_genomes = readDNAStringSet("data/all.fasta")
names_phased = names(phased_genomes) 
phased_genome_names = unlist(lapply(str_split(names_phased,":"), function(x){x[1]}))
name_final_genomes = str_replace_all(names(final_genomes),"-","_")
qc_lineages$g1 = unlist(lapply(str_split(qc_lineages$pangolin_phased,","), function(x){x[1]}))
qc_lineages$g2 = unlist(lapply(str_split(qc_lineages$pangolin_phased,","), function(x){x[2]}))
string_set = c()
qc_lineages_phased = data.frame()
string_set_2 = c()
phased_genome_names_curate = str_replace_all(names(phased_genomes_curate),"_[0-9]$","")
two_lineages = c()
for(i in 1:nrow(qc_lineages)){
  row = qc_lineages[i,]
  if(row$g1 != row$g2){
    string_tmp = phased_genomes[phased_genome_names == row$merged_id]
    names(string_tmp) = str_replace_all(names(string_tmp),":","_")
    string_set = c(string_set, string_tmp)
    tmp_row = rbind(row,row)
    tmp_row$lineage[1] = row$g1
    tmp_row$lineage[2] = row$g2
    tmp_row$strain = paste0(tmp_row$merged_id,":",1:2,sep="")
    
    string_tmp_all = phased_genomes_curate[phased_genome_names_curate == row$merged_id]
    string_set_2 = c(string_set_2,string_tmp_all)
    two_lineages = c(two_lineages, string_tmp_all)
    
    #tmp_row$strain[2] =  paste0(tmp_row$merged_id,":2",sep="")
  }
  else{
    string_set = c(string_set, final_genomes[name_final_genomes == row$merged_id])
    string_set_2 = c(string_set_2, final_genomes[name_final_genomes == row$merged_id])
    
    tmp_row = row
    tmp_row$strain = tmp_row$sample_name_fasta
    
  }
  qc_lineages_phased = rbind(qc_lineages_phased, tmp_row)
}
keep_genomes = do.call("c",string_set)
keep_genomes_2 = do.call("c",string_set_2)

two_lineages = do.call("c",two_lineages)
writeXStringSet(two_lineages, filepath = "data/phased_genomes_two_lineages.fasta")
xstringset = readDNAStringSet(("data/gisaid/sequences_2020-07-10_23-53.fasta"))
library(data.table)
gisaid = fread("data/gisaid/metadata_2020-07-10_23-53.tsv")
new_pango_gisaid = read.table("data/lineages_merged.csv",sep=",",header=T)
gisaid = gisaid %>% left_join(new_pango_gisaid,by=c("strain"="taxon"))
gisaid$date_fix= ymd(gisaid$date)
gisaid$pangolin_lineage = gisaid$lineage
gisaid$probability = gisaid$probability

#qc_input2[!(qc_input2$merged_id %in% bad_samples$merged_id),]

bad_samples = read.csv("data/long_sample_list_august19th.csv")
bad_sample_names = c(str_replace_all(bad_samples$merged_id,"_","-"),paste(bad_samples$merged_id,"_2",sep=""),paste(bad_samples$merged_id,"_1",sep=""))
duplicate_sent_genomes  = read.csv("data/Duplicate Sent Genomes_Deidentified.csv")
dd = duplicate_sent_genomes %>% inner_join(qc_lineages,by=c("id"="uid")) %>% group_by(group) %>% filter(pass) %>% mutate(n=n())  %>% filter(n >=2) 

dd = dd %>% select(merged_id,id,group)
g1 = (qc_lineages %>% left_join(dd,by=c("merged_id"="merged_id")) %>% group_by(group) %>% mutate(n2=sum(!is.na(group)))) %>% filter(n2 > 1) %>% filter(row_number() == 1) %>% ungroup
g2 = (qc_lineages %>% left_join(dd,by=c("merged_id"="merged_id")) %>% group_by(group) %>% mutate(n2=sum(!is.na(group)))) %>% filter(n2== 0) %>% ungroup
qc_lineages_phased = bind_rows(g1,g2)
qc_lineages_phased$pangolin_lineage = qc_lineages_phased$lineage
qc_lineages_phased$genbank_accession = "?"
qc_lineages_phased$paper_url = "?"
#qc_lineages_phased$strain = qc_lineages_phased$sample_name_fasta
qc_lineages_phased$location = "Los Angeles County"
qc_lineages_phased$strain = str_replace_all(qc_lineages_phased$strain,":","_")
#colnames(gisaid)[!colnames(gisaid) %in% colnames(qc_lineages_phased)]
qc_lineages_phased = qc_lineages_phased %>% filter(!(strain %in% bad_sample_names))
qc_lineages_phased$date_submitted = ymd(qc_lineages_phased$date_submitted)
gisaid_phased = qc_lineages_phased %>% select(colnames(qc_lineages_phased)[colnames(qc_lineages_phased) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
gisaid_phased$date = gisaid_phased$date_fix
gisaid_phased$strain = str_replace_all(gisaid_phased$strain,":","_")
gisaid_phased$la_county = ifelse(gisaid_phased$location == "Los Angeles County","Los Angeles County","Other")
#gisaid_phased = gisaid_phased[gisaid_phased$probability > .8,]

#write.table(as.data.frame(gisaid_phased),file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased.tsv",col.names=T,row.names=F,quote=F,sep="\t")
#aa = c(keep_genomes[qc_lineages_phased$sample_name_fasta],xstringset)
#writeXStringSet(aa,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased.fasta")
##### 
pangolin_new_lineages =read.csv("data/lineages.csv")
qc_lineages_phased2 = qc_lineages_phased %>% left_join(pangolin_new_lineages,by=c("strain"="taxon")) %>% mutate(pangolin_lineage= ifelse(is.na(lineage.y),pangolin_lineage,lineage.y),probability= ifelse(is.na(lineage.y),probability.x,probability.y))
#qc_lineages_phased2 %>%rename(pangolin_lineage,pcox)
gisaid_phased = qc_lineages_phased2 %>% select(colnames(qc_lineages_phased2)[colnames(qc_lineages_phased2) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
gisaid_phased$date = gisaid_phased$date_fix
gisaid_phased$strain = str_replace_all(gisaid_phased$strain,":","_")
gisaid_phased = gisaid_phased %>% filter(!(strain %in% bad_sample_names))
gisaid_phased$la_county = ifelse(gisaid_phased$location == "Los Angeles County","Los Angeles County","Other")

# Current run
# TODO: get the list 
#current_run = read.table("~/Dropbox/COVID19/ncov_new/again/ncov/results/north-america_usa_los-angeles_update_hets/current_sub.txt",header=F)

#gisaid_phased = gisaid_phased %>% filter(strain %in% current_run$V1)

#write.table(gisaid_phased,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased_curate.tsv",col.names=T,row.names=F,quote=F,sep="\t")
aa = c(keep_genomes_2,xstringset)
#aa = aa[gisaid_phased$strain]
#writeXStringSet(aa,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased_curate.fasta")

# switch lineages ##
subset_genomes = final_genomes[qc_lineages$sample_name_fasta]
qc_lineages$pangolin_lineage = qc_lineages$lineage
qc_lineages$genbank_accession = "?"
qc_lineages$paper_url = "?"
qc_lineages$strain = qc_lineages$sample_name_fasta
qc_lineages$location = "Los Angeles County"
colnames(gisaid)[!colnames(gisaid) %in% colnames(qc_lineages)]
#qc_lineages
qc_lineages$date_submitted = ymd(qc_lineages$date_submitted)
gisaid_local = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
aa = c(subset_genomes,xstringset)
gisaid_local$date = gisaid_local$date_fix
#write.table(gisaid_local,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final.tsv",col.names=T,row.names=F,quote=F,sep="\t")
#writeXStringSet(aa,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final.fasta")
### Only LA with Wuhan and LA ####
#
#calif_wuhan = gisaid %>% filter(location == "Los Angeles County" | grepl("Hu-1",strain))
#out_seq = c(subset_genomes,xstringset[calif_wuhan$strain])
#gisaid_la = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid)]) 
#write.table(gisaid_la ,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la.tsv",col.names=T,row.names=F,quote=F,sep="\t")
#writeXStringSet(out_seq,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la.fasta")


#writeXStringSet(xstringset[aaa$strain], file="~/Dropbox/COVID19/analysis_for_epi_paper/data/la_icahn.fasta")
#writeXStringSet(keep_genomes,) ###
### For any sample with both lineages write out corresponding lines in the output fasta ###
### For other samples just output their consensus ### 

### Create a big merged one with gisaid  #### 
### And a unmerged one of just LA genomes ####

### Now look at the phased genomes ### 

#df_phased %>% left_join(qc_lineages,by=c("sample_id"="sample_name_fasta"))
df_phased =data.frame(sample_id=str_replace_all(unlist(lapply(str_split(names(keep_genomes),":"), function(x){x[1]})),"_","-"),fasta_name=names(keep_genomes)) %>% left_join(qc_lineages,by=c("sample_id"="sample_name_fasta"))
df_phased$strain = df_phased$fasta_name
calif_wuhan = gisaid %>% filter(grepl("Hu-1",strain))
out_seq = c(keep_genomes,xstringset[calif_wuhan$strain])
gisaid_la_hets = df_phased %>% select(colnames(df_phased)[colnames(df_phased) %in% colnames(gisaid)]) %>% bind_rows(calif_wuhan)
#write.table(gisaid_la_hets ,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la._het.tsv",col.names=T,row.names=F,quote=F,sep="\t")
#writeXStringSet(out_seq,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la_hets.fasta")


#### Het phylogenies ###





gisaid = fread("data/gisaid/metadata_2020-07-10_23-53.tsv")
#gisaid = gisaid %>% left_join(california,by=c("strain"="seqName"))
gisaid$date_fix = ymd(gisaid$date)
new_pango_gisaid = read.table("data/lineages_merged.csv",sep=",",header=T)
gisaid$week = week(gisaid$date_fix)
gisaid_strain = gisaid %>% left_join(new_pango_gisaid,by=c("strain"="taxon"))
gisaid_strain$pangolin_lineage = gisaid_strain$lineage
good_genomes  = read.table("data/gisaid/good_genomes.txt",header=F)
gisaid_strain = gisaid_strain %>% inner_join(good_genomes,by=c("strain"="V1"))
california = gisaid_strain %>% select(strain, pangolin_lineage,date,date_fix, division, country,location,GISAID_clade) %>% filter(division == "California")
california$date = ymd(california$date)
california$week = week(california$date)
qc_lineages_filt = qc_lineages %>% select(strain, date_fix,pangolin_lineage,week)
la=  california %>%  filter(location == "Los Angeles County") %>% select(strain,date_fix,pangolin_lineage,week)
all_la = bind_rows(la,qc_lineages_filt)
need_colors = names(which(table(all_la$pangolin_lineage) >= 5))
color_map = viridisLite::viridis(length(need_colors) + 1)
color_map = RColorBrewer::brewer.pal(length(need_colors) + 1  , name ="Set3")
names(color_map) = c(need_colors,"other")
all_la$pango_new = ifelse(all_la$pangolin_lineage %in% names(color_map),all_la$pangolin_lineage,"other")
gisaid_strain$pango_new =  ifelse(gisaid_strain$pangolin_lineage %in% names(color_map),gisaid_strain$pangolin_lineage,"other")
gisaid_local = gisaid_local %>% filter(strain != "USA/CA-ALSR-0513-SAN/2020")

qc_lineages$pangolin_lineage = qc_lineages$lineage
qc_lineages$genbank_accession = "?"
qc_lineages$paper_url = "?"
qc_lineages$strain = qc_lineages$sample_name_fasta
qc_lineages$location = "Los Angeles County"
#colnames(gis)[!colnames(gisaid) %in% colnames(qc_lineages)]
gisaid_local = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid_strain)]) %>% bind_rows(gisaid_strain)
### Remove this sample with the incorrect date
gisaid_local = gisaid_local %>% filter(strain != "USA/CA-ALSR-0513-SAN/2020")
gisaid_phased_filt = gisaid_phased %>% filter(strain != "USA/CA-ALSR-0513-SAN/2020")
gisaid_phased_filt = gisaid_phased_filt %>% filter(strain %in% good_genomes$V1)
gisaid_phased_filt = gisaid_phased_filt %>% mutate(location=ifelse(location == "San Diego County","San Diego",location))
gisaid_phased_filt = gisaid_phased_filt %>% mutate(location=ifelse(location == "San Francisco County","San Francisco",location))
gisaid_phased_filt = gisaid_phased_filt %>% mutate(month=month(date_fix),week=week(date_fix))
#gisaid_local = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
lineages_manual = read.csv("data/lineages.csv")
#lineages_manual %>% write.table(file="data/S4.tsv",col.names=T,row.names=F,quote=F)
lineages_manual  =lineages_manual %>% mutate(sample=str_replace_all(str_replace_all(taxon,"_[0-2]$",""),"_","-")) %>% group_by(sample) %>% summarise(g1=lineage[1],g2=lineage[2])
lineages_manual$ordered_lineages = apply(cbind(lineages_manual$g1,lineages_manual$g2), 1, function(x){paste(sort(x),sep=",",collapse = ",")})
## ##
