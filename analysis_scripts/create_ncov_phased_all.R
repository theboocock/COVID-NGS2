library(Biostrings)
library(stringr)
source("R/load_datasets.R")
phased_genomes = readDNAStringSet("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/outputs/consensus/merged/phased/sars2/all_3_none.fasta")
phased_genomes_curate = readDNAStringSet("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/phased/consensus/all.fasta")
final_genomes = readDNAStringSet("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/outputs/final/merged/all.fasta")
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
new_pango_gisaid = read.table("~/KRUGLYAK/covid_runs/rerun_all_new_july20/pangolin_gisaid/outputs/pangolin/lineages_merged.csv",sep=",",header=T)
gisaid = gisaid %>% left_join(new_pango_gisaid,by=c("strain"="taxon"))
gisaid$date_fix= ymd(gisaid$date)
gisaid$pangolin_lineage = gisaid$lineage
gisaid$probability = gisaid$probability

qc_lineages_phased$pangolin_lineage = qc_lineages_phased$lineage
qc_lineages_phased$genbank_accession = "?"
qc_lineages_phased$paper_url = "?"
#qc_lineages_phased$strain = qc_lineages_phased$sample_name_fasta
qc_lineages_phased$location = "Los Angeles County"
qc_lineages_phased$strain = str_replace_all(qc_lineages_phased$strain,":","_")
colnames(gisaid)[!colnames(gisaid) %in% colnames(qc_lineages_phased)]

gisaid_phased = qc_lineages_phased %>% select(colnames(qc_lineages_phased)[colnames(qc_lineages_phased) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
gisaid_phased$date = gisaid_phased$date_fix
gisaid_phased$strain = str_replace_all(gisaid_phased$strain,":","_")
write.table(gisaid_phased,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased.tsv",col.names=T,row.names=F,quote=F,sep="\t")
aa = c(keep_genomes,xstringset)
writeXStringSet(aa,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased.fasta")
##### 
pangolin_new_lineages =read.csv("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/phased/pangolin/lineages.csv")
qc_lineages_phased2 = qc_lineages_phased %>% left_join(pangolin_new_lineages,by=c("strain"="taxon")) %>% mutate(pangolin_lineage= ifelse(is.na(lineage.y),pangolin_lineage,lineage.y),probability= ifelse(is.na(lineage.y),probability.x,probability.y))
qc_lineages_phased2 %>%rename(pangolin_lineage,pcox)
gisaid_phased = qc_lineages_phased2 %>% select(colnames(qc_lineages_phased2)[colnames(qc_lineages_phased2) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
gisaid_phased$date = gisaid_phased$date_fix
gisaid_phased$strain = str_replace_all(gisaid_phased$strain,":","_")
# Current run
# TODO: get the list 
#current_run = read.table("~/Dropbox/COVID19/ncov_new/again/ncov/results/north-america_usa_los-angeles_update_hets/current_sub.txt",header=F)

#gisaid_phased = gisaid_phased %>% filter(strain %in% current_run$V1)

write.table(gisaid_phased,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased_curate.tsv",col.names=T,row.names=F,quote=F,sep="\t")
linnaa = c(keep_genomes_2,xstringset)
#aa = aa[gisaid_phased$strain]
writeXStringSet(aa,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/phased_curate.fasta")

# switch lineages ##

subset_genomes = final_genomes[qc_lineages$sample_name_fasta]


qc_lineages$pangolin_lineage = qc_lineages$lineage
qc_lineages$genbank_accession = "?"
qc_lineages$paper_url = "?"
qc_lineages$strain = qc_lineages$sample_name_fasta
qc_lineages$location = "Los Angeles County"
colnames(gisaid)[!colnames(gisaid) %in% colnames(qc_lineages)]
gisaid_local = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid)]) %>% bind_rows(gisaid)
aa = c(subset_genomes,xstringset)
gisaid_local$date = gisaid_local$date_fix
write.table(gisaid_local,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final.tsv",col.names=T,row.names=F,quote=F,sep="\t")
writeXStringSet(aa,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final.fasta")
### Only LA with Wuhan and LA ####
#calif_wuhan = gisaid %>% filter(location == "Los Angeles County" | grepl("Hu-1",strain))
out_seq = c(subset_genomes,xstringset[calif_wuhan$strain])
gisaid_la = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid)]) 
write.table(gisaid_la ,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la.tsv",col.names=T,row.names=F,quote=F,sep="\t")
writeXStringSet(out_seq,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la.fasta")


writeXStringSet(xstringset[aaa$strain], file="~/Dropbox/COVID19/analysis_for_epi_paper/data/la_icahn.fasta")
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
write.table(gisaid_la_hets ,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la._het.tsv",col.names=T,row.names=F,quote=F,sep="\t")
writeXStringSet(out_seq,file="~/Dropbox/COVID19/ncov_new/again/ncov/data/merged/final_only_la_hets.fasta")


#### Het phylogenies ###
