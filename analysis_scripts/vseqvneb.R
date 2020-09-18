### Script will run all the main analysis for the interactive analysis used for the LA epi paper  ## 
library(tidyverse)
library(ggrepel)
library(vcfR)
library(magrittr)
library(lubridate)
library(cowplot)
source("R/utils.R")
source("R/load_datasets.R")
## Create plots of date distribution of samples and the sequecing method used ##
good_uids = gisaid_phased_filt%>% filter(submitting_lab == "UCLA") %>% inner_join(qc_input2,by=c("strain"="sample_name_fasta"))
removed = qc_lineages %>% filter(!(sample_name_fasta %in% gisaid_phased_filt$strain)) %>% select(sample_name_fasta)
qc_input2 = qc_input2 %>% mutate(pass_het_filter = ifelse(sample_name_fasta %in% removed$sample_name_fasta,FALSE,TRUE)) %>% mutate(pass = pass & pass_het_filter)
qc_input2 = qc_input2 %>% mutate(duplication_rate =  1 - qc_input2$mapped_sars2_dedup /qc_input2$mapped_sars2_dup)
qc_input2 %>% select(sample_name_fasta,uid,ct,date_fix,total_read_count,mapped_sars2_dedup,sars2_percent_dedup,mapped_sars2_dup,mapped_sars2_dup,sars2_percent_dup,duplication_rate,coverage_3x,coverage_5x,pass) %>%
  dplyr::rename(sample_name=sample_name_fasta,date_collected=date_fix, read_count=total_read_count,coverage_3x_dedu=coverage_3x,coverage_5x_dup=coverage_5x) %>% write.table(file = "paper/tables//S1.tsv",col.names=T, row.names=F,quote = F,sep="\t")
### ####
g3 = (qc_input2 %>% left_join(dd,by=c("merged_id"="merged_id")) %>% group_by(group) %>% mutate(n2=sum(!is.na(group)))) %>% filter(n2 > 1) %>% filter(row_number() > 1) %>% ungroup
#g2 = (qc_input2 %>% left_join(dd,by=c("merged_id"="merged_id")) %>% group_by(group) %>% mutate(n2=sum(!is.na(group)))) %>% filter(n2== 0) %>% ungroup

gg= bind_rows(g1,g2)
qc_input2 = qc_input2 %>%  mutate(pass_het_filter = ifelse(sample_name_fasta %in% removed$sample_name_fasta | merged_id %in% g3$merged_id,FALSE,TRUE)) %>% mutate(pass = pass & pass_het_filter) 
gisaid_phased_filt %>% filter(submitting_lab == "UCLA")%>% summarise(min(date_fix),max=max(date_fix))


qc_input2 %>%group_by(lib_type) %>% dplyr::summarise(coverage_3x=mean(coverage_3x),total_reads=median(total_read_count),n=n(), passing=sum(pass))

qctable(qc_input2$library_type)
table(qc_input2$library_type, qc_input2$pass)
table(qc_input2$library_type, qc_input2$sars2_percent_dup > .25)
summary(lm(qc_input2$sars2_percent_dup ~ qc_input2$ct + qc_input2$library_type + log(qc_input2$total_read_count)))
drop1(lm(qc_input2$sars2_percent_dup ~ qc_input2$ct + qc_input2$library_type),test="F")
p1 = qc_input2 %>% ggplot(aes(y=coverage_3x,x=ct,color=lib_type)) + geom_point() + scale_color_brewer(palette = "Set1",name="Library type")  +
  xlab("Cycling threshold (Ct)") + ylab("Proportion of SARS-CoV-2 genome at >= 3x") + theme_bw() + theme(legend.position = "none") + geom_hline(yintercept = 0.8)
#p2 = qc_input2 %>% ggplot(aes(y=log2(coverage_3x),x=ct)) + geom_point() +  scale_color_brewer(palette = "Set1") +  facet_wrap(~library_type) + xlab("Cycling threshold (Ct)") + ylab("Log2 coverage") + theme_bw()
#p3 = qc_input2 %>% ggplot(aes(y=log2(mapped_sars2_dedup/29903),x=ct)) + geom_point() +  scale_color_brewer(palette = "Set1") +  facet_wrap(~library_type) + xlab("Cycling threshold (Ct)") + ylab("Log2 coverage") + theme_bw()
p2 = qc_input2 %>% ggplot(aes(y=sars2_percent_dup,x=ct,color=lib_type)) + geom_point() +  scale_color_brewer(palette = "Set1",name="Library type")  +  xlab("Cycling threshold (Ct)") + ylab("Proportion of reads mapping to SARS-CoV-2") + theme_bw()
#p3  = qc_input2 %>% ggplot(aes(y=log2(map),x=ct,color=lib_type)) + geom_point() +  scale_color_brewer(palette = "Set1",name="Library type")  +  xlab("Cycling threshold (Ct)") + ylab("Proportion of reads mapping to SARS-CoV-2") + theme_bw()
p3 = qc_input2 %>% ggplot(aes(y=duplication_rate,x=lib_type,color=lib_type)) + geom_boxplot() + geom_jitter() +  scale_color_brewer(palette = "Set1",name="Library type") + theme_bw() + ylab("Duplication rate") + xlab("Library Type")
#p4 = qc_input2 %>% ggplot(aes(y=log2(mapped_sars2_dedup),x=log2(mapped_sars2_dup),color=lib_type)) + geom_point()+  scale_color_brewer(palette = "Set1",name="Library type") + theme_bw() + ylab("Log2 SARS2 mapped (duplicates present)") + xlab("Log 2 SARS2 mapped (no duplicates present)") + geom_abline()
p4 = qc_input2 %>% ggplot(aes(y=sars2_percent_dedup,x=ct,color=lib_type)) + geom_point() +  scale_color_brewer(palette = "Set1",name="Library type")  +  xlab("Cycling threshold (Ct)") + ylab("Proportion of reads mapping to SARS-CoV-2") + theme_bw()
leg = get_legend(p2 + theme(text=element_text(size=20)))
pf = cowplot::plot_grid(p1 + theme(text=element_text(size=20)),p2+ theme(text=element_text(size=20)) +  theme(legend.position = "none") ,p3+theme(text=element_text(size=20)) +  theme(legend.position = "none"), p4+theme(text=element_text(size=20)) + theme(legend.position = "none") ,labels = "AUTO",rel_widths = c(5,5),label_size = 24)
plot_grid(pf,leg, rel_widths = c(1,.15))

ggsave(file="paper/figures/S1.png",dpi=300,width=16,height = 16)
p1 = qc_input2 %>% ggplot(aes(x=date_fix,fill=lib_type)) + geom_histogram(binwidth = 7) + theme_bw() +  scale_x_date(date_labels  = "%b-%d",date_breaks = "4 weeks") + theme(text=element_text(size=16)) + ylab("Count") + 
  xlab("Date collected") + scale_fill_brewer(palette = "Set1",name="Library type")
p2 = qc_input2 %>% ggplot(aes(x=date_fix,fill=pass)) + geom_histogram(binwidth = 7)+ theme_bw() +  scale_x_date(date_labels  = "%b-%d",date_breaks = "4 weeks") + theme(text=element_text(size=16)) + ylab("Count") +
  xlab("Date collected") + scale_fill_brewer(palette = "Set2", name="Passing QC")
p3 = qc_input2 %>% ggplot(aes(y=(ct),x=pass)) + geom_boxplot(outlier.shape = NA) +  scale_color_brewer(palette = "Set1")  + xlab("Pass QC") + ylab("Cycling threshold (Ct)") + geom_jitter() + facet_wrap(~lib_type)  + theme_bw() + theme(text=element_text(size=16)) 
p4 = qc_input2 %>% ggplot(aes(y=coverage_3x,x=ct)) + geom_point() + scale_color_brewer(palette = "Set1") +  facet_wrap(~lib_type)  +
  xlab("Cycling threshold (Ct)") + ylab("Coverage 3x") + theme_bw()
plot_grid(p1,p2,p3,ncol=3,align="h",labels = "AUTO",label_size = 24)
ggsave(file="paper/figures/S2.png",dpi=300,width=16,height = 8)
library(viridisLite)
df_neb_amp = data.frame()
split_rows = split(qc_input2,f=qc_input2$uid)
for (row in split_rows){
  if(nrow(row) ==2){
    neb = row[row$library_type == "neb",]
    amp = row[row$library_type == "amp",]
    df_tmp = data.frame(merged_id=neb$merged_id, sars2_percent_neb=neb$sars2_percent_dup,coverage_3x_neb=neb$coverage_3x,coverage_mean_neb_dedup=neb$coverage_mean_dedup,
                        coverage_mean_neb_dup=neb$coverage_mean_no_dedup,
                        sars2_percent_amp=amp$sars2_percent_dup,coverage_3x_amp=amp$coverage_3x,
                        coverage_mean_amp_dedup=amp$coverage_mean_dedup,
                        coverage_mean_amp_dup=amp$coverage_mean_no_dedup
                        ,ct=neb$ct,merged_id_amp=amp$merged_id,pass=neb$pass & amp$pass,uid=neb$uid,bin=neb$bin_ct,total_neb=neb$total_read_count,total_amp=amp$total_read_count)
    df_neb_amp = rbind(df_neb_amp,df_tmp)
  }
  
}
write.table(df_neb_amp[df_neb_amp$pass,],file="data/qc/samples_that_worked_in_both_neb_and_amp_july10.csv",quote=F,row.names=F,col.names=T)
merged_vcf_filt= read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20//uid_type/outputs/quasi_species/merged/all.vcf")
p1 = df_neb_amp  %>%  ggplot(aes(y=log(sars2_percent_neb),log(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("V-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2")
cor(log(df_neb_amp$sars2_percent_neb),log(df_neb_amp$sars2_percent_amp))
t.test(log(df_neb_amp$sars2_percent_neb),log(df_neb_amp$sars2_percent_amp))
t.test((df_neb_amp$sars2_percent_neb),(df_neb_amp$sars2_percent_amp))
df_neb_amp  %>%  ggplot(aes(y=(sars2_percent_neb),(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2") + facet_wrap(~bin,scales="free",ncol=5)
ggsave("figures/NEBVRTSEQpercent",width = 16,height = 8,units = "in",dpi = 150)
p2 = df_neb_amp %>%  ggplot(aes(y=(sars2_percent_neb),(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2")
cowplot::plot_grid(p1,p2)
ggsave("figures/NEBVRTSEQ.png",width = 16,height = 8,units = "in",dpi = 150)
p3 = df_neb_amp  %>%filter(coverage_3x_neb > .8 & coverage_3x_amp > .8)  %>% ggplot(aes(y=log(sars2_percent_neb),log(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2")
df_filt = df_neb_amp %>% filter(coverage_3x_neb > 0.8 & coverage_3x_amp > 0.8)
qc_filt_double_ups = qc_input2[qc_input2$uid %in% df_filt$uid,]

ar_list = get_concordance_list_ar(qc_filt_double_ups,merged_vcf_filt)
genotype_comparison = bind_rows(ar_list$genotype_lists,.id="UID")
gt_comparison = genotype_comparison %>% filter(amp != "./." & neb != "./.") %>% mutate(equal_var = neb == amp)
## Including hets.
gt_comparison %>% summarise(equal_var=sum(equal_var)/length(equal_var))
# 0.9993648
# 0.9993648
gt_comparison %>% filter(neb != "0/1" & amp != "0/1" )%>% summarise(equal_var=sum(equal_var)/length(equal_var),n=n())
gt_comparison %>% filter(neb != "0/1" & amp != "0/1" )%>%  summarise(equal_var=sum(equal_var)/length(equal_var),n=n())
gt_comparison %>% filter(neb != "0/1" & amp != "0/1" )%>% summarise(equal_var=sum(equal_var)/length(equal_var),n=n())
gt_comparison %>% filter(neb == "0/1" | amp == "0/1" )
gt_comparison %>% filter(neb == "0/1" | amp == "0/1" ) %>% filter(!equal_var)
gt_comparison %>% filter(neb == "0/1" | amp == "0/1" ) %>% filter(equal_var)
### All validated in IGV!
df_filt = df_neb_amp %>% filter(coverage_5x_neb > 0.8 & coverage_5x_amp > 0.8)
qc_filt_double_ups = qc_input[qc_input$uid %in% df_filt$uid,]
concordance_list = get_concordance_list(qc_filt_double_ups,merged_vcf_filt)
### SNPS look good between methods ###

gt_comparison %>% filter(!equal_var)
gt_comparison_new = gt_comparison %>% filter(neb == "0/1"|amp=="0/1") %>%  left_join(qc_filt_double_ups,by=c("UID"="uid"))
#gt_comparison_new = gt_comparison %>% filter(!equal_var) %>% left_join(qc_filt_double_ups,by=c("UID"="uid"))
## Look at discordant variants
write.table(gt_comparison_new, file="data/qc/qc_variants_all_hets.tsv",quote=F,row.names=F,col.names=T,sep="\t")
gt_comparison %>% filter(neb == "0/1"|amp=="0/1") %>% mutate(manual_check=ifelse(!equal_var,"YES","NA"))  %>% ggplot(aes(y=ar_amp,ar_neb,color=manual_check,label=snp_names)) + geom_point() +
  theme_bw() + facet_wrap(~UID) + geom_abline() + scale_color_brewer(palette = "Set1") + xlab("NEB") + ylab("Amplicon") + ylim(c(0,1)) + xlim(c(0,1)) + geom_label_repel()### ALL the het sites look great.####
### TODO: remove het sites from the consensus fasta files ###

#### Histogram of all collected samples #### 

qc_input2 %>% filter(pass) %>% ggplot(aes(x=date_fix)) + geom_histogram(binwidth = 7) + theme_bw() +  scale_x_date(date_labels  = "%b-%d",date_breaks = "4 weeks") + theme(text=element_text(size=16)) + ylab("Count") + 
  xlab("Date collected") + scale_fill_brewer(palette = "Set1",name="Library type")
ggsave(file="paper/figures/S4.png",dpi=300,width=16,height = 8)

