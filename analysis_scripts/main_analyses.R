### Script will run all the main analysis for the interactive analysis used for the LA epi paper  ## 
library(tidyverse)
library(ggrepel)
library(vcfR)
library(lubridate)
source("R/utils.R")
source("R/load_datasets.R")

#colnames(qc_input2)

## Create plots of date distribution of samples and the sequecing method used ##
#p1 = qc_input2 %>% ggplot(aes(x=date_fix,fill=lib_type)) + geom_histogram(binwidth = 7) + theme_bw() +  scale_x_date(date_labels  = "%b-%d",date_breaks = "4 weeks") + theme(text=element_text(size=24)) + ylab("Count") + 
#  xlab("Date collected") + scale_fill_brewer(palette = "Set1")
## xlab("Date collected") + scale_fill_brewer(palette = "Set2")
##p3=  qc_input2 %>% ggplot(aes(y=(ct),x=pass)) + geom_boxplot() +  scale_color_brewer(palette = "Set1")  + xlab("Pass QC") + ylab("Cycling threshold (Ct)") + geom_jitter() + facet_wrap(~lib_type)  + theme_bw() +theme(text=element_text(size=24))
#plot_grid(p1,p2,p3,ncol=3,width=)



#table(qc_input2$library_type)
#table(qc_input2$library_type, qc_input2$pass)
#table(qc_input2$library_type, qc_input2$sars2_percent_dup > .25)
#summary(lm(qc_input2$sars2_percent_dup ~ qc_input2$ct + qc_input2$library_type + log(qc_input2$total_read_count)))
#drop1(lm(qc_input2$sars2_percent_dup ~ qc_input2$ct + qc_input2$library_type),test="F")

#qc_input2 %>% group_by(library_type) %>% summarise(mean=mean(sars2_percent_dup),median=median(sars2_percent_dup))
#c_input2$coverage_mean_dedup = qc_input2$mapped_sars2_dedup/29903
#p1 = qc_input2 %>% ggplot(aes(y=coverage_3x,xi=ct)) + geom_point() + scale_color_brewer(palette = "Set1") +  facet_wrap(~library_type)  +
#  xlab("Cycling threshold (Ct)") + ylab("Coverage 3x") + theme_bw()
#p2 = qc_input2 %>% ggplot(aes(y=log2(coverage_3x),x=ct)) + geom_point() +  scale_color_brewer(palette = "Set1") +  facet_wrap(~library_type) + xlab("Cycling threshold (Ct)") + ylab("Log2 coverage") + theme_bw()
#p3 = qc_input2 %>% ggplot(aes(y=log2(mapped_sars2_dedup/29903),x=ct)) + geom_point() +  scale_color_brewer(palette = "Set1") +  facet_wrap(~library_type) + xlab("Cycling threshold (Ct)") + ylab("Log2 coverage") + theme_bw()

#p4 = qc_input2 %>% ggplot(aes(y=sars2_percent_dup,x=ct)) + geom_point() +  scale_color_brewer(palette = "Set1")  + facet_wrap(~library_type) +  xlab("Cycling threshold (Ct)") + ylab("SARS-CoV-2 read percentage") + theme_bw()
#p5= qc_input2 %>% ggplot(aes(y=sars2_percent_dedup,x=ct)) + geom_point() +  scale_color_brewer(palette = "Set1")  + facet_wrap(~library_type) +  xlab("Cycling threshold (Ct)") + ylab("SARS-CoV-2 read percentage") + theme_bw()
#p6 = qc_input2 %>% ggplot(aes(y=(ct),x=pass)) + geom_boxplot() +  scale_color_brewer(palette = "Set1")  + xlab("Pass QC") + ylab("Cycling threshold (Ct)") + geom_jitter() + facet_wrap(~library_type)  + theme_bw() +  coord_flip()
#cowplot::plot_grid(p1,p2,p3,p4,p5,p6)#p(~library_type)
#(180e6/exp((( 5.48 + seq(15,40,by=5) * 0.37))))

library(viridisLite)
#qc_input2 %>% mutate(for30x = 30/coverage_mean_dedup * qc_input2$total_read_count) %>% filter(!is.infinite(for30x)) %>% ggplot(aes(y=log(for30x),x=log(total_read_count),color=ct)) + geom_point() + facet_wrap(~library_type) + geom_abline() + scale_color_viridis_c()
#qc_input2 %>% mutate(for30x = 30/coverage_mean_dedup * qc_input2$total_read_count) %>% filter(!is.infinite(for30x))  %>% do(broom::tidy(lm(log(for30x) ~ ct + library_type,.)))
#qc_input2 = qc_input2 %>%filter(coverage_mean_dedup < 30)
#qc_input2 %>% filter(ct < 20) %>% ggplot(aes(x=log2(total_read_count),y=coverage_3x)) + geom_point()

df_neb_amp = data.frame()
#(lm(log(qc_input2$for30x) ~ qc_input2$ct))

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
#### Write out the samples tha pass both analyses.###
write.table(df_neb_amp[df_neb_amp$pass,],file="data/qc/samples_that_worked_in_both_neb_and_amp_july10.csv",quote=F,row.names=F,col.names=T)
merged_vcf_filt= read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july12//uid_type/outputs/quasi_species/merged/all.vcf")
#merged_vcf_filt2= read.vcfR("~/HOFFMAN/rerun_all_new_july12//uid_type/outputs/bcftools_vcfs/merged/filt/sars2/all.vcf.gz")
#png("figur")
p1 = df_neb_amp  %>%  ggplot(aes(y=log(sars2_percent_neb),log(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2")
cor(log(df_neb_amp$sars2_percent_neb),log(df_neb_amp$sars2_percent_amp))
t.test(log(df_neb_amp$sars2_percent_neb),log(df_neb_amp$sars2_percent_amp))
t.test((df_neb_amp$sars2_percent_neb),(df_neb_amp$sars2_percent_amp))

#df_neb_amp  %>%  ggplot(aes(y=log(sars2_percent_neb),log(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2") + facet_wrap(~bin,scales="free")
df_neb_amp  %>%  ggplot(aes(y=(sars2_percent_neb),(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2") + facet_wrap(~bin,scales="free",ncol=5)
ggsave("figures/NEBVRTSEQpercent",width = 16,height = 8,units = "in",dpi = 150)
p2 = df_neb_amp %>%  ggplot(aes(y=(sars2_percent_neb),(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2")
cowplot::plot_grid(p1,p2)
ggsave("figures/NEBVRTSEQ.png",width = 16,height = 8,units = "in",dpi = 150)
p3 = df_neb_amp  %>%filter(coverage_3x_neb > .8 & coverage_3x_amp > .8)  %>% ggplot(aes(y=log(sars2_percent_neb),log(sars2_percent_amp))) + geom_point() + geom_abline() + theme_bw() +ylab("NEB") + xlab("RT-seq") + theme(text=element_text(size=24)) + ggtitle("% SARS-CoV-2")

#summary(lm(df_neb_amp$sars2_percent_neb ~ df_neb_amp$sars2_percent_amp))
#drop1(lm(log(df_neb_amp$sars2_percent_neb) ~ log(df_neb_amp$sars2_percent_amp)),test="F")
#drop1(lm((df_neb_amp$sars2_percent_neb) ~ (df_neb_amp$sars2_percent_amp)),test="F")
#summary(lm(log(df_neb_amp$sars2_percent_neb) ~ log(df_neb_amp$sars2_percent_amp)))

#table(qc_input2$pass,qc_input2$library_type)
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


genomes_in  = qc_lineages_phased2[!is.na(qc_lineages_phased2$lineage.y),]$merged_id

gt_comparison


gt_comparison %>% filter(!equal_var)
gt_comparison_new = gt_comparison %>% filter(neb == "0/1"|amp=="0/1") %>%  left_join(qc_filt_double_ups,by=c("UID"="uid"))
#gt_comparison_new = gt_comparison %>% filter(!equal_var) %>% left_join(qc_filt_double_ups,by=c("UID"="uid"))
## Look at discordant variants
write.table(gt_comparison_new, file="data/qc/qc_variants_all_hets.tsv",quote=F,row.names=F,col.names=T,sep="\t")
gt_comparison %>% filter(neb == "0/1"|amp=="0/1") %>% mutate(manual_check=ifelse(!equal_var,"YES","NA"))  %>% ggplot(aes(y=ar_amp,ar_neb,color=manual_check,label=snp_names)) + geom_point() +
  theme_bw() + facet_wrap(~UID) + geom_abline() + scale_color_brewer(palette = "Set1") + xlab("NEB") + ylab("Amplicon") + ylim(c(0,1)) + xlim(c(0,1)) + geom_label_repel()### ALL the het sites look great.####
### TODO: remove het sites from the consensus fasta files ###
### So if we check the pangolin lineages as expected all of them are exactly the same ####


qc_lineages %>% inner_join(qc_lineages_phased2 %>% filter(!duplicated(merged_id)) %>% filter(!is.na(lineage.y)) , by=c("sample_name_fasta"="sample_name_fasta"))

qc_lineages_phased2 %>%  filter(!duplicated(merged_id)) %>% filter(!is.na(lineage.y)) 
get_concordance_list_ar(qc_lineages_phased2 %>%  filter(!duplicated(merged_id)) %>% filter(!is.na(lineage.y)) ,merged_vcf_filt)
#qc_input2 = qc_input2 %>% mutate(qc_lineage=ifelse(is.na(status),"failed","passed_qc"))
#qc_input2 %>% pivot_wider(uid,names_from=library_type,names_sep =".",values_from=c(lineage,qc_lineage,UFbootstrap)) %>% filter(!is.na(lineage.amp) &  !is.na(lineage.neb))
#### 100% concordance in SNPs so far. 
#qc_input2 %>% ggplot(aes(x=lineage)) + geom_bar()
#qc_lineages = qc_input2 %>% filter(pass) %>% filter(!duplicated(uid)) 
#library(lubridate)
### 
#qc_lineages$date_fix = mdy(qc_lineages$date)
#qc_lineages$week =  week(qc_lineages$date_fix)
#qc_lineages$month = month(qc_lineages$date_fix)
#qc_lineages %>%ggplot(aes(x=lineage)) + geom_bar() + facet_wrap(~month)
#qc_lineages %>% filter(nextstrain_clade != " -- ")  %>% ggplot(aes(x=nextstrain_clade)) + geom_bar() + facet_wrap(~month)

#qc_lineages %>% ggplot(aes(x=week,fill=nextstrain_clade)) + geom_bar(position = "stack",width=.5) + theme_bw() + scale_fill_brewer(palette = "Set3")
#qc_lineages %>% ggplot(aes(x=week,fill=lineage)) + geom_bar(position = "stack",width=.5) + theme_bw() +  scale_fill_brewer(palette = "Set3")
#qc_lineages %>% filter(nextstrain_clade != " -- ")  %>% ggplot(aes(x=gisaid_clade)) + geom_bar() + facet_wrap(~month)
#summary(lm(qc_lineages$week ~ qc_lineages$nextstrain_clade))

#### We can look at the annovar results ####
#### ####

#### ####


#### 

#not_unique_hets %>% filter(SAMPLE %in% c("LA-0193","LA-0061","LA-0172","LA-0197","LA-0198"))#
#gt_comparison %>% filter()
#aaa = get_gt_rbind_df(qc_lineages,merged_vcf_filt) 

a#aa = aaa %>% filter(sample_id %in%  c("LA-0193","LA-0061","LA-0172","LA-0197","LA-0198"))

#aaa %>% pivot_longer(cols="gt",names_from = sample_id)
#
#gt = aaa %>% pivot_wider(names_from = sample_id,values_from=gt)
#gt[apply(gt[,2:6],1,function(x){ length(unique(x))}) > 1,]

#vcf_phased = read.vcfR("~/HOFFMAN/test.vcf.gz")
                                    
#(vcf_phased@gt[,colnames(vcf_phased@gt) == "LA-0061"])

#merged_vcf_filt2= read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july12//uid_type/outputs/quasi_species/merged/all_from_merged.vcf")
#merged_vcf_filt = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july12//uid_type/outputs/quasi_species/merged/all_from_merged.vcf")

#merged_vcf_filt@gt[,colnames(merged_vcf_filt@gt) %in% qc_lineages$sample_name_fasta]
#cbind(merged_vcf_filt2@gt[,colnames(merged_vcf_filt2@gt) == "LA-0041"] ,merged_vcf_filt@gt[,colnames(merged_vcf_filt@gt) == "LA-0041"] )#json_in =fromJSON(file="~/Dropbox/COVID19/ncov_new/again/ncov/results/north-america_usa_los-angeles_update/clades.json")

