source("R/load_datasets.R")
qc_lineages %>%ggplot(aes(x=date_fix)) + geom_histogram(binwidth = 7) + theme_bw() + theme(text=element_text(size=24)) + scale_x_date(date_labels  = "%b-%d",date_breaks = "2 weeks") + ylab("Count") + xlab("Date collected") 
ggsave("figures/1A.png",width=12,height=6)
library(viridisLite)

### Only places reportinng > 50 genomes, before week 17.
gisaid_longhua_guo = gisaid_local  %>% filter(week < 17) %>% filter(country == "USA")  
write.table(gisaid_longhua_guo,file="data/usa_gisaid.txt",quote=F,row.names=F,col.names=T,sep="\t")
write.table(gisaid_local,file="data/gisaid_all.txt",quote=F,row.names=F,col.names=T, sep="\t")

gisaid_local$pango_new =  ifelse(gisaid_local$pangolin_lineage %in% names(color_map),gisaid_local$pangolin_lineage,"other")
gisaid_local %>% group_by(region,pango_new) %>% summarise(n = n()) %>%
     mutate(freq = n / sum(n))

gisaid_local$month = month(gisaid_local$date_fix)
df_b143 =gisaid_local %>% filter(pangolin_lineage == "B.1.43") %>% group_by(location) %>% mutate(n=n()) %>% filter(n > 10) %>% filter(location != "")

gisaid_local %>% group_by(location, week) %>% summarise(b143=sum(pangolin_lineage == "B.1.43")) %>% filter(b143 >= 0.0) %>% filter(location %in% df_b143$location) %>% ggplot(aes(y=b143,x=week, color=location)) + geom_point()

gisaid_local %>% mutate(USA = ifelse(country == "USA","USA","Other"))  %>% group_by(USA, week)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=week, color=USA)) + geom_point() + geom_line()

gisaid_local %>% filter(division %in% c("California","Washington","Utah","Arizona","Oregon"))  %>% group_by(division, week)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=week, color=division))  + geom_line()

gisaid_local %>% filter(division %in% c("California","Washington","Utah","Arizona","Oregon"))  %>% group_by(division, month)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")/length(pangolin_lineage)) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=month, color=division))  + geom_line()

gisaid_local %>% filter(division %in% c("California","Washington")) %>% filter(location %in% df_b143$location) %>% group_by(location, month)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")/length(pangolin_lineage)) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=month, color=location))  + geom_line()


gisaid_local %>% filter(division %in% c("California","Washington")) %>% filter(location %in% df_b143$location) %>% group_by(location, month)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=month, fill=location))  + geom_bar(stat="identity")

gisaid_local %>% filter(division %in% c("California","Washington")) %>% filter(location %in% df_b143$location) %>% group_by(location, month)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")/length(pangolin_lineage)) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=month, fill=location))  + 
  geom_bar(stat="identity", position = position_dodge(preserve = 'single')) + scale_fill_brewer(palette = "Set1")

gisaid_local %>% filter(division %in% c("California","Washington")) %>% filter(location %in% df_b143$location) %>% group_by(location, month)%>% summarise(b143=sum(pangolin_lineage == "B.1.43")) %>% filter(b143 >= 0.0) %>% ggplot(aes(y=b143,x=month, fill=location))  + 
  geom_bar(stat="identity", position = position_dodge(preserve = 'single')) + scale_fill_brewer(palette = "Set1")



p1 = all_la %>%filter(week < 17)  %>% ggplot(aes(x=date_fix,fill=pango_new)) + geom_histogram(binwidth =7) +scale_map  +theme_bw() + theme(text=element_text(size=24)) + ggtitle("LA County") + xlab("Date")

aa = gisaid_local %>% filter(location!="")  %>% filter(week < 17) %>% filter(country == "USA")  %>% group_by(location) %>% summarise(n=n()) %>% filter(n > 50) 
gisaid_local  %>% filter(week < 17) %>% filter(country == "USA") %>% filter(location %in% aa$location) %>%  group_by(location,pango_new) %>% summarise(n = n()) %>%
mutate(freq = n / sum(n))%>% ggplot(aes(y=freq,x=location,fill=pango_new)) + geom_bar(stat="identity") + scale_fill_manual(values = color_map) + theme(axis.text.x=element_text(angle=90, hjust=1))


need_colors = names(which(table(all_la$pangolin_lineage) >= 5))
color_map = viridisLite::viridis(length(need_colors))
color_map = c(color_map,"grey50")
names(color_map) = c(need_colors,"other")
scale_map = scale_fill_manual(name="Lineage",values=c(color_map))
color_map = RColorBrewer::brewer.pal(length(need_colors) + 1  , name ="Set3")
names(color_map) = c(need_colors,"other")
scale_map2 = scale_fill_manual(name="Lineage",values = color_map)

gisaid_strain$pango_new = factor(gisaid_strain$pango_new)
#gisaid_strain %>% filter(location=="King County") %>% mutate(pango_new2 = color_map[match(pango_new,names(color_map))])  
p1 = all_la %>%filter(week < 17)  %>% ggplot(aes(x=date_fix,fill=pango_new)) + geom_histogram(binwidth= 10) +scale_map2 +theme_bw() + theme(text=element_text(size=24))  + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("LA County") + xlab("Date")


## Do something to compare the dates??? so it isn't clipping the last week. ##
p1 = all_la %>%filter(week < 17)  %>% ggplot(aes(x=date_fix,fill=pango_new)) + geom_histogram(binwidth =7) +scale_map2  +theme_bw() + theme(text=element_text(size=24))  + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("LA County") + xlab("Date")
p2 = gisaid_strain %>% filter(location=="King County")  %>% filter(week < 17) %>% ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("King County") + xlab("Date")
p3 = gisaid_strain %>% filter(location=="Santa Clara County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("Santa Clara County") + xlab("Date")
p4 =  gisaid_strain %>% filter(division=="New York") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("New York State") + xlab("Date")
p5 = gisaid_strain %>% filter(location=="Grand Princess Cruise Ship") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("Grand Princess") + xlab("Date")
p6 = gisaid_strain %>% filter(location=="Solano County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("Solano County") + xlab("Date")
p7 = gisaid_strain %>% filter(division == "New Jersey") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("New Jersey") + xlab("Date")


p8= gisaid_strain %>% filter(division == "Connecticut") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("Connecticut") + xlab("Date")
p9 =  gisaid_strain %>% filter(division == "Louisiana") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("Louisiana") + xlab("Date")
p10 =  gisaid_strain %>% filter(location == "San Diego") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("San Diego") + xlab("Date")
p11 = gisaid_strain %>% filter(location == "San Francisco") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("San Francisco") + xlab("Date")
p12 = gisaid_strain %>% filter(location == "Whatcom County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("Whatcom County") + xlab("Date")
p13 = gisaid_strain %>% filter(location == "San Luis Obispo County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("San Luis Obispo County") + xlab("Date")


pf = cowplot::plot_grid(p2 + theme(legend.position = "none"),p3 +  theme(legend.position = "none"),p4 +  theme(legend.position = "none"),p5 +  theme(legend.position = "none"),p6+  theme(legend.position = "none"),p7+theme(legend.position = "none"),
                        p8 +theme(legend.position = "none"),p9 + theme(legend.position = "none"),p10+theme(legend.position = "none"),p11 + theme(legend.position = "none"),p12 + theme(legend.position = "none"),p13 + theme(legend.position = "none") ,ncol=2,align="hv")
cowplot::plot_grid(p1,pf,ncol=2)
ggsave("figures/S4.png",width=20,height=14,dpi=150)
#p5 =  gisaid_strain %>% filter(location=="Sonoma County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("New York State") + xlab("Date")



#p6 =  gisaid_strain %>% filter(location=="San Mateo County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map2 +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("New York State") + xlab("Date")

#p7 =  gisaid_strain %>% filter(location=="San Benito County") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + ggtitle("New York State") + xlab("Date")


gisaid_strain %>% filter(authors == "Xianding Deng et al") %>%  filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-17"))) + facet_wrap(~location,scales="free_y")

 =  gisaid_strain %>% filter(division=="Louisiana") %>% filter(week < 17) %>%  ggplot(aes(x=date_fix,fill=pango_new)) +  geom_histogram(binwidth =7) +scale_map +theme_bw() + theme(text=element_text(size=24)) + scale_x_date(limits=as.Date(c("2020-02-15","2020-04-26"))) + ggtitle("New York State") + xlab("Date")
main = plot_grid(p1 + theme(legend.position = "none") ,p2 + theme(legend.position = "none") ,p3 + theme(legend.position = "none") ,p4 + theme(legend.position = "none"),nrow=4)
plot_grid(main,get_legend(p1),rel_widths   = c(1,.4))
ggsave("figures/S4.png",width=12,height=12)
all_la$pango_new = factor(all_la$pango_new)
p1 = all_la %>%filter(week >= 17)  %>% ggplot(aes(x=date_fix,fill=pango_new)) + geom_histogram(binwidth =7) +scale_map2  +theme_bw() + theme(text=element_text(size=24))  + ggtitle("LA County") + xlab("Date")
#p1 = all_la %>%filter(week >= 17)  %>% ggplot(aes(x=date_fix,fill=pango_new)) + geom_histogram(binwidth =7) +scale_map  +theme_bw() + theme(text=element_text(size=24))  + ggtitle("LA County") + xlab("Date")
p1
ggsave("figures/S5.png",width=12,height=6)

#p2 = all_la %>% mutate(is_b1=grepl("B.1",pangolin_lineage)) %>% ggplot(aes(x=date_fix,fill=is_b1)) + geom_histogram(binwidth =7) + scale_fill_brewer(palette = "Set3") + theme_bw() + theme(text=element_text(size=24))
#all_la %>% ggplot(aes(x=date_fix,fill=pangolin_lineage)) + geom_histogram(binwidth =7) + scale_fill_viridis_d() +theme_bw() + theme(text=element_text(size=24))

cowplot::plot_grid(p1,p2,nrow=2)
library(Biostrings)
all_la$week = week(all_la$date_fix)
all_la %>% group_by(pangolin_lineage) %>% summarise(count=n())

gisaid_strain %>% mutate(is_a=grepl("^A.1$",pangolin_lineage)) %>% filter(location=="King County") %>%  ggplot(aes(x=date_fix,fill=is_a)) +  geom_histogram(binwidth =7) + scale_fill_viridis_d() +theme_bw() + theme(text=element_text(size=24))
gisaid_strain %>% mutate(is_b=grepl("^B$",pangolin_lineage)) %>% filter(location=="Santa Clara County") %>%  ggplot(aes(x=date_fix,fill=is_b)) +  geom_histogram(binwidth =7) + scale_fill_viridis_d() +theme_bw() + theme(text=element_text(size=24))

all_la %>% mutate(is_a=grepl("^A.1$",pangolin_lineage)) %>%   ggplot(aes(x=date_fix,fill=is_a)) + geom_histogram(binwidth =7) + scale_fill_viridis_d() +theme_bw() + theme(text=element_text(size=24))
gisaid_strain %>% mutate(is_b=grepl("^B$",pangolin_lineage)) %>% filter(location=="Santa Clara County") %>%  ggplot(aes(x=date_fix,fill=is_b)) +  geom_histogram(binwidth =7) + scale_fill_viridis_d() +theme_bw() + theme(text=element_text(size=24))



#gisaid %>% filter(week < 17) %>% filter(division=="Washington") %>% ggplot(aes(x=date_fix,fill=pangolin_lineage)) + geom_histogram(binwidth = 7)
#gisaid %>% filter(week < 17) %>% filter(location=="Santa Clara County") %>% ggplot(aes(x=date_fix,fill=pangolin_lineage)) + geom_histogram(binwidth = 7)

qc_lineages %>% mutate(lin_bp= grepl("B.1",lineage)) %>%  ggplot(aes(x=week,fill=lin_bp)) + geom_bar(position = "stack",width=.5) + theme_bw() +  scale_fill_brewer(palette = "Set1")


#summary(lm(qc_lineages$month ~ qc_lineages$nextstrain_clade))
source("R/utils.R")
merged_vcf_filt = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20///uid_type/outputs/final/merged/quasi_species/all.vcf")

gt_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
gt_matrix_raw = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)

### Remove all heterozygote sites ###
gt_matrix[gt_matrix == "0/1" | gt_matrix == "0/2"] = "./."
# Convert to numeric ###
#rownames(gt_matrix)
qc_lineages[qc_lineages$month >=5,]
gt_matrix_may= gt_matrix[,which(colnames(gt_matrix) %in% qc_lineages$sample_name_fasta)]
missing = apply(gt_matrix_may,1, function(x){sum(x == "./.")/length(x)})
freq =  apply(gt_matrix_may,1, function(x){(2*sum(x != "0/0"))/(2*length(x[x!="./."]))})
snp_df_summary = data.frame(snp_names=rownames(gt_matrix_may),missing=missing,freq=freq,pos=as.numeric(merged_vcf_filt@fix[,2]),ref=merged_vcf_filt@fix[,4],alt=merged_vcf_filt@fix[,5])
snp_df_summary %>% filter(missing < 0.1) %>% mutate(label_snp_new = ifelse(freq > .1,snp_names, "")) %>% ggplot(aes(x=pos,y=freq,label=label_snp_new)) + 
  geom_point() + geom_label_repel() + theme_bw() + ylim(c(0,1))

missing = apply(gt_matrix,1, function(x){sum(x == "./.")/length(x)})
freq =  apply(gt_matrix,1, function(x){(2*sum(x != "0/0"))/(2*length(x[x!="./."]))})
snp_df_summary = data.frame(snp_names=rownames(gt_matrix),missing=missing,freq=freq,pos=as.numeric(merged_vcf_filt@fix[,2]),ref=merged_vcf_filt@fix[,4],alt=merged_vcf_filt@fix[,5])
snp_df_summary %>% filter(missing < 0.1) %>% mutate(label_snp_new = ifelse(freq > .1,snp_names, "")) %>% ggplot(aes(x=pos,y=freq,label=label_snp_new)) + 
  geom_point() + geom_label_repel() + theme_bw() + ylim(c(0,1))
genotypes_df = as.data.frame(t(gt_matrix))
genotypes_df = merge(genotypes_df,qc_lineages,by.x="row.names",by.y="sample_name_fasta")
### Read in annnovar for function ####

p1 = genotypes_df %>%  filter(`23403_A_G` != "./.") %>% group_by(week,`23403_A_G`) %>% summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%ggplot(aes(y=freq, x=week,fill=`23403_A_G`)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Set1")
p2 = genotypes_df %>% filter(`23403_A_G` != "./.") %>%ggplot(aes(x=date_fix,fill=`23403_A_G`)) + geom_histogram(binwidth = 7) + scale_fill_brewer(palette = "Set1") +theme_bw() + theme(text=element_text(size=24))
#co
#genotypes_df %>%  filter(`23403_A_G` != "./.") %>% group_by(month,`23403_A_G`) %>% summarise(n = n()) %>%
#  mutate(freq = n / sum(n)) %>%ggplot(aes(y=freq, x=month,fill=`23403_A_G`)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Set1")

genotypes_df %>% filter(`23403_A_G` != "./.") %>%ggplot(aes(y=ct,x=`23403_A_G`)) + geom_boxplot() + scale_fill_brewer(palette = "Set1") + theme_bw()
spike_filt = genotypes_df %>% filter(`23403_A_G` != "./.") 
summary(lm(spike_filt$ct[spike_filt$library_type == "neb"] ~ spike_filt$`23403_A_G`[spike_filt$library_type == "neb"]))

s#ummary(lm(spike_filt$ct[spike_filt$library_type == "amp"] ~ spike_filt$`23403_A_G`[spike_filt$library_type == "amp"]))
#lm(genotypes_df$)

p1=het_summary %>% left_join(snp_df_summary,by=c("snp_names"="snp_names")) %>%  mutate(label_snp_new = ifelse(freq.x> .01 | freq.y > 0.05,snp_names, "")) %>% ggplot(aes(y=freq.x,x=pos.x,label=label_snp_new)) + geom_point() + geom_label_repel() +ggtitle("Het frequency")
p2 = het_summary %>% left_join(snp_df_summary,by=c("snp_names"="snp_names")) %>%  mutate(label_snp_new = ifelse(freq.x> .01 | freq.y > 0.05,snp_names, "")) %>% mutate(freq_y_name=ifelse(freq.y > 0.5,1-freq.y,freq.y)) %>% ggplot(aes(y=freq_y_name,x=pos.x,label=label_snp_new)) + geom_point() + geom_label_repel() + ggtitle("Allele frequency")
cowplot::plot_grid(p1,p2,nrow=2)

gt_matrix_hets = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
gt_matrix_hets = gt_matrix_hets[(apply(gt_matrix_hets,1,function(x){any(x== "0/1" | x== "0/2")})),]
het_df = as.data.frame(t(gt_matrix_hets))
het_freq =  apply(gt_matrix_hets,1, function(x){(2*sum(x == "0/1" | x== "0/2"))/(2*length(x[x!="./."]))})
het_missing =  apply(gt_matrix_hets,1, function(x){sum(x == "./.")/length(x)})
het_count = apply(gt_matrix_hets,1, function(x){(sum(x == "0/1" | x== "0/2"))})

het_sample_count= apply(gt_matrix_hets,2,function(x){sum(x == "0/1" | x == "0/2")})
atleast_one_het = het_sample_count >= 1
het_summary_sample = data.frame(het_sample_count=het_sample_count,atleast_one=atleast_one_het,sample_id=colnames(gt_matrix_hets)) 
het_summary_sample %>% mutate(merged_label= ifelse(het_sample_count > 1, sample_id,"")) %>% ggplot(aes(x=het_sample_count,label=merged_label)) + geom_histogram() + theme_bw() + xlab("Intra-patient variable count") + ylab("Count")

het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=atleast_one_het)) + geom_bar() + scale_fill_brewer(palette = "Set1")
het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=factor(het_sample_count))) + geom_bar(width=.9) + scale_fill_brewer(palette = "Set1")

unique_hets= lapply(1:ncol(gt_matrix_hets),function(x){ print(x);
  idx_hets = which(gt_matrix_hets[,x] == "0/1" | gt_matrix_hets[,x] == "0/2")
 # print(idx_hets)
  #print(idx_hets)
  if(length(idx_hets) > 0){
    gts=gt_matrix_hets[idx_hets,x]
    #rint(gts)
    xx= sapply(1:length(idx_hets), function(y) { any(gts[y] == gt_matrix_hets[idx_hets[y],-x]) })
    
    
    #p#rint(xx)
    #print(sum(xx,na.rm=T))
    #print()
    print(xx)
    print(gts[xx])
    df_out = data.frame(snp_names=names(idx_hets[!xx]))
   # colnames(df_out) =
    return(df_out)
  }
    #gt_matrix_hets[idx_hets[y],y]
})
names(unique_hets) = colnames(gt_matrix_hets)
unique_het_sites = bind_rows(unique_hets,.id="SAMPLE")
#View(unique_het_sites)


not_unique_hets= lapply(1:ncol(gt_matrix_hets),function(x){ print(x);
  idx_hets = which(gt_matrix_hets[,x] == "0/1" | gt_matrix_hets[,x] == "0/2")
  # print(idx_hets)
  #print(idx_hets)
  if(length(idx_hets) > 0){
    gts=gt_matrix_hets[idx_hets,x]
    #rint(gts)
    xx= sapply(1:length(idx_hets), function(y) { any(gts[y] == gt_matrix_hets[idx_hets[y],-x]) })
    #p#rint(xx)
    #print(sum(xx,na.rm=T))
    #print()
    print(xx)
    print(gts[xx])
    df_out = data.frame(snp_names=names(idx_hets[xx]))
    # colnames(df_out) =
    return(df_out)
  }
  #gt_matrix_hets[idx_hets[y],y]
})
names(not_unique_hets) = colnames(gt_matrix_hets)
not_unique_hets = bind_rows(not_unique_hets,.id="SAMPLE")

aaa = qc_lineages %>% filter(g1 != g2) %>% left_join(not_unique_hets, by=c("strain"="SAMPLE"))
unique(aaa$snp_names)

unique_het_sites = bind_rows(unique_hets,.id="SAMPLE")
#View(unique_het_sites)
het_df = merge(het_df,qc_lineages,by.x="row.names",by.y="sample_name_fasta")


# Example het site in time
het_df %>% filter(`28032_A_G,T` != "./.") %>%ggplot(aes(x=week,fill=`28032_A_G,T`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
het_df %>% filter(`28881_G_A` != "./.") %>%ggplot(aes(x=week,fill=`4711_T_C`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
p1 = het_df %>% filter(`28881_G_A` != "./.") %>%ggplot(aes(x=week,fill=`28881_G_A`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
p2 = het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=atleast_one_het)) + geom_bar() + scale_fill_brewer(palette = "Set1")
cowplot::plot_grid(p1,p2,nrow=2)
#het_df %>% filter(`28881_G_A` != "./.") %>%ggplot(aes(x=week,fill=`28883_G_C`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
### Samples with hets. ###


### count sample with unique hets


#het_df %>% filter(`28032_A_G,T` != "./.") %>%ggplot(aes(x=week,fill=`23403_A_G`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
het_common = het_df[,colnames(het_df) %in% c("week","merged_id","lineage","nextstrain_clade",names( which(het_freq > 0.05)))]

het_common[apply(het_common,1, function(x){any((x== "0/1" | x == "0/2") & x!= "./.")}),]



het_summary = data.frame(snp_names=rownames(gt_matrix_hets),missing=het_missing,freq=het_freq,count=het_count)
het_summary$pos = as.numeric(unlist(lapply(str_split(het_summary$snp_names,"_"), function(x){x[[1]]})))
#het_df$pos = as.numeric(unlist(lapply(str_split(het_df$snp_names,"_"), function(x){x[[1]]})))
het_summary %>% filter(missing < 0.1) %>% mutate(label_snp_new = ifelse(freq > .01,snp_names, "")) %>% ggplot(aes(x=pos,y=freq,label=label_snp_new)) + 
  geom_point() + geom_label_repel() + theme_bw() 

het_summary %>% filter(missing < 0.1) %>% mutate(label_snp_new = ifelse(freq > .01,snp_names, "")) %>% ggplot(aes(x=pos,y=freq,label=label_snp_new)) + 
  geom_point() + geom_label_repel() + theme_bw() 

het_summary %>% filter(missing < 0.1) %>% ggplot(aes(x=count))  + geom_histogram() + theme_bw()

#### So for het sites in each sample we need to check whether it is unique to our sample set ####


####



gt_matrix[gt_matrix == "0/0"] = 0
gt_matrix[gt_matrix == "1/1" | gt_matrix == "2/2"] = 1
gt_matrix[gt_matrix == "./."] = NA
gt_matrix = apply(gt_matrix,2,as.numeric)
#princomp(gt_matrix)
#cor(as.numeric(het_common))
relatedness_matrix = cor(gt_matrix,method="spearman",use ="pairwise.complete.obs")

princomp(cor(gt_matrix,method="spearman",use ="pairwise.complete.obs"))
x = (softImpute(gt_matrix,rank.max = 30))
#gt_matrix %>% gt_matrix
xx = complete(gt_matrix,x)

xxxx = cor(t(xx[(apply(xx,1,function(x){sum(x>0.1)}) > 5) & apply(xx,1,sd) !=0,]))
colnames(xxxx) = snp_df_summary$snp_names[apply(xx,1,function(x){sum(x>0.1)}) > 5 & apply(xx,1,sd) !=0]
pheatmap::pheatmap(xxxx,cluster_cols = F,cluster_rows = F)
xxxxx=  cor(t(xx[(apply(gt_matrix,1,function(x){sum(x>0.1,na.rm=T)}) > 5) & apply(xx,1,sd) !=0,]),use = "pairwise.complete.obs")
colnames(xxxxx) = snp_df_summary$snp_names[apply(gt_matrix,1,function(x){sum(x>0.1,na.rm=T)}) > 5 & apply(xx,1,sd) !=0]
pheatmap::pheatmap(xxxxx,cluster_cols = F,cluster_rows = F)

#rownames(xxxx) = 
pc = data.frame(sample_lineage=colnames(gt_matrix),pc1=x$v[,1],pc2=x$v[,2])
pheatmap::pheatmap(xxxx,cluster_cols = F,cluster_rows = F)
pc %>% left_join(qc_lineages,by=c("sample_lineage"="sample_name_fasta")) %>% ggplot(aes(y=pc1,x=pc2,color=nextstrain_clade)) + geom_point() + theme_bw()
pc %>% left_join(qc_lineages,by=c("sample_lineage"="sample_name_fasta")) %>% ggplot(aes(y=pc1,x=pc2,color=lineage)) + geom_point() + theme_bw()



#### Read in matrix ### 

mat = read.table(("~/Dropbox/COVID19/ncov_new/again/ncov//distance.txt"))
mat_colnames = read.table("~/Dropbox/COVID19/ncov_new/again/ncov/colnames.txt")
mat_rowsnames = read.table("~/Dropbox/COVID19/ncov_new/again/ncov/rownames.txt")
rownames(mat) = mat_colnames[,1]
colnames(mat) = mat_rowsnames[,1]
### OLDEST SAMPLE ####

gt_matrix_raw[,"LA-0193"][gt_matrix_raw[,"LA-0193"] != "0/0"]

mat = mat[),]

###

match_all = data.frame()
for(i in 1:ncol(mat)){
  #mat[,-i]
  idx_in_rows = which(rownames(mat) %in% colnames(mat)[i])
  mat_tmp = mat[-idx_in_rows,]
  name = colnames(mat_tmp)[i]
 # which(mat[-idx_in_rows,i] == min( mat[-idx_in_rows,i]))
  distances = mat_tmp[which(mat_tmp[,i] == min(mat_tmp[,i])),i]
  name_strains=rownames(mat_tmp)[which(mat_tmp[,i] == min(mat_tmp[,i]))]
  xx = data.frame(name=name,distances=paste(distances,sep=",",collapse = ","),name_strains=paste(name_strains,sep=",",collapse = ","))
  xx#distances =mat[-idx_in_rows,i][which(mat[-idx_in_rows,i] == min( mat[-idx_in_rows,i]))]
  
  match_all = rbind(match_all,xx)
}


match_all$close = as.vector(sapply(match_all$distances,function(x){length(str_split(x,",")[[1]])}))

#### Los samples ###

library(data.table)
#install.packages("data.table")
gisaid = fread("~/Downloads/metadata_2020-07-10_23-53.tsv")
california = read.delim("data/qc/california.csv",sep=";")

gisaid = gisaid %>% left_join(california,by=c("strain"="seqName"))
gisaid$date_fix = ymd(gisaid$date)


#### 


gisaid %>% filter(country== "USA") %>% filter(GISAID_clade %in% c("GH","G","S","GR")) 
#qc_lineages$date_fix

california = gisaid %>% select(strain, pangolin_lineage,date, division, country,location,GISAID_clade,clade,QCStatus) %>% filter(division == "California" & QCStatus == "Pass")
california$date = ymd(california$date)
california$week = week(california$date)
p1 =california %>% filter(clade != "-") %>% filter(location == "Los Angeles County") %>% ggplot(aes(x=week, fill=clade)) + geom_bar(width=.9) + xlim(c(0,26)) + scale_fill_brewer(palette = "Set3") + ggtitle("Other LA genomes")
p2 =qc_lineages %>% filter(nextstrain_clade != " -- ") %>% ggplot(aes(x=week,fill=nextstrain_clade)) + geom_bar(width=.9) + xlim(c(0,26)) + scale_fill_brewer(palette = "Set3") + ggtitle("Our LA genomes") + theme()
p3 =california %>% filter(clade != "-") %>%  filter(location != "Los Angeles County") %>% ggplot(aes(x=week, fill=clade)) + geom_bar(width=.9) + xlim(c(0,26)) + scale_fill_brewer(palette = "Set3") + ggtitle("Other californian genomes")
cowplot::plot_grid(p2 +theme(legend.justification = c(0,1)) + theme_bw() ,p1 + theme(legend.justification = c(0,1)) + theme_bw() ,p3 + theme(legend.justification = c(0,1)) + theme_bw(),nrow=3,align = "v")

library(Biostrings)

usa_feb = gisaid %>% filter(country== "USA") %>% filter(month(date_fix) < 3)
xstringset = readDNAStringSet(("data/gisaid/sequences_2020-07-10_23-53.fasta"))
writeXStringSet(xstringset[usa_feb$strain],file="data/gisaid/early_usa.fasta")



names(xstringset)

nextclade = read.delim("data/gisaid/nextclade.csv",sep=";")
nextclade %>% left_join(gisaid,by=c("seqName"="strain")) %>% filter(grepl("20",clade))
  
  mm = read.csv("~/Downloads/metadata.csv")
  
  
  wri
  
### Check for strais with multiple lineages
  
  
### 
  
pango_merge = read.delim("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/outputs/pangolin/lineages_merged.csv",sep="\t")
taxon_list = unlist(lapply(str_split(pango_merge$taxon_p,":"), function(x){x[2]}))
pango_merge$phase = taxon_list
#aaa = pango_merge %>% left_join(qc_lineages,by=c("taxon"="sample_name_fasta"))

pang$p1 = unlist(lapply(str_split(aaa$nextstrain_phased,","), function(x){x[1]}))
aaa$p2 = unlist(lapply(str_split(aaa$nextstrain_phased,","), function(x){x[2]}))

pango_merge$g1 = unlist(lapply(str_split(pango_merge$pangolin_phased,","), function(x){x[1]}))
pango_merge$g2 = unlist(lapply(str_split(pango_merge$pangolin_phased,","), function(x){x[2]}))
gt = pango_merge %>% filter(g1 != g2)

ggg = read.table("~/Downloads/test.txt",stringsAsFactors = F)

#gt = (aaa %>% inner_join(ggg,by=c("uid"="V1")))
library(Biostrings)
ystringset = readDNAStringSet("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/outputs/consensus/merged/phased/sars2/all_3_none.fasta")
ncov= read.delim("~/KRUGLYAK/covid_runs/rerun_all_new_july12/uid_type/outputs/ncov/merged/qc_report.tsv")
###
idx_to_keep = c()
fx = data.frame()

for (seqname in gt$taxon){
  idx_tmp = grep(str_replace(seqname,"-","_"), names(ystringset))
  idx_to_keep = c(idx_to_keep, idx_tmp)
  p1 = gt$g1[gt$taxon == seqname]
  p2 = gt$g2[gt$taxon == seqname]
  fx = rbind(fx,data.frame(seqname, new_names=names(ystringset)[idx_tmp],pangolin_lineage=c(p1,p2)))
}
#fx %>% aaa
ncov.2 = ncov %>% left_join(fx,by=c("strain"="seqname")) %>% filter((strain %in% gt$taxon.y) | strain == "Wuhan/Hu-1/2019" | strain == "Wuhan/WH01/2019") 
ncov.2$new_names[(nrow(ncov.2)-1):nrow(ncov.2)] = ncov.2$strain[(nrow(ncov.2)-1):nrow(ncov.2)]

ncov.2$old_strain = ncov.2$strain
ncov.2$strain = ncov.2$new_names
ncov.2$pangolin_lineage = ncov.2$pangolin_lineage.x
ncov.2$pangolin_lineage[is.na(ncov.2$pangolin_lineage.x)] = ncov.2$pangolin_lineage.y[is.na(ncov.2$pangolin_lineage.x)]
writeXStringSet(c(ystringset[idx_to_keep],xstringset[names(xstringset)[grepl("Wuhan/Hu",names(xstringset))]],xstringset[names(xstringset)[grepl("Wuhan/WH01",names(xstringset))]]), file="../ncov_new/again/ncov/data/merged//het_sites.fasta")
write.table(ncov.2,file="../ncov_new/again/ncov/data/merged/het.tsv",sep="\t",col.names=T,row.names=F,quote=F)


a = read.csv("~/KRUGLYAK/covid_runs/rerun_all_new_july20/downsample/outputs/pangolin/raw.csv")
b = read.csv("~/KRUGLYAK/covid_runs/rerun_all_new_july20/downsample/outputs/pangolin/imputted.csv")
coverage = str_split(a$taxon,"_")
a$coverage = as.numeric(unlist(lapply(coverage,function(x){x[3]})))
a$strain = str_split(a$taxon,"_")
a$strain = as.character(unlist(lapply(a$strain,function(x){x[2]})))

b$strain = str_split(b$taxon,"_")
b$strain = as.character(unlist(lapply(b$strain,function(x){x[2]})))
coverage = str_split(b$taxon,"_")
b$coverage = as.numeric(unlist(lapply(coverage,function(x){x[3]})))
a = a %>% group_by(strain) %>% mutate(lk = (lineage[coverage==max(coverage)])) %>% mutate(correct =lk == lineage)
b = b %>% group_by(strain) %>% mutate(lk = (lineage[coverage==max(coverage)])) %>% mutate(correct =lk == lineage)

#merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(coverage.x),group=strain.x,label=coverage.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") + theme_bw() + theme(text=element_text(size=20)) + geom_label_repel()
#merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(correct.x),group=strain.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") 

merge(a,b,by=1) %>% ggplot(aes(y=UFbootstrap.x,x=UFbootstrap.y,color=factor(correct.x),group=strain.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") 


merge(a,b,by=1) %>% ggplot(aes(y=UFbootstrap.x,x=UFbootstrap.y,color=factor(correct.y),group=strain.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") 


merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(correct.y),group=strain.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") 
merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(correct.x),group=strain.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") 


p1 = merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(correct.x),group=strain.x,label=coverage.x)) + geom_point() + scale_color_brewer(palette = "Set1") + geom_abline() +ylab("Imputation probability")  +xlab("Raw probability") + theme_bw() + theme(text=element_text(size=20))

p2 = merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(correct.y),group=strain.x,label=coverage.x)) + geom_point() + scale_color_brewer(palette = "Set1") + geom_abline() +ylab("Imputation probability")  +xlab("Raw probability") + theme_bw() + theme(text=element_text(size=20))

cowplot::plot_grid(p1,p2,ncol=2)

merge(a,b,by=1) %>% ggplot(aes(y=probability.y,x=probability.x,color=factor(correct),group=strain.x,label=coverage.x)) + geom_point() + scale_color_viridis_d() + geom_abline() +ylab("Strain assignment") + theme_bw() + theme(text=element_text(size=20))


#merge(a,b,by=1)  %>% group_by(strain) %>% mutate(pangolin_lineage)
#merge(a,b,by=1) %>% summarise(sum(lineage.x ==lineage.y))
###

a$imputted = F
b$imputted = T
a = bind_rows(a,b)
coverage = str_split(a$taxon,"_")
a$coverage = as.numeric(unlist(lapply(coverage,function(x){x[3]})))
a %>% filter(coverage <= 10) %>% ggplot(aes(x=probability,fill=factor(coverage))) + geom_histogram() + facet_wrap(~coverage*imputted)

a %>% filter(coverage <= 10) %>% ggplot(aes(x=probability,fill=factor(correct))) + geom_histogram() + facet_wrap(~coverage*imputted)


a %>% ggplot(aes(x=probability,fill=factor(correct))) + geom_histogram() + facet_wrap(~imputted)

### library
