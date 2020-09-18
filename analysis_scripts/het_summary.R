source("R/utils.R")
source("R/load_datasets.R")
library(ggrepel)
merged_vcf_filt = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20///uid_type/outputs/quasi_species/merged/all.vcf")
merged_vcf_filt = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20///uid_type/outputs/quasi_species/merged/all_from_merged.vcf")

test_genomes = c("LA_0171","LA_0260","LA_0081","LA_0258","LA_0193","LA_0239")

m1  = cbind(merged_vcf_filt@fix[,2],merged_vcf_filt@fix[,4],merged_vcf_filt@fix[,5],merged_vcf_filt@gt[,colnames(merged_vcf_filt@gt) == "LA-0061"]
)
#merged_vcf_filt@gt[,colnames(merged_vcf_filt@gt) == "LA-0061"]
m1[!grepl("0/0",m1[,4]),]
merged_vcf_filt3 = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/test2.vcf")
qc_input_raw = read.delim("~/KRUGLYAK/covid_runs/rerun_all_new_july20/uid_type/outputs/qc_report/qc_report.tsv")
##### 
gt_matrix = get_gt_matrix(qc_input_raw,merged_vcf_filt3,qc_input_raw[qc_input_raw$uid == "U0011",]$sample)
ar_list = list() 
i = 1
b_1_1= c(4711,28881,28882,28883)
b_1_43= c(1059,379,25563)
for(sample in qc_input_raw[qc_input_raw$uid == "U0011",]$sample){
  ar_matrix = get_gt_matrix(qc_input_raw,merged_vcf_filt3, sample)
  ar_matrix$snp_name = rownames(ar_matrix)
  ar_matrix$ar = as.numeric(ar_matrix$ar)
  ar_list[[i]] = ar_matrix
  i = i + 1
}
names(ar_list) =qc_input_raw[qc_input_raw$uid == "U0011",]$sample
het_sites_ar = bind_rows(ar_list,.id="sample")
positions = as.factor(as.numeric(unlist(lapply(str_split(het_sites_ar$snp_name,"_"), function(x){x[[1]][1]}))))
het_sites_ar$positions = positions
het_sites_ar %>% filter(sample == "LA-0171")
#het_sites_ar %>% filter(positions %in% c(b_1_1,b_1_43)) %>% mutate(ar_b_1_43=ifelse(positions %in% b_1_43,ar,1-ar)) %>% filter(dp!=0)%>% 
#  ggplot(aes(y=ar_b_1_43,x=factor(pos),group=sample,label=dp)) + geom_point() + geom_line()+ facet_wrap(~sample) + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
 # geom_text_repel() + theme_bw() + theme(text=element_text(size=24))
het_sites_ar %>% filter(positions %in% c(b_1_1,b_1_43)) %>% mutate(takara=ifelse(sample %in% 135,paste("TAKARA",sample,sep="_"),paste("NEB",sample,sep="_")))%>% mutate(ar_b_1_43=ifelse(positions %in% b_1_43,ar,1-ar)) %>% filter(dp!=0) %>%
  ggplot(aes(y=ar_b_1_43,x=factor(positions),group=sample,label=dp,color=takara)) + geom_point() + geom_line() + geom_text_repel() + theme_bw()


het_sites_ar %>% mutate(takara=ifelse(sample %in% 135,paste("TAKARA",sample,sep="_"),paste("NEB",sample,sep="_")))%>% filter(dp >= 2) %>%
  ggplot(aes(y=ar,x=factor(positions),group=sample,label=dp,color=takara)) + geom_point() + geom_line() + geom_text_repel() + theme_bw()

gt_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
gt_matrix_raw = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
ar_list = list() 
i = 1
for(sample in qc_input2$sample_name_fasta){
  ar_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt, sample)
  ar_matrix$snp_name = rownames(ar_matrix)
  ar_matrix$ar = as.numeric(ar_matrix$ar)
  ar_list[[i]] = ar_matrix

  i = i + 1
}
names(ar_list) = qc_input2$sample_name_fasta
het_sites_ar = bind_rows(ar_list,.id="sample")
positions = as.factor(as.numeric(unlist(lapply(str_split(het_sites_ar$snp_name,"_"), function(x){x[[1]][1]}))))
het_sites_ar$positions = positions
### LA-00171
het_28881 = het_sites_ar[het_sites_ar$positions == 28881,][het_sites_ar[het_sites_ar$positions == 28881,]$ar > .1 & het_sites_ar[het_sites_ar$positions == 28881,]$ar <.9,]


dim(het_28881) 
(plate_map %>% inner_join(het_28881,by=c("sample"="sample")))
snp_names = het_sites_ar %>% filter(gt == "0/1") %>% filter(sample == "LA-0171")
test_genomes = c("LA-0171","LA-0260","LA-0081","LA-0258","LA-0193","LA-0239","LA-0181")
het_sites_ar %>% filter(sample %in% test_genomes) %>% filter(snp_name %in% snp_names)
plate_map = read.csv("plate_map_1.csv")
plate_map2 = read.csv("data/plate_map_remap2.csv")

plate_map2 %>% inner_join(qc_lineages,by=c("uid"="uid"))  %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 4711) %>% ggplot(aes(y=ar,x=sample_name_fasta)) + geom_point() + ylim(c(0,1))
aaa = plate_map2 %>% inner_join(qc_lineages,by=c("uid"="uid"))  %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 4711) 
#qc_input_raw[qc_input_raw$sample == 34,
qc_input_raw[qc_input_raw$sample == 345,]
qc_input_raw[qc_input_raw$sample == 346,]
qc_input_raw[qc_input_raw$sample == 347,]

qc_input_raw[qc_input_raw$sample == 364,]
qc_input_raw[qc_input_raw$sample == 349,]
qc_input_raw[qc_input_raw$sample == 350,]


plate_map %>% filter(uid %in% qc_lineages$uid[qc_lineages$g1 != qc_lineages$g2])
het_sites_ar %>% filter(sample %in% test_genomes) %>% mutate(alt_cov=ar*dp) %>% filter(snp_name %in% snp_names$snp_name) %>% ggplot(aes(y=ar,x=positions,color=sample))   + geom_point() + facet_wrap(~sample)
((plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 28882)) %>% ggplot(aes(y=ar,x=z,label=coverage_mean_dedup)) + geom_point() + geom_label_repel()

p1 = ((plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(positions == 4711)) %>% ggplot(aes(y=ar,x=y,label=coverage_mean_dedup)) + geom_point(size=ar) + geom_label_repel()
p2 = ((plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 4711)) %>% ggplot(aes(y=ar,x=z,label=coverage_mean_dedup)) + geom_point(size=ar) + geom_label_repel()
p3= ((plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 4711)) %>% ggplot(aes(y=ar,x=y,label=merged_id)) + geom_point() + geom_label_repel()
p4= ((plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 4711)) %>% ggplot(aes(y=ar,x=z,label=merged_id)) + geom_point() + geom_label_repel()


p5 = ((plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(pos == 28882)) %>% ggplot(aes(y=ar,x=z,label=pangolin_lineage)) + geom_point() + geom_label_repel()


p1 = ((plate_map %>% inner_join(qc_input2,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample"))  %>% filter(pass) %>% filter(positions == 4711)) %>% filter(lib_type=="NEB") %>% ggplot(aes(y=ar,x=y,label=merged_id)) + geom_point() + geom_label_repel()



p2 = ((plate_map %>% inner_join(qc_input2,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample"))%>% filter(positions == 4711)) %>% filter(lib_type=="NEB") %>% mutate(ar_all = ifelse(ar < 0.96 & ar > .05, T,F)) %>% ggplot(aes(y=z,x=y,label=lineage,color=ar_all)) + geom_point() + geom_label_repel()
p3 = ((plate_map %>% inner_join(qc_input2,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample"))%>% filter(positions == 4711)) %>% filter(lib_type=="NEB")  %>% mutate(ar_all = ifelse(ar < 0.96 & ar > .05, T,F))%>%  ggplot(aes(y=z,x=y,label=merged_id,color=ar_all)) + geom_point() + geom_label_repel()
p4 = ((plate_map %>% inner_join(qc_input2,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample"))%>% filter(positions == 4711)) %>% filter(lib_type=="NEB") %>% mutate(ar_all = ifelse(ar < 0.96 & ar > .05, T,F))%>%  ggplot(aes(y=z,x=y,label=uid,color=ar_all)) + geom_point() + geom_label_repel()


duplicate_sent_genomes  = read.csv("data/Duplicate Sent Genomes_Deidentified.csv")
dd = duplicate_sent_genomes %>% inner_join(qc_lineages,by=c("id"="uid")) %>% group_by(group) %>% filter(pass) %>% mutate(n=n())  %>% filter(n >=2) 

qc_lineages %>% filter(lib_type == "NEB") %>% mutate(diff_lineages=ifelse(g1 != g2, "phased","unphased"))  %>% mutate(label=ifelse(merged_id %in% aaa$merged_id, sample_name_fasta,"")) %>% ggplot(aes(y=sars2_percent_dedup,x=ct,color=diff_lineages,label=label)) + geom_point()  + geom_text_repel()


for (pairs in unique(dd$group)){
  aa = pivot_wider(het_sites_ar %>% filter(sample %in% dd$sample_name_fasta[dd$group == pairs]),id_cols=snp_name, names_from = sample, values_from= gt)
  print(aa[aa[,3] != aa[,2],])
}



#(plate_map %>% inner_join(qc_lineages,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(positions == 4711) %>% 
 aaa =  plate_map %>% full_join(qc_input2,by=c("uid"="uid")) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample")) %>% filter(positions == 4711) %>%inner_join(het_summary_sample,by=c("sample_name_fasta"="sample_id"))
#aaa$het_sample_count
 aaa %>% mutate(greater_than_4 = ifelse(het_sample_count > =2,"bad","good")) %>% filter()
 
aaa[aaa$het_sample_count >=4,] %>% select(merged_id,  uid, ct,het_sample_count, date_fix, lib_type) %>% write.csv(file="data/long_sample_list_august19th.csv",row.names=F)
plot 

cowplot::plot_grid(p1,p2,p3,p4)
plot_list=list()
j = 1
all_het_sites = rbind(singleton_het_sites,het_found_in_one,not_unique_hets)
for(i in 1:length(qc_lineages$sample_name_fasta)){
  sample_in = qc_lineages$sample_name_fasta[i]
  tmp_df = het_sites_ar %>% filter(sample ==sample_in)  %>% filter(gt =="0/1") %>%  mutate(ar_folded=ifelse(ar > 0.5,1-ar,ar))
  if(nrow(tmp_df) > 3){
    p_all = tmp_df %>% mutate(SE=sqrt(ar_folded * (1-ar_folded)/dp)) %>% mutate(lower=ar_folded - 2 * SE, upper= ar_folded + 2 * SE)  %>% ggplot(aes(y=ar_folded,x=factor(positions),group=sample,label=dp)) +
     geom_point()+ geom_line() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,.5)) + theme_bw() + geom_text_repel()  + ylab("Minor Allele Frequency") + xlab("Variant position") + ggtitle(sample_in) +
       theme(axis.text.x=element_text(angle=90, hjust=1))  + theme(text=element_text(size=18))
    plot_list[[j]] = p_all
    j = j + 1
  }
}
pos = as.numeric(unlist(lapply(str_split(ar_matrix$snp_name,"_"), function(x){x[[1]][1]})))
het_sites_ar$pos = pos
het_sites_ar %>% filter(sample %in% qc_lineages_phased2$sample_name_fasta[!is.na(qc_lineages_phased2$lineage.y)]) %>% filter(gt == "0/1")  %>% group_by(snp_name) 
qc_lineages_phased$ordered_lineages = apply(cbind(qc_lineages_phased$g1,qc_lineages_phased$g2), 1, function(x){paste(sort(x),sep=",",collapse = ",")})

bbb =  ((plate_map %>% inner_join(qc_input2,by=c("uid"="uid"))) %>% inner_join(het_sites_ar,by=c("sample_name_fasta"="sample"))  %>% filter(positions == 4711)) %>% filter(lib_type=="NEB")
p1 = het_sites_ar %>% filter(snp_name %in% names(table(not_unique_hets$snp_names)[table(not_unique_hets$snp_names) >2]))  %>% filter(positions %in% c(b_1_1,b_1_43))%>% filter(sample %in% bbb$sample_name_fasta) %>%
left_join(lineages_manual,by=c("sample"="sample"))  %>%
  mutate(ar_b_1_43=ifelse(positions %in% b_1_43,ar,1-ar))%>% ggplot(aes(y=ar_b_1_43,x=factor(positions),group=sample,label=dp)) + geom_point() + geom_line()+ facet_wrap(~sample,ncol=8,nrow=6) + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  geom_text_repel() + theme_bw() + theme(text=element_text(size=24)) + scale_color_brewer(palette = "Set1",name="Lineage pair") + ylab("B.1.43 Allele frequency") +  theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,1)) + xlab("Position")
#ggsave("")
ggsave("qc_plate_long.png",width=16,height=12)
het_sites_ar %>% filter(snp_name %in% names(table(not_unique_hets$snp_names)[table(not_unique_hets$snp_names) >2]))  %>% filter(pos %in% c(b_1_1,b_1_43))%>% filter(sample %in% qc_lineages_phased2$sample_name_fasta[!is.na(qc_lineages_phased2$lineage.y)]) %>%
  left_join(lineages_manual,by=c("sample"="sample")) %>% 
  mutate(ar_b_1_43=ifelse(pos %in% b_1_43,ar,1-ar))%>% ggplot(aes(y=ar_b_1_43,x=factor(pos),group=sample,label=dp)) + geom_point() + geom_line()+ facet_wrap(~sample) + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  geom_text_repel() + theme_bw() + theme(text=element_text(size=24)) + scale_color_brewer(palette = "Set1",name="Lineage pair") + ylab("B.1.43 Allele frequency") +  theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,1)) + xlab("Position")


a=qc_lineages %>% filter(sample_name_fasta %in% qc_lineages_phased2$sample_name_fasta[is.na(qc_lineages_phased2$lineage.y)])  %>% filter(pangolin_lineage == "B.1.43") %>% head(n=2)
b=qc_lineages %>% filter(sample_name_fasta %in% qc_lineages_phased2$sample_name_fasta[is.na(qc_lineages_phased2$lineage.y)])  %>% filter(pangolin_lineage == "B.1.1") %>% head(n=2)
a = rbind(a,b)

p2 = het_sites_ar %>% filter(snp_name %in% names(table(not_unique_hets$snp_names)[table(not_unique_hets$snp_names) >2])) %>% filter(pos %in% c(b_1_1,b_1_43)) %>% filter(sample %in% a$sample_name_fasta) %>% 
  mutate(ar_b_1_43=ifelse(pos %in% b_1_43,ar,1-ar))%>% ggplot(aes(y=ar_b_1_43,x=factor(pos),group=sample,label=dp)) + geom_point() + geom_line()+ facet_wrap(~sample) + theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  geom_text_repel() + theme_bw() + theme(text=element_text(size=24)) + scale_color_brewer(palette = "Set1",name="Lineage pair") + ylab("B.1.43 Allele frequency") +  theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,1)) + xlab("Position")

plot_grid(p1,p2,labels="AUTO",label_size = 20)
ggsave("paper/figures//S9.png",height=16,width=16,dpi=300)
#het_sites_ar %>% filter(snp_name %in% names(table(not_unique_hets$snp_names)[table(not_unique_hets$snp_names) >2])) %>% filter(sample %in% qc_lineages_phased2$sample_name_fasta[!is.na(qc_lineages_phased2$lineage.y)]) %>% 
#  mutate(ar_folded=ifelse(ar > 0.5,1-ar,ar))%>% ggplot(aes(y=ar,x=snp_name,group=sample)) + geom_point()+ geom_line() + facet_wrap(~sample) + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,1))


genomes_in=read.csv("data/evann_samples/evan_samples.csv")

bad_plate  = aa %>% select(merged_id,ct,uid,date_fix)
bad_plate$other_id = str_replace_all(bad_plate$merged_id,"_","-")

aa = bad_plate %>% mutate(ifelse(other_id %in% genomes_in$seqid,"early_samples_flu_retest","later_samples_flu_retest"))
write.table(aa, file="data/flu_retest.csv",quote=F, row.names = F)

#qc_lineages %>% filter(g1==g2) %>% inner_join(genomes_in,by=c("sample_name_fasta"="seqid"))



gt_matrix[gt_matrix == "0/1" | gt_matrix == "0/2"] = "./."
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
genotypes_df %>% filter(`23403_A_G` != "./.") %>%ggplot(aes(y=ct,x=`23403_A_G`)) + geom_boxplot() + scale_fill_brewer(palette = "Set1") + theme_bw()
spike_filt = genotypes_df %>% filter(`23403_A_G` != "./.") 
summary(lm(spike_filt$ct[spike_filt$library_type == "neb"] ~ spike_filt$`23403_A_G`[spike_filt$library_type == "neb"]))
p1=het_summary %>% left_join(snp_df_summary,by=c("snp_names"="snp_names")) %>%  mutate(label_snp_new = ifelse(freq.x> .01 | freq.y > 0.05,snp_names, "")) %>% ggplot(aes(y=freq.x,x=pos.x,label=label_snp_new)) + geom_point() + geom_label_repel() +ggtitle("Het frequency")
p2 = het_summary %>% left_join(snp_df_summary,by=c("snp_names"="snp_names")) %>%  mutate(label_snp_new = ifelse(freq.x> .01 | freq.y > 0.05,snp_names, "")) %>% mutate(freq_y_name=ifelse(freq.y > 0.5,1-freq.y,freq.y)) %>% ggplot(aes(y=freq_y_name,x=pos.x,label=label_snp_new)) + geom_point() + geom_label_repel() + ggtitle("Allele frequency")
cowplot::plot_grid(p1,p2,nrow=2)


merged_vcf_filt = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20///uid_type/all_samples_new_filt.vcf")

gt_matrix_hets = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
gt_matrix_hets = gt_matrix_hets[(apply(gt_matrix_hets,1,function(x){any(x== "0/1" | x== "0/2")})),]
het_df = as.data.frame(t(gt_matrix_hets))
het_freq =  apply(gt_matrix_hets,1, function(x){(2*sum(x == "0/1" | x== "0/2"))/(2*length(x[x!="./."]))})
het_missing =  apply(gt_matrix_hets,1, function(x){sum(x == "./.")/length(x)})
het_count = apply(gt_matrix_hets,1, function(x){(sum(x == "0/1" | x== "0/2"))})

het_sample_count= apply(gt_matrix_hets,2,function(x){sum(x == "0/1" | x == "0/2")})
atleast_one_het = het_sample_count >= 1
het_summary_sample = data.frame(het_sample_count=het_sample_count,atleast_one=atleast_one_het,sample_id=colnames(gt_matrix_hets)) 

p1 = het_summary_sample %>% mutate(merged_label= ifelse(het_sample_count > 1, sample_id,"")) %>% ggplot(aes(x=het_sample_count,label=merged_label)) + geom_histogram() + theme_bw() + xlab("Intra-patient variant sites per patient") + ylab("Count") + theme(text=element_text(size=20))
plot_grid(p1)
ggsave("paper/figures/S5.png",dpi=300,width=8,height=6)

p4 = het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=atleast_one_het)) + geom_bar() + scale_fill_brewer(palette = "Set1", name="Samples with more than\n1 intra-patient variant")  + xlab("Collection Week")+ theme_bw() + theme(text=element_text(size=20))
p5 = het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=lib_type,fill=atleast_one_het)) + geom_bar() + scale_fill_brewer(palette = "Set1", name="Samples with more than\n1 intra-patient variant")  + xlab("Library Type")+ theme_bw() + theme(text=element_text(size=20))
#p6 =  het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=factor(het_sample_count)))+ geom_bar()  + xlab("Collection Week")+ theme_bw() + theme(text=element_text(size=20)) 
p6 = het_summary_sample %>% mutate(het_sample_count2=ifelse(het_sample_count>1 & het_sample_count <5,"2-4",het_sample_count)) %>%
  mutate(het_sample_count2=ifelse(het_sample_count>1 & het_sample_count <5,"2-4",het_sample_count)) %>% 
  mutate(het_sample_count2=ifelse(het_sample_count>4,"5-16",het_sample_count)) %>%
  left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=factor(het_sample_count2))) + geom_bar(width=.9) + scale_fill_brewer(name="Number of intra-patient variable sites",palette = "Set2")  + ylab("Count")  + xlab("Collection Week")+ theme_bw() + theme(text=element_text(size=20))

legend = get_legend(p4)
pf = plot_grid(p4 + theme(legend.position = "none"),p5 + theme(legend.position = "none"), labels = "AUTO",ncol=2,label_size = 20)
plot_grid(pf,legend,rel_widths = c(1,.3),p6,labels=c("","","C"),label_size = 20)
ggsave("paper/figures//S6.png",width=16, height=12,dpi=300)

unique_hets= lapply(1:ncol(gt_matrix_hets),function(x){
  idx_hets = which(gt_matrix_hets[,x] == "0/1" | gt_matrix_hets[,x] == "0/2")
  if(length(idx_hets) > 0){
    gts=gt_matrix_hets[idx_hets,x]
    xx= sapply(1:length(idx_hets), function(y) { any(gts[y] == gt_matrix_hets[idx_hets[y],-x]) })
    xxx = sapply(1:length(idx_hets), function(y) { 
      tmp = gt_matrix_hets[idx_hets[y],-x]
      tmp = tmp[tmp != "./."] 
      length(table(tmp)) == 1 
      })
    df_out = data.frame(snp_names=names(idx_hets[!xx]))
    df_out = data.frame(snp_names=names(idx_hets[!xx & xxx]))
    return(df_out)
  }
})
names(unique_hets) = colnames(gt_matrix_hets)
singleton_het_sites = bind_rows(unique_hets,.id="SAMPLE")
unique_hets= lapply(1:ncol(gt_matrix_hets),function(x){ 
  idx_hets = which(gt_matrix_hets[,x] == "0/1" | gt_matrix_hets[,x] == "0/2")
  if(length(idx_hets) > 0){
    gts=gt_matrix_hets[idx_hets,x]
    xx= sapply(1:length(idx_hets), function(y) { any(gts[y] == gt_matrix_hets[idx_hets[y],-x]) })
    xxx = sapply(1:length(idx_hets), function(y) { 
      tmp = gt_matrix_hets[idx_hets[y],-x]
      tmp = tmp[tmp != "./."] 
      length(table(tmp)) == 1 
    })
    df_out = data.frame(snp_names=names(idx_hets[!xx]))
    df_out = data.frame(snp_names=names(idx_hets[!xx & !xxx]))
    return(df_out)
  }
})
names(unique_hets) = colnames(gt_matrix_hets)
het_found_in_one = bind_rows(unique_hets,.id="SAMPLE")
not_unique_hets= lapply(1:ncol(gt_matrix_hets),function(x){ 
  idx_hets = which(gt_matrix_hets[,x] == "0/1" | gt_matrix_hets[,x] == "0/2")
  if(length(idx_hets) > 0){
    gts=gt_matrix_hets[idx_hets,x]
    xx= sapply(1:length(idx_hets), function(y) { any(gts[y] == gt_matrix_hets[idx_hets[y],-x]) })
    df_out = data.frame(snp_names=names(idx_hets[xx]))
    return(df_out)
  }
})

names(not_unique_hets) = colnames(gt_matrix_hets)
not_unique_hets = bind_rows(not_unique_hets,.id="SAMPLE")
sample_ins = het_summary_sample$sample_id[!het_summary_sample$atleast_one ]
for(snp in unique(not_unique_hets$snp_names)){
  snp_idx = which(rownames(gt_matrix_hets) == snp) 
  col_idx = gt_matrix_hets[,which(colnames(gt_matrix_hets) %in% sample_ins)]
}

gt_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
gt_matrix_raw = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
ar_list = list() 
i = 1
for(sample in qc_lineages$sample_name_fasta){
  ar_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt, sample)
  ar_matrix$snp_name = rownames(ar_matrix)
  ar_matrix$ar = as.numeric(ar_matrix$ar)
  ar_list[[i]] = ar_matrix
  
  i = i + 1
}
names(ar_list) = qc_lineages$sample_name_fasta
het_sites_ar = bind_rows(ar_list,.id="sample")
positions = as.factor(as.numeric(unlist(lapply(str_split(het_sites_ar$snp_name,"_"), function(x){x[[1]][1]}))))
het_sites_ar$positions = positions

het_sites_ar %>% filter(snp_name %in% unique(not_unique_hets$snp_names)) %>% group_by(snp_name) %>% filter(sample %in% sample_ins)  %>% summarise(ar=sum(gt[gt!="./."]=="0/0"), ac=sum(gt[gt!="./."] != "0/0"),af=sum(gt[gt!="./."] != "0/0")/length(gt[gt!="./."]))
het_sites_ar %>% filter(snp_name %in% unique(not_unique_hets$snp_names)) %>% group_by(snp_name) %>% filter(sample %in% sample_ins)  %>% summarise(reference_count=sum(gt[gt!="./."]=="0/0"), alternative_count=sum(gt[gt!="./."] != "0/0"),alt_freq=sum(gt[gt!="./."] != "0/0")/length(gt[gt!="./."]),missing_cout=sum(gt=="./.")) %>%
  write.table(file="paper/tables/S2.tsv",quote=F,row.names=F,col.names=T,sep="\t")

het_sites_ar %>% filter(snp_name %in% singleton_het_sites$snp_names) %>% group_by(snp_name) %>% filter(sample %in% sample_ins)  %>% summarise(reference_count=sum(gt[gt!="./."]=="0/0"), alternative_count=sum(gt[gt!="./."] != "0/0"),alt_freq=sum(gt[gt!="./."] != "0/0")/length(gt[gt!="./."]), missing_count=sum(gt=="./."))  %>%
  write.table(file="paper/tables//S3.tsv",quote=F,row.names = F,col.names=T,sep="\t")

het_df = merge(het_df,qc_lineages,by.x="row.names",by.y="sample_name_fasta")


# Example het site in time
het_df %>% filter(`28032_A_G,T` != "./.") %>%ggplot(aes(x=week,fill=`28032_A_G,T`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
het_df %>% filter(`28881_G_A` != "./.") %>%ggplot(aes(x=week,fill=`4711_T_C`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
p1 = het_df %>% filter(`28881_G_A` != "./.") %>%ggplot(aes(x=week,fill=`28881_G_A`)) + geom_bar() + scale_fill_brewer(palette = "Set1")
p2 = het_summary_sample %>% left_join(qc_lineages, by=c("sample_id"="sample_name_fasta")) %>% ggplot(aes(x=week,fill=atleast_one_het)) + geom_bar() + scale_fill_brewer(palette = "Set1")
cowplot::plot_grid(p1,p2,nrow=2)
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
#### Let's create the pairwise boxplot####
### Let's use the original all.vcf
merged_vcf_filt = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20///uid_type/outputs/quasi_species/merged/all.vcf")

gt_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
#gt_matrix_raw = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
ar_list = list() 
i = 1
for(sample in qc_lineages$sample_name_fasta){
  ar_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt, sample)
  ar_matrix$snp_name = rownames(ar_matrix)
  ar_matrix$ar = as.numeric(ar_matrix$ar)
  ar_list[[i]] = ar_matrix
  
  i = i + 1
}
names(ar_list) = qc_lineages$sample_name_fasta
het_sites_ar = bind_rows(ar_list,.id="sample")
positions = as.factor(as.numeric(unlist(lapply(str_split(het_sites_ar$snp_name,"_"), function(x){x[[1]][1]}))))
het_sites_ar$positions = positions
extract_het_sites_between = het_sites_ar %>%filter(sample %in% qc_lineages_phased2$sample_name_fasta[is.na(qc_lineages_phased2$lineage.y)] ) #%>% filter(sample %in% qc_lineages$sample_name_fasta[grepl("B.1",qc_lineages$pangolin_phased)])
pairs2= c()
for(i in 1:100){
  
 s1 = sample(unique(extract_het_sites_between$sample)  ,2,replace=F)
 # lin = qc_lineages  %>% filter(sample_name_fasta == s1) %>% select(pangolin_lineage)
#  lineages = lin$pangolin_lineage
 # tmp_ex = extract_het_sites_between %>% filter(sample %in% qc_lineages$sample_name_fasta[qc_lineages$pangolin_lineage != lineages])
#  s2 = sample(unique(tmp_ex$sample)  ,1,replace=F)
  tmp_ex = extract_het_sites_between %>% filter(sample %in% s1)
 # tmp_ex2 = extract_het_sites_between %>% filter(sample%in%s1)
  #g1 = bind_rows(tmp_ex,tmp_ex2)
  snp_differences = tmp_ex %>% filter(gt != "0/1" & gt != "./." & gt != "0/2")  %>% pivot_wider(names_from=sample,values_from=gt,id_cols = snp_name)
  pairs2 = c(pairs2,sum(snp_differences[,2] != snp_differences[,3],na.rm=T))
  #print(snp_differences[which(snp_differences[,2] != snp_differences[,3]),])
}

pairs_phased = c()

for(sample_in  in unique(qc_lineages_phased2$sample_name_fasta[!is.na(qc_lineages_phased2$lineage.y)])){
  pairs_phased = c(pairs_phased,nrow(het_sites_ar %>% filter(sample %in% sample_in) %>% filter(gt == "0/1" | gt == "0/2")))
}

pairs_phased2  = c()
sample_df = data.frame(sample=qc_lineages_phased2$sample_name_fasta[is.na(qc_lineages_phased2$lineage.y)]) %>% inner_join(het_summary_sample[het_summary_sample$atleast_one,],by=c("sample"="sample_id")) 
for(sample_in  in unique(sample_df$sample)){
  pairs_phased2 = c(pairs_phased2,nrow(het_sites_ar %>% filter(sample %in% sample_in) %>% filter(gt == "0/1" | gt == "0/2")))
}

d2 = data.frame(Phased="Coinfection",pairs=pairs_phased)
d1 = data.frame(Phased="Random pairs (different lineages)",pairs=pairs)
d4= data.frame(Phased="Random pairs", pairs=pairs2)
d3 = data.frame(Phased="Same lineage",pairs=pairs_phased2)
d1 = rbind(d1,d2)
d1 = rbind(d1,d3)
d1 = rbind(d1,d4)
#d1$Phased = levels(d1$Phased)
d1$Phased = factor(d1$Phased,levels=c("Coinfection","Same lineage", "Random pairs","Random pairs (different lineages)"))
d1 %>% ggplot(aes(y=pairs,x=Phased)) + geom_boxplot(width=0.5,outlier.shape = NA) + geom_jitter(height = 0,width=.25) +theme_bw() +xlab("") + ylab("Pairwise differences") + theme_bw() + theme(text=element_text(size=26))
ggsave("figures/S9.png",width=16,height=12,dpi=150)





### #### 

