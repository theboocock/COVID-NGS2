library(treeio)
library(ggtree)
library(tidyverse)
library(cowplot)
library(castor)

source("R/utils.R")
source("R/load_datasets.R")


#### Write out tables ###
#write.table(gisaid_phased_filt,file="data/gisaid_all.txt",quote=F,row.names=F,col.names=T, sep="\t")
#gisaid_phased_filt$probability[match(lineages_manual$taxon, gisaid_phased_filt$strain)] = lineages_manual$probability
#gisaid_phased_filt$pangolin_lineage[match(lineages_manual$taxon, gisaid_phased_filt$strain)] = lineages_manual$lineage
#gisaid_phased_filt %>% filter(submitting_lab == "UCLA") %>%   select(strain,pangolin_lineage,probability) %>% write.table("paper/tables/S2.tsv", quote=F,row.names=F,col.names=T,sep="\t")
gisaid_phased_filt %>% filter(submitting_lab != "UCLA") %>%   select(strain,pangolin_lineage,probability) %>% write.table("paper/tables/S3.tsv", quote=F,row.names=F,col.names=T,sep="\t")
gisaid_phased_filt %>% filter(submitting_lab != "UCLA") %>% arrange(date_fix) %>% filter(strain != "USA/CA-ALSR-0513-SAN/2020") %>% filter(grepl("B.1",pangolin_lineage)) %>% filter(country == "USA") %>% head(n=20)  %>% write.table("paper/tables/S7.tsv", quote=F,row.names=F, col.names=T,sep="\t")
gisaid_phased_filt %>% select(-GISAID_clade,-lineage,-probability)%>%  filter(pangolin_lineage == "B.1.43") %>% select(strain,region,division,country,location,submitting_lab) %>% write.table("paper/tables/S4.tsv",quote=F,row.names=F,col.names=T, sep="\t")

#### Evann samples #### 
#evan_in = read.csv("data/evann_samples/evan_samples.csv")
#evan_in %>% filter(patient != "") %>% left_join(gisaid_phased_filt,by=c("seqid"="strain")) %>% select()

### B.1.43 frequenncies

#gisaid_phased_filt %>% select(-GISAID_clade,-lineage,-probability)%>%  filter(pangolin_lineage == "B.1.43") %>% filter(country == "USA") %>% group_by(division) %>% dplyr::summarise(n=n()) %>% mutate(prop=n/sum(n))


## Load the trees into R ##
tree_all = read.tree("data/nextstrain_ncov_north-america_usa_los-angeles_hets_curate_three_timetree.nwk")
diverged_tree = read.tree("data/nextstrain_ncov_north-america_usa_los-angeles_hets_curate_three_tree.nwk")
diverged_table = as_tibble(diverged_tree)
all_la = (gisaid_phased_filt %>% filter(location == "Los Angeles County"))
all_la_strains = all_la$strain
tree_tbl = as_tibble(tree_all)
tree_tbl = tree_tbl %>% mutate(la=ifelse(label %in% all_la$strain,"LA","O"))
tree_tbl_la = tree_tbl %>% filter(label %in% all_la$strain)

diverged_tbl_la = diverged_table %>% filter(label %in% all_la$strain)

parents = tree_tbl$node[grepl("NODE",tree_tbl$label) | grepl("travel",tree_tbl$label)]
local_events = data.frame()
for(node_in in parents){
  tmp = (tree_tbl  %>% child(node_in))
  tmp = tmp[tmp$label %in% tree_all$tip.label,]
  if(nrow(tmp) > 2){
    tree_tbl$la[tree_tbl$node == node_in] = "O"
  }
}
library(castor)
pars = ifelse(tree_all$tip.label  %in% all_la$strain, 1,2)
max_parsimony = asr_max_parsimony(tree = tree_all,tip_states = pars)
node_states = max.col(max_parsimony$ancestral_likelihoods,ties.method = "first")
tree_tbl$la[4069:nrow(tree_tbl)] = ifelse(node_states == 1,"LA","O")
count_intro = 0
intros_nodes = c()
for(node_in in parents){
  is_la = (tree_tbl  %>% filter(node ==node_in))
  is_la = is_la$la
  tmp = (tree_tbl  %>% child(node_in))
  if(is_la =="O"){
    count_intro = count_intro + sum(tmp$la == "LA")
    tmp2 = tmp[tmp$la == "LA",]
    if(nrow(tmp2) > 0){
      intros_nodes = c(intros_nodes,tmp2$node)
    }
  }
}
out_intro = data.frame()
samples_in = c()
for(node_in in intros_nodes){
  tmp = tree_tbl %>% offspring(node_in) 
  if(nrow(tmp) == 0){
    tmp = data.frame(node=node_in,total_la=1,total=1)
    samples_in = c(samples_in, tree_tbl$label[tree_tbl$node == node_in])
  }else{
    {
      samples_in = c(samples_in, tmp$label)
      tmp = tmp %>% dplyr::summarise(node=node_in,total_la = sum(label %in% all_la_strains),total=n())
    }
  }
  out_intro = rbind(out_intro,tmp)
}
out_intro2 = out_intro %>% arrange(-total)
dim(out_intro2)
##complete=T
#i# = 1
# removed = c()
#while(complete){
#  #if(local_event_sort)
#  if(is.na(out_intro2$node[i])){complete=F}
#  removed = c(removed,(idx_rm$node))
#  out_intro2 = out_intro2 %>% filter(!(node %in% idx_rm$node))
#  # print(dim(local_event_tmp))
#  i = i + 1
#}
### Identify clades of strains ###
### read.json tree ###
library(rjson)
tree_annot_all = fromJSON(file = "data/branch_lengths.json")
introduction_clades = data.frame()

for(node_in in out_intro2$node){
  
  node_label = tree_tbl$label[tree_tbl$node == node_in]
  if(grepl("_travel_history",node_label)){
    node_label = str_replace_all(node_label,"_travel_history","")
  }
  date_tree = (tree_annot_all$nodes[[node_label]]$date)
  tmp = tree_tbl %>% offspring(node_in) %>% filter(label %in% all_la_strains) %>% left_join(gisaid_phased_filt,by=c("label"="strain"))
  if(nrow(tmp) == 0){
    tmp = tree_tbl %>%filter(node == node_in) %>% left_join(gisaid_phased_filt,by=c("label"="strain"))   %>% mutate(total_la_strains=1,node_date=date_tree) %>% select(node, total_la_strains, label, pangolin_lineage, node_date) %>% dplyr::rename(strain=label)
  }else{
    tmp = tmp %>% dplyr::summarise(node=node_in,total_la_strains=n(),strain=label[1], pangolin_lineage = names(sort(table(pangolin_lineage),decreasing = T))[1],node_date=date_tree)
  }
  #  tmp$date_tree = ymd()
  introduction_clades = rbind(introduction_clades,tmp)
}
#introduction_clades$week = week(introduction_clades$node_date)
introduction_clades = introduction_clades %>% mutate(node_date_fix = ymd(node_date)) %>% mutate(week=week(node_date_fix), month=month(node_date_fix))
#introduction_early= introduction_clades %>% filter(week <= 16)
#introduction_filt_early= introduction_clades %>% filter(week <= 16) %>% filter(node %in% out_intro2$node[out_intro2$total_la > 1])
### Region coloring ####
list_strains_in_tree = gisaid_phased_filt %>% filter(strain %in% tree_all$tip.label)
tree_tbl %>% filter(!(label %in% gisaid_phased_filt$strain))

tip_list = list()
#names(tip_list)= unique(list_strains_in_tree$region)
#tip_list[[1]] = tree_tbl_la$label
#names(tip_list)[7]= "LA County"
for(i in 1:length(unique(list_strains_in_tree$region))){
  tmp_lis = list_strains_in_tree$strain[list_strains_in_tree$region == unique(list_strains_in_tree$region)[i]]
  
  # tmp_lis = tmp_lis[!tmp_lis %in% tree_tbl_la$label]
  tip_list[[i]] = tmp_lis
}
tip_list2= list()
tip_list2[["North America"]] = tip_list[[1]]
tip_list2[["Asia"]] = tip_list[[6]]
tip_list2[["Europe"]] = tip_list[[5]]
tip_list2[["Other"]] = c(tip_list[[2]],tip_list[[3]], tip_list[[4]])

names(tip_list)= c(unique(list_strains_in_tree$region))
#tip_list[[i + 1]] = tree_tbl_la$label
#names(tip_list)[7]= "LA County"
#groups = tree_all$tip.label

#gisaid_local = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid_strain)]) %>% bind_rows(gisaid_strain)gisaid_local = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid_strain)]) %>% bind_rows(gisaid_strain)
#gisaid_local2 = qc_lineages %>% select(colnames(qc_lineages)[colnames(qc_lineages) %in% colnames(gisaid_strain)]) %>% bind_rows(gisaid)

annotate_regions_df = gisaid_local2[match(tree_all$tip.label,gisaid_local2$strain),]
strains_annotate = gisaid_local$region
names(strains_annotate) = gisaid_local$strain
aa_2222 = groupOTU(tree_all,tip_list2,group_name = "Region")
scale = RColorBrewer::brewer.pal(3,name = "Set1")[-c(2)]
scale[3] = "#377eb8"
scale=c(scale,"grey80")
scale[1] = "#984ea3"
#scale = c("#000000", scale)
##scale = paste(scale,"70",sep="")
#scale = c(scale[1:3],"#000000",scale[4:6])
### 

onea = ggtree(aa_2222,aes(color=Region),mrsd='2020-06-29')+ scale_color_manual(values = scale)+ geom_tippoint(mapping = aes(subset = node %in% tree_tbl_la$node), color="#000000",size=3,shape=17) + theme_tree2() + theme(text=element_text(size=30))  

aa_222_div = groupOTU(diverged_tree,tip_list2,group_name="Region")
alpha = 0.4
sonea = ggtree(aa_222_div,aes(color=Region))+ scale_color_manual(values = scale)+ geom_tippoint(mapping = aes(subset = node %in% tree_tbl_la$node), color="#000000",size=1.5,shape=17) + theme_tree2() + theme(text=element_text(size=24)) 
oneb  =ggtree(aa_2222,aes(color="grey80"),mrsd='2020-06-29')  %>% scaleClade(5450, scale=10)  %>% scaleClade(4202,scale=20) %>% scaleClade(5264,scale=20)  %>%
  scaleClade(5700,scale=20) + scale_color_manual(values = c("grey80")) + theme_tree2()  +
  geom_hilight(node=5450,fill="#377eb8",alpha=alpha) +  geom_hilight(node=4202,fill="#4daf4a",alpha=alpha) +  geom_hilight(node=5264,fill="#984ea3",alpha=alpha) + 
  geom_hilight(node=5700,fill="#ff7f00",alpha=alpha) +theme(text=element_text(size=30),legend.position = "none") + geom_tippoint(mapping = aes(subset = node %in% tree_tbl_la$node),shape=17,size=3,color="#000000") +   geom_nodepoint(mapping = aes(subset = node %in% out_intro2$node[out_intro2$total_la > 1]), color="#e41a1c",size=5,shape=8)
soneb  = ggtree(aa_222_div,aes(color="grey80"))  %>% scaleClade(5450, scale=10)  %>% scaleClade(4202,scale=20) %>% scaleClade(5264,scale=20)  %>%
  scaleClade(5700,scale=20) + scale_color_manual(values = c("grey80")) + theme_tree2()  +
  geom_hilight(node=5450,fill="#377eb8",alpha=alpha) +  geom_hilight(node=4202,fill="#4daf4a",alpha=alpha) +  geom_hilight(node=5264,fill="#984ea3",alpha=alpha) + 
  geom_hilight(node=5700,fill="#ff7f00",alpha=alpha) +theme(text=element_text(size=30),legend.position = "none") + geom_tippoint(mapping = aes(subset = node %in% tree_tbl_la$node),shape=17,size=3,color="#000000") +   geom_nodepoint(mapping = aes(subset = node %in% out_intro2$node[out_intro2$total_la > 1]), color="#e41a1c",size=5,shape=8)
plot_grid(sonea,soneb + theme(legend.position = "none") ,labels="AUTO", label_size = 24) 
ggsave("paper/figures//S4.png",width=16,height=12,dpi=150)

cowplot::plot_grid(onea,oneb,nrow=1,ncol=2,labels="AUTO",label_size = 30)

ggsave("figures/1.png",width=40,height=20,dpi=300)
ggsave("figures/1_b.pdf",width=24,height=18/2,dpi=300)

onea + theme(legend.position = "none")
ggsave("figures/1_a.pdf",width=24/2,height=18/2,dpi=300)
oneb
ggsave("figures/1_b.pdf",width=24/2,height=18/2,dpi=300)


cowplot::plot_grid(sonea,soneb,nrow=1,ncol=2,labels="AUTO",label_size = 30)
ggsave("figures/S4.pdf",width=24/2,height=18/2,dpi=300)

#gisaid_local = gisaid_phased_filt
#gisaid_local$week = week(gisaid_local$date_fix)
#gisaid_local$pangolin_simple = factor(gisaid_local$pangolin_simple)

#gisaid_local = gisaid_local %>% mutate(pangolin_simple2=ifelse(pangolin_lineage %in% need_colors_next,pangolin_lineage,"Other B.1"))
#gisaid_local$pangolin_simple2 = factor(gisaid_local$pangolin_simple2)
need_colors_next = c(names(which(table(introduction_clades$pangolin_lineage) >=3)))
introduction_clades = introduction_clades %>% mutate(pangolin_simple=ifelse(pangolin_lineage %in% need_colors_next,pangolin_lineage,"Other"))
introduction_clades$pangolin_simple = factor(introduction_clades$pangolin_simple)
introduction_clades = introduction_clades %>% mutate(pangolin_simple2=ifelse(pangolin_lineage %in% need_colors_next,pangolin_lineage,"Other B.1"))
introduction_clades$pangolin_simple2 = factor(introduction_clades$pangolin_simple2)
introduction_early= introduction_clades %>% filter(week <= 16)
introduction_filt_early= introduction_clades %>% filter(week <= 16) %>% filter(node %in% out_intro2$node[out_intro2$total_la > 1])


gisaid_local = gisaid_phased_filt %>% mutate(pangolin_simple=factor(ifelse(pangolin_lineage %in% need_colors_next,pangolin_lineage,"Other")))
gisaid_local = gisaid_local %>% mutate(pangolin_simple2=factor(ifelse(pangolin_lineage %in% need_colors_next,pangolin_lineage,"Other B.1")))


#need_colors_next = c(names(which(table(introduction_clades$pangolin_lineage) >=4 )),"B.1.43")
week_filter = 17
color_map = RColorBrewer::brewer.pal(length(need_colors_next) + 2  , name ="Set3")
color_map=color_map[-c(2)]
names(color_map) = c(need_colors_next,"Other")
scale_map2 = scale_fill_manual(name="Lineage",values = color_map)
p1 = introduction_clades %>% filter(week < week_filter)  %>% ggplot(aes(x=node_date_fix,fill=pangolin_simple)) + geom_histogram(binwidth = 7) + scale_map2 + 
  theme_bw() + scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20")))  + theme(text=element_text(size=24))  + ggtitle("LA County Introductions") + xlab("Date") + ylab("Count")
p2 = gisaid_local %>% filter(location=="King County")  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle("King County") + xlab("Date") +
  xlab("Date") + ylab("Count")
locale = "King County"
p2 = gisaid_local %>% filter(location==locale)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")
locale = "Santa Clara County"
p3 = gisaid_local %>% filter(location==locale)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")
locale= "New York"
p4 = gisaid_local %>% filter(division==locale)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")
locale = "Europe"
p5 = gisaid_local %>% filter(region==locale)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")
locale = "Italy"
p6= gisaid_local %>% filter(country==locale)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")
locale = "Iran"
p7 = gisaid_local %>% filter(locale==country)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")

locale = "Wuhan"
p8 = gisaid_local %>% filter(locale==location)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")
locale = "Louisiana"
p9 = gisaid_local %>% filter(locale==division)  %>% filter(week < week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple))  + 
  scale_map2+  geom_histogram(binwidth =7)  +theme_bw() + theme(text=element_text(size=24)) +scale_x_date(limits=as.Date(c("2020-01-15","2020-04-20"))) + ggtitle(locale) + xlab("Date") +
  xlab("Date") + ylab("Count")


pf = cowplot::plot_grid(p2 + theme(legend.position = "none"),p3 +  theme(legend.position = "none"),p4 +  theme(legend.position = "none"),p5 +  theme(legend.position = "none"),p6+  theme(legend.position = "none"),p7+theme(legend.position = "none"),
                        p8+theme(legend.position = "none"), p9+theme(legend.position = "none"),ncol=2,align="hv")

cowplot::plot_grid(onea,oneb,p1,pf,nrow=2,ncol=2,labels="AUTO",label_size = 30)

ggsave("figures/1.png",width=24,height=18,dpi=300)
#ggsave("figures/1_a.png",width=24,height=18,dpi=300)
#pf
#ggsave("figures/1_d.pdf",width=24/2,height=18/2,dpi=300)
#gisaid_
#p1
#ggsave("figures/1_c.pdf",width=24/2,height=18/2,dpi=300)


#color_map = RColorBrewer::brewer.pal(length(need_colors_next) + 2  , name ="Set3")
#names(color_map) = c(need_colors_next,"Other")
#color_map = RColorBrewer::brewer.pal(length(need_colors_next) + 2  , name ="Set3")
#color_map = color_map[-c(2)]
#names(color_map) = c(need_colors_next,"Other B.1")
#scale_map3 = scale_fill_manual(name="Lineage",values = color_map)

#pitro = introduction_clades %>% filter(week >= week_filter) %>% ggplot(aes(x=node_date_fix,fill=(pangolin_simple2)))+ geom_histogram(binwidth = 7)  + scale_map3 + theme_bw()  + scale_x_date(limits=as.Date(c("2020-04-20","2020-06-30"))) + theme(text=element_text(size=20))  + xlab("") +
#  ylab("Count") + ggtitle("LA County Introductions")
#gisaid_local$week = week(gisaid_local$date_fix)
#pp_count = gisaid_local %>%  filter(location == "Los Angeles County") %>% filter(week >= week_filter) %>% ggplot(aes(x=date_fix,fill=pangolin_simple2)) + geom_histogram(binwidth = 7)  + scale_map3 + theme_bw()+ theme(text=element_text(size=20)) + xlab("Date") + ylab("Count") + ggtitle("LA County Genomes")
#pp_count_leg = get_legend(pp_count)
#get
#pg = cowplot::plot_grid(pitro + theme(legend.position = "none"),pp_count  +  theme(legend.position = "none"),nrow=2,labels = "AUTO",label_size = 20)
#library(cowplot)
#plot_grid(pg, pp_count_leg,rel_widths = c(1,.4))
#ggsave("paper/figures//S6.png", width=8,height=8,dpi=150)

#cowplot::plot_grid(p1,p2)
#ggsave("figures/2_new.png",width = 10,height=8)
# geom_tippoint(subset=.(group==1))
three_a = gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% ggplot(aes(x=date_fix,fill=Lineage)) + geom_histogram(binwidth = 14)  + 
  theme_bw()  + theme(text=element_text(size=30))   + xlab("Date") + ylab("Count") +  scale_fill_manual(values = c("grey75","#e41a1c")) +   scale_x_date(limits=as.Date(c("2020-01-15","2020-06-30"))) 
la = gisaid_phased_filt %>% filter(location == "Los Angeles County")
month_numeric <-lubridate::month(la$date_fix)
month_label <- lubridate::month(la$date_fix,label=T)
month_numeric = month_numeric[(!duplicated(month_label))]
month_label = unique(month_label)


la_sum = la %>%  filter(location == "Los Angeles County") %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(week,Lineage) %>% dplyr::summarise(n=n()) %>% mutate(freq = n/sum(n),total=(n)) %>% filter(Lineage == "B.1.43") 
p_la = la %>% filter(location == "Los Angeles County") %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(month,Lineage) %>% dplyr::summarise(n=n()) %>% mutate(freq = n/sum(n),total=(n)) %>% filter(Lineage == "B.1.43") %>%
  ggplot(aes(y=freq,x=month,label=total)) + geom_line() + geom_point() + geom_text_repel(size=10) + scale_x_continuous(labels=month_label,breaks=month_numeric) + theme_bw() + theme(text=element_text(size=30)) + xlab("Month") + ylab("B.1.43 frequency") + ylim(c(0,.38))
p_la = la %>% filter(submitting_lab == "UCLA") %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(month,Lineage) %>% dplyr::summarise(n=n()) %>% mutate(freq = n/sum(n),total=(n)) %>% filter(Lineage == "B.1.43") %>%
  ggplot(aes(y=freq,x=month,label=total)) + geom_line() + geom_point() + geom_text_repel(size=10) + scale_x_continuous(labels=month_label,breaks=month_numeric) + theme_bw() + theme(text=element_text(size=30)) + xlab("Month") + ylab("B.1.43 frequency") + ylim(c(0,.38))


states_map = map_data("state")
b141 = gisaid_phased_filt %>%  mutate(division=tolower(division)) %>% filter(pangolin_lineage == "B.1.43") %>% group_by(division) %>% dplyr::summarise(n=n()) 
states <- map_data("state")
colnames(states)[5] <- "State"
library(stringi)
library(plyr)
states$State <- (stri_trans_totitle(states$State))
df <- data.frame(state.x77,
                 State = state.name,
                 state_lower=tolower(state.name),
                 Abbrev = state.abb,
                 Region = state.region,
                 Division = state.division
)  
library(stringi)
library(scatterpie)
df2 <- merge(states,df,by="State")
df2 <- df2[order(df2$order),]
mid_range <- function(x) mean(range(x,na.rm=TRUE))
centres <- ddply(df2, .(state_lower),
                 colwise(mid_range,.(lat,long,Population)))

b141 = b141 %>% right_join(centres,by=c("division"="state_lower")) %>% inner_join(df,by=c("division"="state_lower"))
b142 = b141 %>% filter(!is.na(n))
pa = ggplot() + geom_map(data=b141,aes(fill=n,map_id=division), map=states_map)  +   expand_limits(x = states_map$long, y = states_map$lat) + geom_text(data=b141,aes(label=n,x=long,y=lat),size=8,color="white") + scale_fill_viridis_c(end=.5,name="B.1.43 Count")  + coord_map() +
  geom_label_repel(data=b142,aes(label=State,x=long,y=lat),size=4,nudge_x = .5,nudge_y=3,segment.color = 'transparent') + theme_map() + theme(text=element_text(size=20))
county_map = map_data("county")
freq_cali_wash = gisaid_phased_filt  %>%mutate(location=ifelse(location == "San Diego County","San Diego",location)) %>% mutate(location=ifelse(location == "San Francisco County","San Francisco",location))%>%  mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% filter(division == "Washington" | division == "California") %>% group_by(division,location,Lineage) %>% 
  dplyr::summarise(total=n()) %>% group_by(division, location) %>% dplyr::mutate(freq=total/sum(total)) %>% dplyr::ungroup() %>% dplyr::mutate(count_lower= tolower(unlist(lapply(str_split(location," County"), function(x){x[[1]][1]}))))
county = county_map  %>% filter(region == "california" | region == "washington")
colnames(states)[5] <- "County"
county$County <- (stri_trans_totitle(county$subregion))
mid_range <- function(x) mean(range(x,na.rm=TRUE))
centres_county <- ddply(county, .(subregion),
                        colwise(mid_range,.(lat,long)))
freq_cali_data = freq_cali_wash %>% dplyr::inner_join(centres_county,by=c("count_lower"="subregion")) %>% mutate(count_upper=toupper(count_lower)) 
freq_cali_data= pivot_wider(freq_cali_data,id_cols=c("count_lower","lat","long"),names_from="Lineage",values_from=c("freq","total"))
freq_cali_data =freq_cali_data %>% mutate(County=str_to_title(count_lower))
aaa =freq_cali_wash %>% dplyr::full_join(centres_county,by=c("count_lower"="subregion")) %>% mutate(count_upper=toupper(count_lower)) 
aaa = aaa %>% mutate(state_lower = tolower(division)) %>%dplyr::full_join(centres,by=c("state_lower"="state_lower")) %>% filter(location == "")
aaa = pivot_wider(aaa,id_cols=c("division","lat.y","long.y"),names_from="Lineage",values_from=c("freq","total")) %>% mutate(label=paste(division,"Unspecified"))
cali_county = county %>% filter(region == "california")
pb = ggplot() + geom_polygon(data=cali_county,aes(group=subregion,x=long,y=lat),fill="white",color="black") +   geom_scatterpie(data=freq_cali_data,aes(group=count_lower,x=long,y=lat,r=(total_B.1.43)/80),cols=c("freq_Other","freq_B.1.43"))  +
  geom_label_repel(data=subset(freq_cali_data,!is.na(total_B.1.43)),aes(label=County,x=long,y=lat),nudge_x = -2.5,nudge_y=-.5,size=5) +  geom_label_repel(data=subset(freq_cali_data,!is.na(total_B.1.43)),aes(label=total_B.1.43,x=long,y=lat),nudge_x = .8,nudge_y=.75,size=8,segment.color = 'transparent')  + coord_map() + theme_map() +
  geom_scatterpie(data=aaa,aes(group=division,x=long.y,y=lat.y,r=(total_B.1.43)/80),cols=c("freq_Other","freq_B.1.43")) +   geom_label_repel(data=aaa,aes(label=label,x=long.y,y=lat.y),nudge_x = 2,nudge_y = -.4,size=5) +   theme(legend.position="none") + 
  geom_label_repel(data=aaa,aes(label=total_B.1.43,x=long.y,y=lat.y),nudge_x = .4,nudge_y = .5,segment.color = 'transparent',size=8)  + xlim(c(min(cali_county$long),max(cali_county$long))) + ylim(c(min(cali_county$lat),max(cali_county$lat))) + scale_fill_manual(values = c("#e41a1c","grey75"))
washington = county %>% filter(region == "washington")

pc = ggplot() + geom_polygon(data=county,aes(group=subregion,x=long,y=lat),fill="white",color="black") +   geom_scatterpie(data=freq_cali_data,aes(group=count_lower,x=long,y=lat,r=(total_B.1.43)/80),cols=c("freq_Other","freq_B.1.43"))  +
  geom_label_repel(data=subset(freq_cali_data,!is.na(total_B.1.43)),aes(label=County,x=long,y=lat),nudge_x = -2.5,nudge_y=-.5,size=5) +  geom_label_repel(data=subset(freq_cali_data,!is.na(total_B.1.43)),aes(label=total_B.1.43,x=long,y=lat),nudge_x = .8,nudge_y=.2,size=8,segment.color = 'transparent')  + coord_map() + theme_map() +
  geom_scatterpie(data=aaa,aes(group=division,x=long.y,y=lat.y,r=(total_B.1.43)/80),cols=c("freq_Other","freq_B.1.43")) +   geom_label_repel(data=aaa,aes(label=label,x=long.y,y=lat.y),nudge_x = 4,nudge_y = -.2,size=5) + 
  geom_label_repel(data=aaa,aes(label=total_B.1.43,x=long.y,y=lat.y),nudge_x = .6,nudge_y = .5,segment.color = 'transparent',size=8)  +  xlim(c(min(washington$long),max(washington$long))) + ylim(c(min(washington$lat),max(washington$lat)+.5)) + scale_fill_manual(values = c("#e41a1c","grey75")) +   theme(legend.position="none")


pa =  ggplot() + geom_polygon(data=washington,aes(group=subregion,x=long,y=lat),fill="white",color="black",size=1) + coord_map()+ theme_map()
pbb = ggplot() + geom_polygon(data=cali_county,aes(group=subregion,x=long,y=lat),fill="white",color="black",size=1) + coord_map()+ theme_map()
ggplot() + geom_polygon(data=states,aes(group=County,x=long,y=lat),fill="white",color="black") + coord_map()+ theme_map()


#plot_grid(pa,pbb)
#ggsave("figures/cali_wash.png",width=16,height=12,dpi=300)

#pg = plot_grid(pb + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),pc + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
#plot_grid(pa + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),pg,three_a,p_la,labels="AUTO",ncol=2,rel_heights = c(1,.8),label_size=30)
plot_grid(three_a,p_la,labels="AUTO",ncol=2,rel_heights = c(1,.8))

ggsave("figures/2cd.pdf",width=16,height=8,dpi=150)
#ggplot() + geom_map(data=counnty)
#gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(month,Lineage) %>% summarise(n=n()) %>% mutate(freq = n/sum(n),total=sum(n)) %>% filter(Lineage == "B.1.43") %>%
#ggplot(aes(y=freq,x=month,label=total)) + geom_line() + geom_point() + geom_text_repel() + theme_bw()

#locales = c("Whatcom County","San Diego","Contra Costa County","Skagit County", "Imperial County","Orange County CA","San Luis Obispo County","Ventura County")

#p1 = gisaid_phased_filt %>% filter(location %in% locales) %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(month,Lineage) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq = n/sum(n),total=(n)) %>% filter(Lineage == "B.1.43") %>%
#  ggplot(aes(y=freq,x=month,label=total)) + geom_line() + geom_point() + geom_text_repel() + scale_x_continuous(labels=month_label,breaks=month_numeric) + theme_bw() + theme(text=element_text(size=30)) + xlab("Month") + ylab("B.1.43 frequency")  + ylim(c(0,.25))
#p2 = gisaid_phased_filt %>% filter(location %in% locales) %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(month,location,Lineage) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq = n/sum(n),total=(n)) %>% filter(Lineage == "B.1.43") %>%
#  ggplot(aes(y=freq,x=month,label=total, group=location)) + geom_line() + geom_point() + geom_text_repel(size=10) + scale_x_continuous(labels=month_label,breaks=month_numeric) + theme_bw() + theme(text=element_text(size=30)) + xlab("Month") + ylab("B.1.43 frequency") + facet_wrap(~location)

#gisaid_phased_filt %>% filter(location %in% locales) %>%
 # mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>%
#  group_by(month,location,Lineage) %>%
#  dplyr::summarise(n=n()) %>%
 # dplyr::mutate(freq = n/sum(n),total=(n)) 

#plot_grid(p1,p2,rel_widths = c(0.8,1.4),labels="AUTO",label_size = 24)
#ggsave(width=16,height=12,file="paper/figures//S8.png",dpi=150)
#library(forcats)
#gisaid_phased_filt %>% filter(location %in% locales) %>% mutate(Lineage=factor(ifelse(pangolin_lineage == "B.1.43","B.1.43","Other"),levels=c("Other","B.1.43"))) %>% group_by(month,location,Lineage) %>% summarise(n=n()) %>% mutate(freq = n/sum(n),total=sum(n)) %>% filter(Lineage == "B.1.43") %>%
#  ggplot(aes(y=freq,x=month,label=total, group=location)) + geom_line() + geom_point() + geom_text_repel() + scale_x_continuous(labels=month_label,breaks=month_numeric) + theme_bw() + theme(text=element_text(size=30)) + xlab("Month") + ylab("B.1.43 frequency") + facet_wrap(~location)

#p1 = gisaid_local %>% filter(location %in% locales) %>% mutate(b_1_43=ifelse(pangolin_lineage == "B.1.43",T,F))  %>% mutate(Location=location) %>% filter(b_1_43)%>% ggplot(aes(x=date_fix,fill=Location)) + geom_histogram(binwidth = 14)  + 
#  theme_bw()  + theme(text=element_text(size=20))  + ggtitle("B.1.43") + xlab("Date") + ylab("Count") + scale_fill_brewer(palette = "Set2") +   scale_x_date(limits=as.Date(c("2020-01-15","2020-06-30"))) 

#p2 = gisaid_local %>% filter(location %in% locales) %>% mutate(b_1_43=ifelse(pangolin_lineage == "B.1.43",T,F)) %>% mutate(Location=location)  %>% filter(!b_1_43)%>% ggplot(aes(x=date_fix,fill=location)) + geom_histogram(binwidth = 14)  + 
#  theme_bw()  + theme(text=element_text(size=20))  + ggtitle("Other lineages") + xlab("Date") + ylab("Count") + scale_fill_brewer(palette = "Set2") +   scale_x_date(limits=as.Date(c("2020-01-15","2020-06-30"))) 
#leg = get_legend(p1)
#library(forcats)

#pg= plot_grid(p1 + theme(legend.position = "none"),p2 + theme(legend.position = "none"),labels="AUTO",nrow=2,rel_widths = c(1,1),label_size = 24)
#plot_grid(pg, leg,label_size = 24,rel_widths = c(3,1))
#ggsave("paper/figures/S7.png",width=12,height=8,dpi=150)
#a = read.delim("../ncov_new/again/ncov/results/4710.txt", header=F)
#aaa =  a %>% inner_join(gisaid_phased_filt,by=c("V1"="strain"))
#table(aaa$pangolin_lineage)
#aaa %>% filter(pangolin_lineage == "B.1.1") %>% select(V1,date,pangolin_lineage) %>% dplyr::rename(strain=V1) %>% write.table("paper/tables/S10.tsv",sep="\t",col.names=T,row.names=F)

#View(aaa)
#snp = read.delim("~/Dropbox/COVID19/ncov_new/again/ncov/results/sequence-diagnostics.tsv")
#aaa %

#g# = snp %>% inner_join(gisaid_phased_filt,by=c("strain"="strain"))

#g %>% filter(location == "Whatcom County")
#### Create a nextstrain  build for B.1.43 and B.1.1 focused builds ### 