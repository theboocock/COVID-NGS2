library(ggtree)
library(treeio)
library(tidyverse)
tree_b_1_43 = read.tree("data/nextstrain_ncov_la_all_b_one_two_timetree.nwk")
tree_b_1_43_div = read.tree("data/nextstrain_ncov_la_all_b_one_two_tree.nwk")

tree_b_1_43_tbl = as_tibble(tree_b_1_43)
tree_b_1_43_div_tbl= as_tibble(tree_b_1_43_div)
library(rjson)
tree_json = fromJSON(file="~/Dropbox/COVID19/ncov_new/again/ncov/results/la_all_b_one_two//branch_lengths.json")
nt_muts = fromJSON(file="~/Dropbox/COVID19/ncov_new/again/ncov/results/la_all_b_one_two/nt_muts.json")
nodes_mut =tree_b_1_43_tbl %>% offspring(1745) %>% filter(grepl("NODE",label))
nodes_mut = nodes_mut$label
aa = nt_muts$nodes[nodes_mut]
aa = lapply(aa,function(x){data.frame(x$muts)})  
aa = bind_rows(aa,.id = "label") 
aa = aa[aa$x.muts  %in% c("T2244C","A20374G","G25012A","G4657A","C379A","C337T"),]
tree_b_1_43_tbl %>% filter(label == "LA_0193_1")
aa %>% inner_join(tree_b_1_43_tbl,by=c("label"="label"))
tree_df = tree_b_1_43_tbl %>% offspring(2019) %>% inner_join(gisaid_phased_filt,by=c("label"="strain"))
nodes_mut =tree_b_1_43_tbl %>% offspring(2019) %>% filter(grepl("NODE",label))
nodes_mut = nodes_mut$label

tree_df = tree_df  %>% mutate(cali_or_wash = ifelse(division != "California" & division != "Washington","Other",division))
#label_location = list()

tree_df = tree_df  %>% mutate(cali_or_wash_location = ifelse(division != "California" & division != "Washington","Other",location)) 

locations_list = list()
j = 1
for(location_filt in unique(tree_df$cali_or_wash_location)){
  tree_tmp = tree_df %>% filter(cali_or_wash_location == location_filt) %>% select(label)
  locations_list[[j]] =  tree_tmp$label
  j = j + 1
}
names(locations_list) =  unique(tree_df$cali_or_wash_location)
division_list=list()
j = 1 
for(location_filt in unique(tree_df$cali_or_wash)){
  
  tree_tmp = tree_df %>% filter(cali_or_wash == location_filt) %>% select(label)
  division_list[[j]] = tree_tmp$label
  j = j + 1
}
names(division_list) = unique(tree_df$cali_or_wash)

tree_b_1_43 = groupOTU(tree_b_1_43,locations_list,group_name = "Location")
tree_b_1_43_div =  groupOTU(tree_b_1_43_div,locations_list,group_name = "Location")

tree_b_1_43 = groupOTU(tree_b_1_43,division_list,group_name = "State")
tree_b_1_43_div =  groupOTU(tree_b_1_43_div,division_list,group_name = "State")


b1 = gisaid_phased_filt %>% select(-GISAID_clade,-lineage,-probability)%>%  filter(pangolin_lineage == "B.1.43")
tree_b_1_43_div_tbl %>% filter(label %in% b1$strain) %>% ggplot(aes(x=branch.length)) + geom_histogram()

nodes = aa %>% inner_join(tree_b_1_43_tbl,by=c("label"="label"))
aaa = ggtree(tree_b_1_43_div,color="grey50")

#View(aaa$data%>% inner_join(b1,by=c("label"="strain")) %>% filter())
#tree_df = tree_b_1_43_tbl %>% offspring(2106) %>% inner_join(gisaid_phased_filt,by=c("label"="strain"))
#nodes_mut =tree_b_1_43_tbl %>% offspring(2106) %>% filter(grepl("NODE",label))
#nodes_mut = nodes_mut$label

#tree_df = tree_df  %>% mutate(cali_or_wash = ifelse(division != "California" & division != "Washington","Other",division))
#label_location = list()

#tree_df = tree_df  %>% mutate(cali_or_wash_location = ifelse(division != "California" & division != "Washington","Other",location)) 

aaa= ggtree(tree_b_1_43_div)
aaa$date = aaa$data %>% full_join(aa,by=c("label"="label"))
ggtree(tree_b_1_43,color="grey50") %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2)   + scale_color_brewer(palette = "Dark2") + geom_nodelab(aes(label=x.muts))
ggtree(tree_b_1_43_div,color="grey50") %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2)   + scale_color_brewer(palette = "Dark2")+ geom_label(aes(x=branch,label=x.muts))

ggtree(tree_b_1_43,color="grey50")  %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  +
  geom_tippoint(mapping = aes(subset = node %in% tree_df$node[tree_df$location != ""],color=Location),size=2)   + scale_color_brewer(palette = "Paired")+ geom_label(aes(x=branch,label=x.muts))+ theme(text=element_text(size=14))
aa %>% inner_join(tree_b_1_43_tbl,by=c("label"="label"))
ggsave("paper/figures//S9.png",width=16*.75, height=12*.75,dpi=150)

ggtree(tree_b_1_43_div,color="grey50")  %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node[tree_df$location != ""],color=Location),size=2)   + scale_color_brewer(palette = "Paired") # geom_label(aes(x=branch,label=x.muts))


ggtree(tree_b_1_43_div,color="grey50")  %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node[tree_df$location != ""],color=Location),size=2)   + scale_color_brewer(palette = "Paired") +
  geom_highlight(node=2119,fill="#E41A1C",alpha=0.2) + geom_highlight(node=2141,fill="#377EB8",alpha=0.2)  + geom_highlight(node=2167,fill="#4DAF4A",alpha=0.2)  + geom_highlight(node=2211,fill="#984EA3",alpha=0.2)  + geom_highlight(node=2289,fill= "#FF7F00",alpha=0.2) + geom_label(aes(x=branch,label=x.muts))


tree_b_1_43_tbl %>% offspring(2019) %>% inner_join(gisaid_phased_filt,by=c("label"="strain")) %>% select(division)


ggtree(tree_b_1_43_div,color="grey50")  %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node[tree_df$location != ""],color=Location),size=2)   + scale_color_brewer(palette = "Paired") # geom_label(aes(x=branch,label=x.muts))
ggsave(file="figures/test_div.pdf",width=16,height=10,dpi=150)


b_143_strain_list = tree_b_1_43_tbl %>% offspring(2019) %>% inner_join(gisaid_phased_filt,by=c("label"="strain")) %>% filter(location == "Los Angeles County") %>% filter(label != "USA/CA-CSMC31/2020") 
b_143_strain_lockdown = b_143_strain_list %>% filter(week < 17 & week > 8)
xstringset = readDNAStringSet(("data/masked.fasta"))
writeXStringSet(xstringset[b_143_strain_lockdown$label],filepath = "data/b_1_43_lockdown.fasta")
b_143_strain_lockdown %>% select(label, date_fix) %>% write.table("data/b_1_43_lockdown_dates.txt",sep="\t",col.names=FALSE, row.names=FALSE, quote=FALSE)



a = read.delim("../ncov_new/again/ncov/results/4710.txt", header=F)

b_1_1 = gisaid_phased_filt %>% filter(location == "Los Angeles County" | location == "Orange County CA") %>% filter(strain %in% a$V1) 
focus_on_lockdown = gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% filter(week < 17 & week > 8)

writeXStringSet(xstringset[focus_on_lockdown$strain],filepath = "data/focus_on_lockdown.fasta")
focus_on_lockdown %>% select(strain, date_fix) %>% write.table("data/focus_on_lockdown_dates.txt",sep="\t",col.names=FALSE, row.names=FALSE, quote=FALSE)


writeXStringSet(xstringset[b_1_1$strain],filepath = "data/b111.fasta")
b_1_1 %>% select(strain,date_fix) %>% write.table(file = "data/b111.dates.txt", quote = F,sep="\t",row.names=F,col.names = F)

la_strains = gisaid_phased_filt %>% filter(location == "Los Angeles County")

writeXStringSet(xstringset[la_strains$strain],filepath = "data/la.fasta")
la_strains %>% select(strain,date_fix) %>% write.table(file = "data/la.dates.txt", quote = F,sep="\t",row.names=F,col.names = F)

a = gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% group_by(week) %>% mutate(n=n()) %>% filter(n <= 9)
 b =gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% group_by(week) %>% mutate(n=n()) %>% filter(n > 9) %>% dplyr::sample_n(10,replace = FALSE) 
subset_genomes = bind_rows(a,b) 
subset_genomes  %>% ungroup %>% select(strain, date_fix) %>% write.table(file="data/subset_dates.txt",quote=FALSE,row.names = FALSE, col.names = FALSE, sep ="\t")
strains_sub = subset_genomes$strain
(gisaid_phased_filt %>% filter(location == "Los Angeles County") %>% group_by(week) %>% dplyr::sample_n(10,replace = TRUE))
i
#b_143_strain_list = gisaid_phased_filt %>%  filter(location == "Los Angeles County") %>% select(strain)
writeXStringSet(xstringset[strains_sub],filepath = "data/strains_sub.fasta")

writeXStringSet(xstringset[b_143_strain_list$label],filepath = "data/b143.fasta")
###
tree_b_1_43_tbl %>% inner_join(gisaid_phased_filt,by=c("label"="strain")) %>% filter(location == "Los Angeles County")%>% filter(label != "USA/CA-CSMC31/2020")   %>% select(label,date_fix) %>% write.table(quote = F,file="data/beast_dates.txt",sep="\t",col.names = F,row.names = F)

tree_b_1_43_tbl %>% offspring(2019) %>% inner_join(gisaid_phased_filt,by=c("label"="strain")) %>%  filter(label != "USA/CA-CSMC31/2020")  %>% filter(location == "Los Angeles County") %>% select(label,date_fix)  %>% mutate(date_fix=(date_fix)) %>% write.table(quote = F,file="data/beast_dates.txt",sep="\t",col.names = F,row.names = F)
a= (0.8e-03/sqrt(5e-4))^2
b = sqrt(5e-4)^2/(0.8e-03)
#### #### B_1_43 #### ##### 



ggsave(file="figures/test_div.pdf",width=16,height=12,dpi=150)


ggsave(file="figures/test_div.pdf",width=12,height=6,dpi=150)

#tree_b_1_43_tbl %>% filter(label =)



ggtree(tree_b_1_43_div,color="grey50") %>% scaleClade(node=2106,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2) + scale_color_brewer(palette = "Dark2")
ggtree(tree_b_1_43,color="grey50") %>% scaleClade(node=2106,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2)   + scale_color_brewer(palette = "Dark2")


ggtree(tree_b_1_43_div,color="grey50") %>% scaleClade(node=2106,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=Location)) + scale_color_brewer(palette = "Dark2")
ggtree(tree_b_1_43,color="grey50") %>% scaleClade(node=2106,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2)   + scale_color_brewer(palette = "Dark2")

### 
b_1_43_json = fromJSON(file = "~/Dropbox/COVID19/ncov_new/again/ncov/results/north-america_usa_los-angeles_hets_curate_three/")


tree_b_1_43
