library(ggtree)
library(treeio)
library(tidyverse)
tree_b_1_43 = read.tree("data/nextstrain_ncov_la_all_b_one_two_timetree.nwk")
tree_b_1_43_div = read.tree("data/nextstrain_ncov_la_all_b_one_two_tree.nwk")

tree_b_1_43_tbl = as_tibble(tree_b_1_43)
tree_b_1_43_div_tbl= as_tibble(tree_b_1_43_div)
library(rjson)
tree_json = fromJSON(file="data/tree_json_b_1_43.json")
nt_muts = fromJSON(file="data/nt_muts.json")
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
#tree_b_1_43_div_tbl %>% filter(label %in% b1$strain) %>% ggplot(aes(x=branch.length)) + geom_histogram()

nodes = aa %>% inner_join(tree_b_1_43_tbl,by=c("label"="label"))
aaa = ggtree(tree_b_1_43_div,color="grey50")

#View(aaa$data%>% inner_join(b1,by=c("label"="strain")) %>% filter())
#tree_df = tree_b_1_43_tbl %>% offspring(2106) %>% inner_join(gisaid_phased_filt,by=c("label"="strain"))
#nodes_mut =tree_b_1_43_tbl %>% offspring(2106) %>% filter(grepl("NODE",label))
#nodes_mut = nodes_mut$label

#tree_df = tree_df  %>% mutate(cali_or_wash = ifelse(division != "California" & division != "Washington","Other",division))
#label_location = list()

#tree_df = tree_df  %>% mutate(cali_or_wash_location = ifelse(division != "California" & division != "Washington","Other",location)) 

#aaa= ggtree(tree_b_1_43_div)
#aaa$date = aaa$data %>% full_join(aa,by=c("label"="label"))
#ggtree(tree_b_1_43,color="grey50") %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2)   + scale_color_brewer(palette = "Dark2") + geom_nodelab(aes(label=x.muts))
#ggtree(tree_b_1_43_div,color="grey50") %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  + geom_tippoint(mapping = aes(subset = node %in% tree_df$node,color=State),size=2)   + scale_color_brewer(palette = "Dark2")+ geom_label(aes(x=branch,label=x.muts))

ggtree(tree_b_1_43,color="grey50",mrsd = "2020-06-19")  %<+% aa %>% scaleClade(node=2019,scale = 30)+ theme_tree2()  +
  geom_tippoint(mapping = aes(subset = node %in% tree_df$node[tree_df$location != ""],color=Location),size=2)   + scale_color_brewer(palette = "Paired")+ geom_label(aes(x=branch,label=x.muts))+ theme(text=element_text(size=14))
#aa %>% inner_join(tree_b_1_43_tbl,by=c("label"="label"))
ggsave("paper/figures//S5.png",width=16*.75, height=12*.75,dpi=150)

