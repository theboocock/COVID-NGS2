tree_phased = read.tree("data/nextstrain_ncov_north-america_usa_los-angeles_hets_curate_two_timetree.nwk")

#tree_phased= read.tree("~/Downloads/nextstrain_ncov_north-america_usa_los-angeles_hets_curate_tree (1).nwk")

tree_tbl_phased = as.tibble(tree_phased)
tree_tbl_phased[tree_tbl_phased$label == "LA_0193_1",]
tree_tbl_phased[tree_tbl_phased$label == "LA_0193_2",]


tree_tbl_phased[tree_tbl_phased$label == "LA_0172_1",]
tree_tbl_phased[tree_tbl_phased$label == "LA_0172_2",]

a = read.delim("../ncov_new/again/ncov/results/4710.txt", header=F)
aaa =  a %>% inner_join(gisaid_phased_filt,by=c("V1"="strain"))
table(aaa$location)

county= map_data("county")
county = county %>% filter(region == "california")

freq_cali_wash = gisaid_phased_filt  %>%mutate(location=ifelse(location == "San Diego County","San Diego",location)) %>% mutate(location=ifelse(location == "San Francisco County","San Francisco",location))%>%  mutate(Lineage=factor(ifelse(strain %in%aaa$V1, "B.1.1(4711T)","Other"),levels=c("B.1.1(4711T)","Other"))) %>% filter(division == "Washington" | division == "California") %>% group_by(division,location,Lineage) %>% 
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


paa =ggplot() + geom_polygon(data=county, aes(group=subregion,y=lat,x=long),fill="white",color="black") + coord_map() + geom_text(data=freq_cali_data[!is.na(freq_cali_data$`total_B.1.1(4711T)`),],aes(label=`total_B.1.1(4711T)`,y=lat,x=long),size=10) + geom_label_repel(data=freq_cali_data[!is.na(freq_cali_data$`total_B.1.1(4711T)`),],aes(label=County,y=lat,x=long),nudge_x =.5,nudge_y=1.5,size=10) + theme_map() 
ggsave("figures/3A.png",width=8,height=12,dpi=150)
# 487 and 653
la_strains = c(all_la_strains,tree_tbl_phased$label[grep("^LA",tree_tbl_phased$label)])

group1 = tree_phased$tip.label[tree_phased$tip.label %in% la_strains] 
group2 =  tree_phased$tip.label[!tree_phased$tip.label %in% la_strains] 

tree_tbl_phased_only = tree_tbl_phased %>% filter(label %in% qc_lineages_phased2[!is.na(qc_lineages_phased2$lineage.y),]$strain)

la_genomes_all = tree_tbl_phased %>% filter(label %in% la_strains)
#highlight_these = qc_lineages_phased2[!is.na(qc_lineages_phased2$lineage.y),]
tree_phased2 =groupOTU(tree_phased, list(group1,group2))
p2 = ggtree(tree_phased2,mrsd="2020-06-29",color="grey80") %>% scaleClade(6306,scale=40)  %>% scaleClade(6808,scale=100)  + geom_tippoint(mapping = aes(subset = node %in% la_genomes_all$node),shape=17,size=4) + scale_color_manual(values = c("grey80","grey80")) +
  geom_highlight(node=6306,fill="#377eb8",alpha=0.4) +   geom_highlight(node=6808,fill="#4daf4a",alpha=0.4)  + theme_tree2() + theme(text=element_text(size=20))  
p2
ggsave("figures/3b.pdf",height=8,width=12,dpi=300)

tree_phased4 =groupOTU(tree_phased3, list(group1,group2))