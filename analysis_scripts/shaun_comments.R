### ###



superspreader_list = list()
i = 1
for(node_in in out_intro2$node){
  node_label = tree_tbl$label[tree_tbl$node == node_in]
  if(grepl("_travel_history",node_label)){
    node_label = str_replace_all(node_label,"_travel_history","")
  }
  date_tree = (tree_annot_all$nodes[[node_label]]$date)
  tmp = tree_tbl %>% offspring(node_in) %>% filter(label %in% all_la_strains) %>% left_join(gisaid_phased_filt,by=c("label"="strain")) %>% filter(submitting_lab != "UCLA" & !grepl("Cedars",originating_lab)) %>% left_join(gisaid_phased_filt,by=c("label"="strain"))
  if(nrow(tmp) > 0){
  superspreader_list[[i]] = tmp
  i =i +1
  }
#  if(nrow(tmp) > 0){
#    tmp = tree_tbl %>%filter(node == node_in) %>% left_join(gisaid_phased_filt,by=c("label"="strain"))   ) %>% select(node, total_la_strains, label, pangolin_lineage, node_date) %>% dplyr::rename(strain=label)
 ## }else{
  #  tmp = tmp %>% dplyr::summarise(node=node_in,total_la_strains=n(),strain=label[1], pangolin_lineage = names(sort(table(pangolin_lineage),decreasing = T))[1],node_date=date_tree)
#  }
  #  tmp$date_tree = ymd()
 # introduction_clades = rbind(introduction_clades,tmp)
#  i =i +1
}
