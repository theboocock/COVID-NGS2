local_events2 = c()
parents2 = tree_tbl$node
for(node_in in parents2){
  print(node_in)
  tmp = (tree_tbl  %>% offspring(node_in) %>% filter(!grepl("NODE",label)))
  if(nrow(tmp) ==0){
    la_in  = tree_tbl$la[tree_tbl$node == node_in]
    if(la_in == "LA"){
      tmp = data.frame(node=node_in, total_la=1,total=1)  
    }else{
      tmp = data.frame(node=node_in, total_la=0,total=1)  
    }
  }else{
    tmp =  tmp %>% summarise(node=node_in,total_la = sum(label %in% all_la_strains ), total=n())
  }
  if(tmp$total_la >=1 & tmp$total_la/tmp$total == 1)
  {
    local_events2 = rbind(local_events2,tmp)
  }
}

local_event_tmp2 = local_events2 %>% arrange(-total)
#local_event_tmp = local_events %>% arrange(-total)
complete=T
i = 1

while(complete){
  #if(local_event_sort)
  if(is.na(local_event_tmp2$node[i])){complete=F}
  idx_rm = offspring(tree_tbl,local_event_tmp2$node[i]) %>% filter(node %in% local_event_tmp2$node)
  local_event_tmp2 = local_event_tmp2 %>% filter(!(node %in% idx_rm$node))
  # print(dim(local_event_tmp))
  i = i + 1
}

for(node_in in parents){
  print(node_in)
  tmp = (tree_tbl  %>% offspring(node_in) %>% filter(!grepl("NODE",label)) %>% summarise(node=node_in,total_la = sum(label %in% all_la_strains ), total=n()))
  if(tmp$total_la >=1 & tmp$total_la/tmp$total > 1)
  {
    local_events = rbind(local_events,tmp)
  }
}



local_events = data.frame()
for(node_in in parents){
  print(node_in)
  tmp = (tree_tbl  %>% offspring(node_in) %>% filter(!grepl("NODE",label)) %>% summarise(node=node_in,total_la = sum(label %in% all_la_strains ), total=n()))
  if(tmp$total_la >=1 & tmp$total_la/tmp$total > 0.8)
  {
    local_events = rbind(local_events,tmp)
  }
}


local_event_tmp = local_events %>% arrange(-total)
complete=T
i = 1

gisaid_local %>% group_by(region,pangolin_lineage) %>% summarise(n = n()) %>%
  mutate(freq = n / sum(n))

###



while(complete){
  #if(local_event_sort)
  if(is.na(local_event_tmp$node[i])){complete=F}
  idx_rm = offspring(tree_tbl,local_event_tmp$node[i]) %>% filter(node %in% local_event_tmp$node)
  local_event_tmp = local_event_tmp %>% filter(!(node %in% idx_rm$node))
  # print(dim(local_event_tmp))
  i = i + 1
}

diverged_table %>% offspring(6930) %>% filter()

#### Introductions #### 

### #### 

for(i in 1:nrow(local_event_tmp)){
  tmp_node = local_event_tmp$node[i]
  print(ggtree(tree_subset(tree_all,tmp_node,levels_back =4)) + geom_tiplab() + geom_nodelab(aes(label=node))) 
  if( i > 4){
    break
  }
}


pdf("figures/introductons.pdf")
for(i in 1:nrow(local_event_tmp)){
  tmp_node = local_event_tmp$node[i]
  print(ggtree(tree_subset(tree_all,tmp_node,levels_back =4)) + geom_tiplab() + geom_nodelab(aes(label=node))) 
  if( i > 4){
    break
  }
}
dev.off()