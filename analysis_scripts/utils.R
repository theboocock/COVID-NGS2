library(stringr)
library(vcfR)
library(glue)
get_gt = function(x){
  x = str_split(x,":")
  x2 = unlist(lapply(x,function(x){x[1]}))
  #print(x2)
  return(x2)
}

get_gt_ar= function(x){
  x = str_split(x,":")
  gt = unlist(lapply(x,function(x){x[1]}))
  ar = unlist(lapply(x,function(x){x[5]}))
  dp = unlist(lapply(x,function(x){x[3]}))
  #print(x2)
  return(cbind(gt,ar,dp))
}


get_gt_rbind_df = function(qc_input, vcf_input){
  snp_names = paste(vcf_input@fix[,2],vcf_input@fix[,4],vcf_input@fix[,5],sep="_")
  in_list = list()
  i  = 1
  for (sample in qc_input$sample_name_fasta)
  {
    #print(sample)
    x = get_gt_matrix(qc_input, vcf_input, sample_names = sample)
   # print(x)
    in_list[[i]] = data.frame(snp_name = rownames(x), gt=x[,1],ar=x[,2])
    i = i + 1
  }
  names(in_list) = qc_input$sample_name_fasta
  in_list = bind_rows(in_list, .id="sample_id")
  return(in_list)
}

get_gt_matrix = function(qc_input,vcf_input,sample_names){
  genotypes = vcf_input@gt[,colnames(vcf_input@gt) %in% sample_names]
 # print(genotypes)
  
  if ((length(sample_names))==1){
    genotype_list = cbind(apply(t(genotypes), 2, get_gt_ar ))
    genotype_list = data.frame(gt=genotype_list[1,],ar=genotype_list[2,],dp=genotype_list[3,])
    genotype_list$ar = as.numeric(genotype_list$ar)
    genotype_list$dp = as.numeric(genotype_list$dp)
   # print(genotype_list)
    snp_names = paste(vcf_input@fix[,2],vcf_input@fix[,4],vcf_input@fix[,5],sep="_")
    ## Remove any hets ##
    #genotype_list[genotype_list == "0/1" & genotype_list ==]
    
    rownames(genotype_list) = snp_names
    return(genotype_list)
    
  }else{
    print("HERE")
  genotype_list = cbind(apply(genotypes, 2, get_gt ))
  snp_names = paste(vcf_input@fix[,2],vcf_input@fix[,4],vcf_input@fix[,5],sep="_")
  ## Remove any hets ##
  #genotype_list[genotype_list == "0/1" & genotype_list ==]
  
  rownames(genotype_list) = snp_names
  return(genotype_list)
  }
}

get_concordance_list = function(qc_input, vcf_input){
  qc_input = qc_input %>% group_by(uid) %>% mutate(count_samples=n())
  qc_input = qc_input[qc_input$count_samples > 1,]
  #View(input_file)
  print(qc_input)
  in_list = split(qc_input,f=qc_input$uid)
  #print(in_list)
  uids = names(in_list) 
  snp_names = paste(vcf_input@fix[,2],vcf_input@fix[,4],vcf_input@fix[,5],sep="_")
 # print(snp_names)
  concordance_list = list()
  genotype_lists = list()
  for (i in 1:length(in_list)){
    sample_ids = in_list[[i]]
    sample_name_fasta = sample_ids$sample_name_fasta
    genotypes = vcf_input@gt[,colnames(vcf_input@gt) %in% sample_name_fasta]
    genotype_list = cbind(apply(genotypes, 2, get_gt ))
    #rownames(genotype_list) = snp_names 
   # idx_keep = !(apply(genotype_list, 1,function(x){all(x == "0/0")})) & !(apply(genotype_list, 1,function(x){any(x == "./.")}))
    #genotype_list = genotype_list[!(apply(genotype_list, 1,function(x){all(x == "0/0")})),]
    #genotype_list = genotype_list[!(apply(genotype_list, 1,function(x){any(x == "./.")})),]
    genotype_list = data.frame(genotype_list)
    genotype_list$snp_names = snp_names
    colnames(genotype_list) = c(sample_ids$library_type,"snp_names")
    genotype_lists[[i]] = genotype_list
    print(genotype_list)
    matrix_concordance = matrix(nrow=ncol(genotype_list),ncol=ncol(genotype_list))
    diag(matrix_concordance) = 1
    for(m in 1:(ncol(genotype_list)-1)){
      for (j in 2:ncol(genotype_list)){
        x2 = sum(genotype_list[,m] == genotype_list[,j])/nrow(genotype_list)
        matrix_concordance[j,m] = x2
        matrix_concordance[m,j] = x2
      }
    }
    print(matrix_concordance)
    concordance_list[[i]] = matrix_concordance
  }
  names(concordance_list) = names(in_list)
  names(genotype_lists) = names(in_list)  
  return(list(concordance_list=concordance_list,genotype_lists=genotype_lists))
}

get_ad = function(quasi_in, column){
  
  
  is_vcfr = (class(quasi_in) == "vcfR")[1]
  #
  #class(quasi_in) == "vcfR"
  if (is_vcfr){
    format_string = quasi_in@gt[1,1]
    quasi_in = quasi_in@gt
  }else{
    format_string = quasi_in[1,1]
  }
  idx = which(unlist(lapply(str_split(format_string,":"), function(x) { x  == column})))
  if(length(idx) > 0){
    g = apply(quasi_in[,2:ncol(quasi_in)],1,function(x){ 
      y = sapply(x,function(y) { str_split(y,":")[[1]][idx]});
      y
    })
    g = t(g)
    #alapply(vcf_input@gt[,2:ncol(vcf_input@gt)],1,function(x){ print(x);y = str_split(x,":")[[1]]; y[idx]})
    
  }else{
    message(glue("Could not find column: {column}",column=column))
  }
  return(g)
}


get_concordance_list_ar = function(qc_input, quasi_in){
  qc_input = qc_input %>% group_by(uid) %>% mutate(count_samples=n())
  qc_input = qc_input[qc_input$count_samples > 1,]
  #View(input_file)
  print(qc_input)
  in_list = split(qc_input,f=qc_input$uid)
  #print(in_list)
  uids = names(in_list) 
  snp_names = paste(quasi_in@fix[,2],quasi_in@fix[,4],quasi_in@fix[,5],sep="_")
  concordance_list = list()
  genotype_lists = list()
  for (i in 1:length(in_list)){
    sample_ids = in_list[[i]]
    sample_name_fasta = sample_ids$sample_name_fasta
    sample_types = sample_ids$sample_type
    sample_types = str_replace_all(sample_types,"neb_old","neb")
    genotypes = quasi_in@gt[,c(1,which(colnames(quasi_in@gt) %in% sample_name_fasta))]
    genotype_list = cbind(apply(genotypes, 2, get_gt ))
    ar_list = get_ad(genotypes,"AR")
    dp_list = get_ad(genotypes,"DP")
    #genotype_list = genotype_list[!(apply(genotype_list, 1,function(x){all(x == "0/0")})),]
    ##keep_idx = !(apply(genotype_list, 1,function(x){any(x == ".")}))
    #genotype_list = genotype_list
    print(ar_list)
    #rownames(genotype_list) = snp_names[keep_idx]
    genotype_list = data.frame(genotype_list)
    colnames(genotype_list) = sample_types
    genotype_list$snp_name = snp_names
    #print(ar_list)
    ar_colnames = paste("ar",sample_ids$library_type,sep="_")
    ar_list =  apply(ar_list,2,as.numeric)
    dp_colnames = paste("dp",sample_ids$library_type,sep="_")
    dp_list = apply(dp_list,2,as.numeric)
    genotype_list =(cbind(genotype_list,ar_list,dp_list))
    colnames(genotype_list) = c("format",sample_ids$library_type,"snp_names",ar_colnames,dp_colnames)
    
      #genotype_list[,4] = ar_list
    
    genotype_lists[[i]] = genotype_list
    #print(genotype_list)
    matrix_concordance = matrix(nrow=ncol(genotype_list),ncol=ncol(genotype_list))
    diag(matrix_concordance) = 1
    for(m in 1:(ncol(genotype_list)-1)){
      for (j in 2:ncol(genotype_list)){
        x2 = sum(genotype_list[,m] == genotype_list[,j])/nrow(genotype_list)
        matrix_concordance[j,m] = x2
        matrix_concordance[m,j] = x2
      }
    }
    print(matrix_concordance)
    concordance_list[[i]] = matrix_concordance
  }
  names(concordance_list) = names(in_list)
  names(genotype_lists) = names(in_list)  
  return(list(concordance_list=concordance_list,genotype_lists=genotype_lists))
}

calculate_af= function(vcf_input){
  gt_tmp = vcf_input@gt
  print(gttmp)
}