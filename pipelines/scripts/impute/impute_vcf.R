#!/usr/bin/env Rscript

get_gt_matrix = function(qc_input,vcf_input,sample_names){
  genotypes = vcf_input@gt[,colnames(vcf_input@gt) %in% sample_names]
 # print(genotypes)
  
  if ((length(sample_names))==1){
    genotype_list = cbind(apply(t(genotypes), 2, get_gt_ar ))
    genotype_list = data.frame(gt=genotype_list[1,],ar=genotype_list[2,])
    print(genotype_list)
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

get_gt = function(x){
  x = str_split(x,":")
  x2 = unlist(lapply(x,function(x){x[1]}))
  #print(x2)
  return(x2)
}

get_gt_ar= function(x){
  x = str_split(x,":")
  x2 = unlist(lapply(x,function(x){x[1]}))
  x3 = unlist(lapply(x,function(x){x[5]}))
  #print(x2)
  return(cbind(x2,x3))
}


args = commandArgs(trailingOnly=T)

library(vcfR)
library(stringr)
library(softImpute)
vcf_master = args[1]
vcf_output = args[2]
qc_input = args[3]

qc_input = read.delim(qc_input, sep="\t",header=T)

vcf_master = read.vcfR(vcf_master)


#### ####

#### For each sample only impute using the other sample data #####

impute_each_sample = function(sample_focal, gt_matrix,lambda=3,rank.max=3,type="svd",threshold=0.7){
    ### Remove all sample focals.
    sample_names = colnames(gt_matrix)
	keep_idx_in = grep("-",sample_names)
    idx_focal = which(sample_names == sample_focal) 
    reworked_sample_names = str_replace_all(sample_names, "-","_")
    sample_focal_pre= str_split(sample_focal,"_")[[1]][2]
    ### Remove all reference to the focal sample ### 
    reworked_sample_names_split = unlist(lapply(str_split(reworked_sample_names, "_"),function(x){x[2]}))
    keep_idx = which(!(reworked_sample_names_split == sample_focal_pre))
    keep_idx = keep_idx[keep_idx%in% keep_idx_in]
    idx_keep_and_sample_focal = c(keep_idx,idx_focal)
    impute = gt_matrix[,idx_keep_and_sample_focal]
    xxx = softImpute(impute, rank.max=rank.max,lambda=lambda,trace=F,type=type)
    xx = complete(gt_matrix[,idx_keep_and_sample_focal], xxx)
    colnames(xx) = colnames(gt_matrix)[idx_keep_and_sample_focal]
    rownames(xx) = rownames(gt_matrix)
    return(xx[,ncol(xx)])
}


sample_names = colnames(vcf_master@gt)[2:ncol(vcf_master@gt)]
gt_matrix = get_gt_matrix(qc_input, vcf_master,sample_names=sample_names)
gt_matrix[gt_matrix == "0/0"] = -1
#### Convert this to two rows one for each of the variants 

gt_matrix[gt_matrix == "1/1"] = 1
gt_matrix[gt_matrix == "2/2"] = 3
gt_matrix[gt_matrix == "./."] = NA
gt_matrix = apply(gt_matrix,2,as.numeric)

### Focal samples ###
sample_matrix_focal  = matrix(nrow=nrow(gt_matrix), ncol=ncol(gt_matrix))

sample_names_focused = sample_names[!(sample_names %in% qc_input$sample_name_fasta)]
i=1

library(parallel)
print(sample_names_focused)
rows = mclapply(sample_names_focused, function(x){

        sample_gt_matrix = impute_each_sample(x, gt_matrix)
        return(sample_gt_matrix)
},mc.cores=16)

sample_matrix_focal = gt_matrix
for (i in 1:length(rows)){
    #### Ok figure out the column ##
    sample_id = sample_names_focused[i]
    idx_replace = which(colnames(sample_matrix_focal) == sample_id)
    print(idx_replace)
    sample_matrix_focal[,idx_replace] = rows[[i]]
}
#for (sample in sample_names){
#    if (!(sample %in% qc_input$sample_name_fasta)){
#        sample_gt_matrix = impute_each_sample(sample, gt_matrix)
#        sample_matrix_focal[,i] = sample_gt_matrix
#        i=i+1
#    }
#}
write.table(sample_matrix_focal, file="outputs/impute/merged/impute_new_design.txt", quote=F,row.names=T,col.names=T)

xxx = softImpute(gt_matrix,rank.max=30,lambda = 3,trace=T,type="svd")
xx = complete(gt_matrix,xxx)


colnames(gt_matrix) = sample_names
rownames(gt_matrix) = paste(vcf_master@fix[,2],vcf_master@fix[,4],vcf_master@fix[,5],sep="_")
colnames(xx) = sample_names
rownames(xx) = paste(vcf_master@fix[,2],vcf_master@fix[,4],vcf_master@fix[,5],sep="_")   
write.table(gt_matrix, file="outputs/impute/merged/impute_orig.txt", quote=F,row.names=T,col.names=T)

#### create new vcfoutput from xx ######

create_gt_matrix = function(gt_matrix,vcf_in, threshold=0.8){
   fixed_region = ":.:.:.:."
   genotypes = apply(gt_matrix, 2, function(x){
        y= as.character(x)
        y[x > threshold] = paste0("1/1",fixed_region)
        y[x > 2] = paste0("2/2",fixed_region)
        y[x < -threshold] = paste0("0/0",fixed_region)
        y[x>-threshold & x < threshold] = paste0("./.", fixed_region)
        return(y)
   })
  y =  cbind(vcf_in@gt[,1], genotypes) 
  #print(y)
  return(y)
}


write.table(xx,file="outputs/impute/merged/impute.txt", quote=F,row.names=T,col.names=T)

#gt_matrix2 = get_gt_matrix(qc_input, vcf_downsample, sample_names)  
#gt_matrix2= get_gt_matrix(qc_input, vcf_master,sample_names=qc_input$sample_name_fasta)
#gt_matrix2[gt_matrix2 == "0/0"] = -1
#gt_matrix2[gt_matrix2 == "1/1" | gt_matrix2 == "2/2"] = 1
#gt_matrix2[gt_matrix2 == "./."] = NA
##gt_matrix2 = apply(gt_matrix2,2,as.numeric)

gtm = create_gt_matrix(sample_matrix_focal,vcf_master, threshold=0.7)
vcf_master@gt = gtm
write.vcf(vcf_master, file=vcf_output)

