source("R/utils.R")
gt_matrix = get_gt_matrix(qc_input2,merged_vcf_filt2,sample_names = colnames(merged_vcf_filt2@gt)[2:length( colnames(merged_vcf_filt2@gt))])

#### ####
#gt_matrix = gt_matrix[,qc_input2$sample_name_fasta[qc_input2$coverage_3x > .2]]
gt_matrix[gt_matrix == "0/0"] = -1
gt_matrix[gt_matrix == "1/1" | gt_matrix == "2/2"] = 1
gt_matrix[gt_matrix == "./."] = NA
gt_matrix = apply(gt_matrix,2,as.numeric)

xxx = softImpute(gt_matrix,rank.max=30,lambda = 3,trace=T,type="svd")
xx = complete(gt_matrix,xxx)

baddies = qc_input2[qc_input2$coverage_3x > 0.1 & qc_input2$coverage_3x < 0.8,]$sample_name_fasta
#gt_matrix[is.na(gt_matrix)] = as.double("NA")
### ####
xx[xx > - .8 & xx < .8] = NA
pheatmap(gt_matrix,cluster_rows = F,cluster_cols = F)
pheatmap(xx,cluster_rows = F,cluster_cols = F)


plot(apply(xx,2,function(x){sum(is.na(x)/length(x))}) )
pheatmap(xx[,qc_input2[qc_input2$coverage_3x < 0.8,]$sample_name_fasta],cluster_rows = F,cluster_cols = F)
pheatmap(gt_matrix[,qc_input2[qc_input2$coverage_3x < 0.8,]$sample_name_fasta],cluster_rows = F,cluster_cols = F)
plot(apply(xx,2,function(x){sum(is.na(x)/length(x))}),  apply(gt_matrix,2,function(x){sum(is.na(x)/length(x))}),xlab="Coverage after imputation",ylab="Coverage before imputation")
summary_imputation= data.frame(sample=colnames(xx),after=apply(xx,2,function(x){sum(!is.na(x))/length(x)}), before=apply(gt_matrix,2,function(x){sum(!is.na(x))/length(x)}))
ss = summary_imputation %>% filter(before < 0.8 & after > 0.8) %>% left_join(qc_input2, by=c("sample"="sample_name_fasta"))


test_all = read.table("~/KRUGLYAK/covid_runs/rerun_all_new_july20/downsample/outputs/impute/merged/impute_new_design.txt")
downsample_all = read.table("~/KRUGLYAK/covid_runs/rerun_all_new_july20/downsample/outputs/impute/merged/impute.txt")
downsample_all[downsample_all > -.7 & downsample_all < 0.7] = NA
downsample_orig = read.table("~/KRUGLYAK/covid_runs/rerun_all_new_july20/downsample/outputs/impute/merged/impute_orig.txt")
xx = (test_all[!is.na(test_all[,grep("LA_0065",colnames(downsample_all))][,9]),grep("LA_0065",colnames(downsample_all))])
xxxx = downsample_orig[!is.na(downsample_orig[,grep("LA_0003",colnames(downsample_all))][,9]),grep("LA_0003",colnames(downsample_all))]
#xxxx = downsample_orig[,grep("LA_0021",colnames(downsample_all))]
cov = (apply(downsample_all,2,function(x){sum(!is.na(x))/length(x)}))
impute =  (apply(downsample_orig,2,function(x){sum(!is.na(x))/length(x)}))
#xx = (downsample_all[,grep("LA_0021",colnames(downsample_all))])
sum((xx[,3] * xx[,9] > 0)/nrow(xx),na.rm=T )
sapply(1:9, function(x){((sum(xx[,x] * xx[,9] > 0,na.rm=T)/sum(!is.na(xx[,9]))))})
sapply(1:9, function(x){((sum(xxxx[,x] * xxxx[,9] > 0,na.rm=T)/sum(!is.na(xxxx[,9]))))})
# Accuracy
sapply(1:9, function(x){((sum(xx[,x] * xx[,9] > 0,na.rm=T)/nrow(xx)))})
sapply(1:9, function(x){((sum(xxxx[,x] * xxxx[,9] > 0,na.rm=T)/nrow(xxxx)))})


vcf = read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july20/downsample/outputs/impute/merged/merged.vcf.gz")

gt_matrix = get_gt_matrix(qc_lineages, vcf,sample_names = colnames(vcf@gt)[2:length(colnames(vcf@gt))])
#### Convert this to two rows one for each of the variants

impute_each_sample = function(sample_focal, gt_matrix,lambda=3,rank.max=3,type="svd",threshold=0.7){
  ### Remove all sample focals.
  sample_names = colnames(gt_matrix)
  idx_focal = which(sample_names == sample_focal)
  print(idx_focal)
  keep_idx_in = grep("-",sample_names)
  reworked_sample_names = str_replace_all(sample_names, "-","_")
  sample_focal_pre= str_split(sample_focal,"_")[[1]][2]
  ### Remove all reference to the focal sample ###
  reworked_sample_names_split = unlist(lapply(str_split(reworked_sample_names, "_"),function(x){x[2]}))
  keep_idx = which(!(reworked_sample_names_split == sample_focal_pre))
  ###
  keep_idx = keep_idx[keep_idx%in% keep_idx_in]
  idx_keep_and_sample_focal = c(keep_idx,idx_focal)
  impute = gt_matrix[,idx_keep_and_sample_focal]
  print(dim(impute))
  xxx = softImpute(impute, rank.max=rank.max,lambda=lambda,trace=F,type=type)
  xx = complete(gt_matrix[,idx_keep_and_sample_focal], xxx)
  colnames(xx) = colnames(gt_matrix)[idx_keep_and_sample_focal]
  rownames(xx) = rownames(gt_matrix)
  return(xx[,ncol(xx)])
}
gt_matrix[gt_matrix == "0/0"] = -1
gt_matrix[gt_matrix == "1/1"] = 1
gt_matrix[gt_matrix == "2/2"] = 3
gt_matrix[gt_matrix == "./."] = NA
gt_matrix = apply(gt_matrix,2,as.numeric)
sample = "LA_0065_3.0"
idx_sample = which("LA-0065" == colnames(gt_matrix))
xx = impute_each_sample(sample,gt_matrix)

threshold = 0.7

for(lambda in seq(1,20,by=1)){
  xx = impute_each_sample(sample,gt_matrix,lambda = lambda,threshold = .9)
  
  xx[xx<= -threshold]= -1
   xx[xx > threshold] = 1
   xx[xx >= -threshold & xx <= threshold] = NA
   multi = xx * gt_matrix[,idx_sample]
  print(sum(multi,na.rm=T)/sum(!is.na(multi)))
  #sum(xx * gt_matrix[,3] > 0,na.rm=T)/
}

