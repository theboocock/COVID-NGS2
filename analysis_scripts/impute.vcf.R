library(vcfR)

source("R/utils.R")

get_gt_matrix_imputation(impute_vcf){
  
  
}

impute_vcf = function(impute_vcf,qc_input2,lambda=2, rank.max=30,plot=T,threshold=0.8){
  
  
  merged_vcf_filt2= read.vcfR(impute_vcf)
  merged_vcf_filt@gt[,colnames(merged_vcf_filt2@gt) %in% qc_lineages$sample_name_fasta]
  gt_matrix = get_gt_matrix(qc_input2,merged_vcf_filt,sample_names = colnames(merged_vcf_filt2@gt)[2:length( colnames(merged_vcf_filt2@gt))])
  
  gt_matrix[gt_matrix == "0/0"] = -1
  gt_matrix[gt_matrix == "1/1" | gt_matrix == "2/2"] = 1
  gt_matrix[gt_matrix == "./."] = NA
  snp_names = rownames(gt_matrix)
  sample_names = colnames(gt_matrix)
  gt_matrix = apply(gt_matrix,2,as.numeric)
  rownames(gt_matrix) = snp_names
  colnames(gt_matrix) = sample_names
  #princomp(gt_matrix)
  #cor(as.numeric(het_common))
 # relatedness_matrix = cor(gt_matrix,method="spearman",use ="pairwise.complete.obs")
  
#  princomp(cor(gt_matrix,method="spearman",use ="pairwise.complete.obs"))
  soft_impute_matrix = softImpute(gt_matrix,rank.max=30,lambda = 3,trace=T,type="svd")
  gt_matrix_complete = complete(gt_matrix,soft_impute_matrix)
  gt_matrix_complete[gt_matrix_complete > - threshold & gt_matrix_complete < threshold] = NA
  #gt_matrix %>% gt_matrix
  #xx = complete(gt_matrix,x)
  
  
  ## Greater than 10% coverage
  sample_with_lower_coverage = (qc_input2[qc_input2$coverage_3x > .1 & qc_input2$coverage_3x <= 0.8,]$sample_name_fasta)
  
  if (plot)
    pdf("imputation.pdf")
      pheatmap(gt_matrix_complete[,sample_with_lower_coverage],cluster_cols = F, cluster_rows = F,labels_row = "")
      pheatmap(gt_matrix[,sample_with_lower_coverage],cluster_cols = F, cluster_rows = F,labels_row = "")
      
      pheatmap(gt_matrix_complete,cluster_cols = F, cluster_rows = F,labels_row = "")
      pheatmap(gt_matrix,cluster_cols = F, cluster_rows = F,labels_row = "")
      
      #pheatmat(gt_matrix)
    dev.off()
    
  pheatmap::pheatmap(xxxx,cluster_cols = F,cluster_rows = F,labels_row = NULL)
  pheatmap::pheatmap(xxxxx,cluster_cols = F,cluster_rows = F)
  
  
}