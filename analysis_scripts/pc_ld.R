source("R/utils.R")


gt_matrix = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
gt_matrix_raw = get_gt_matrix(qc_lineages,merged_vcf_filt,qc_lineages$sample_name_fasta)
### Remove all heterozygote sites ###
gt_matrix[gt_matrix == "0/1" | gt_matrix == "0/2"] = "./."
gt_matrix[gt_matrix == "0/0"] = 0
gt_matrix[gt_matrix == "1/1" | gt_matrix == "2/2"] = 1
gt_matrix[gt_matrix == "./."] = NA
gt_matrix = apply(gt_matrix,2,as.numeric)
#princomp(gt_matrix)
#cor(as.numeric(het_common))
relatedness_matrix = cor(gt_matrix,method="spearman",use ="pairwise.complete.obs")

princomp(cor(gt_matrix,method="spearman",use ="pairwise.complete.obs"))
x = (softImpute(gt_matrix,rank.max = 30))
#gt_matrix %>% gt_matrix
xx = complete(gt_matrix,x)

xxxx = cor(t(xx[(apply(xx,1,function(x){sum(x>0.1)}) > 5) & apply(xx,1,sd) !=0,]))
colnames(xxxx) = snp_df_summary$snp_names[apply(xx,1,function(x){sum(x>0.1)}) > 5 & apply(xx,1,sd) !=0]
pheatmap::pheatmap(xxxx,cluster_cols = F,cluster_rows = F)
xxxxx=  cor(t(xx[(apply(gt_matrix,1,function(x){sum(x>0.1,na.rm=T)}) > 5) & apply(xx,1,sd) !=0,]),use = "pairwise.complete.obs")
colnames(xxxxx) = snp_df_summary$snp_names[apply(gt_matrix,1,function(x){sum(x>0.1,na.rm=T)}) > 5 & apply(xx,1,sd) !=0]
pheatmap::pheatmap(xxxxx,cluster_cols = F,cluster_rows = F)

#rownames(xxxx) = 
pc = data.frame(sample_lineage=colnames(gt_matrix),pc1=x$v[,1],pc2=x$v[,2])
pheatmap::pheatmap(xxxx,cluster_cols = F,cluster_rows = F)
pc %>% left_join(qc_lineages,by=c("sample_lineage"="sample_name_fasta")) %>% ggplot(aes(y=pc1,x=pc2,color=nextstrain_clade)) + geom_point() + theme_bw()
pc %>% left_join(qc_lineages,by=c("sample_lineage"="sample_name_fasta")) %>% ggplot(aes(y=pc1,x=pc2,color=lineage)) + geom_point() + theme_bw()