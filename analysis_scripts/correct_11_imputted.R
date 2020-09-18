library(vcfR)
merged_vcf_filt= read.vcfR("~/KRUGLYAK/covid_runs/rerun_all_new_july12//uid_type/outputs/quasi_species/merged/all.vcf")
new_gt =read.csv("data/het_genotypes/het_genotypes_aug1.csv")

gt_tmp = merged_vcf_filt@gt

chr_ref_alt = lapply(colnames(new_gt)[2:length(colnames(new_gt))], function(x){str_split(str_replace_all(x,"X",""),"_")})

pos = as.numeric(unlist(lapply(chr_ref_alt,function(x){x[[1]][1]})))
pair_het = c()
pair_het_after = c()

for( i in 1:nrow(new_gt)){
  sample = new_gt$merged_uid[i]
  print(sample)
  sample_in_vcf = str_replace_all(sample,"_","-")
  idx_col = which(colnames(gt_tmp) == sample_in_vcf)
  tmp_sample_vcf = merged_vcf_filt@gt[,idx_col]
  j = 2
  for(p in pos){
    print(p)
    idx_vcf  = which(as.numeric(merged_vcf_filt@fix[,2])== p)
    print(tmp_sample_vcf[idx_vcf])
 #   s#tr_starts()
    str_out = glue("{gt}:.:.:.:.",gt=str_replace_all(new_gt[i,j]," ",""))
    print(str_out)
    gt_tmp[idx_vcf,idx_col] = str_out
    print(gt_tmp[idx_vcf,idx_col])
    j = j + 1
  }
}

merged_vcf_filt2 = merged_vcf_filt
merged_vcf_filt2@gt = gt_tmp
write.vcf(merged_vcf_filt2,file="data/het_genotypes/final.vcf.gz")
