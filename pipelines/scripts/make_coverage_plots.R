#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(magrittr)
library(cowplot)
args = commandArgs(trailingOnly=T)
in_coverage= (args[1])
## log and non log ## ??? ###
out_plot = args[2]
sample_name = args[3]
## coverage at different thresholds also one plot per sample

title <- ggdraw() + draw_label(sample_name, size=40)


in_plot = read.table(in_coverage,header=F)
p1 = in_plot%>% ggplot(aes(x=V2,y=V3)) + geom_point() + theme_bw() + xlab("Position") + ylab("Coverage") + theme(text=element_text(size=20))  
p2 = in_plot%>% ggplot(aes(x=V2,y=log2(V3))) + geom_point() + theme_bw() + xlab("Position") + ylab("Log10 coverage") + theme(text=element_text(size=20))
thresholds = seq(5,100,by=1)

cov_df=data.frame(thresholds=thresholds,percent_covered=sapply(thresholds, function(x){sum(in_plot$V3 > x)/nrow(in_plot)}))
p3= cov_df %>% ggplot(aes(y=percent_covered, x= thresholds)) + geom_point() + theme_bw() + xlab("Minimum coverage") + ylab("Proportion of gennome covered") + theme(text=element_text(size=20))
p4=cowplot::plot_grid(p1,p2,p3,ncol=3)
png(out_plot, width=1600,height=600)
cowplot::plot_grid(title,p4,ncol=1,rel_heights=c(.1,1))
dev.off()
