setwd("C:/Users/gael/charbel_2022/rap_h9/")

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

# .libPaths(c("/shared/projects/dubii2021/gcastel/r_packages_gcastel/","/shared/ifbstor1/software/miniconda/envs/r-4.0.3/lib/R/library/"))

library(dendextend)
library(reticulate)
library(dplyr)
library(tidyr)
library(tibble)
library(reticulate)
library(pcaMethods)
library(devtools)
library(ggplot2)
library(grid)
library(Matrix)
library(stringr)
library("Biobase")
library("BiocGenerics")
library(uwot)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(DESeq2)
library(gplots)
library(DiffBind)
library(parallel)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggpmisc)

timeNow<-function(x=NULL, y=NULL, m=NULL, d=NULL , h=NULL, min=NULL,ymd_h_m=NULL){x=date();
date=unlist(str_split(pattern = "-", Sys.Date()))
time=unlist(str_split(pattern = ":",unlist(str_split(pattern=" ",Sys.time()))[2]))
return(paste(c(date,time), collapse = "_"))
}


chrX_annot=read.table("chrX_annot.gtf", header=F, sep="\t")
chrX_annot_gene=chrX_annot[chrX_annot$V3=='gene',]

bin_tab=read.table("../data/rap_h9_no_igg_norm_scores_per_bin.tab", header=T, sep="\t", comment.char = '')
colnames(bin_tab)=c("chr","start","end","Input_Primed_Rep1","Input_Primed_Rep2","Input_PXGL_Rep1",
                    "Input_PXGL_Rep2","RAP_XIST_Primed_Rep1","RAP_XIST_Primed_Rep2","RAP_XIST_PXGL_Rep1",
                    "RAP_XIST_PXGL_Rep2")  

bin_tab=bin_tab[bin_tab$chr=='chrX',]


rownames(bin_tab)=paste(sep = "_",bin_tab$chr,bin_tab$start,bin_tab$end)
bin_tab$chr=NULL
bin_tab$start=NULL
bin_tab$end=NULL

# removal of Input, same rationale as for IgG 'forced' sequencing, as discussed with Magali.
bin_tab=bin_tab[,!grepl('Input',cn(bin_tab))]

bin_tab=bin_tab[complete.cases(bin_tab),]

levels(as.factor(sub("_.*","",rn(bin_tab))))


#split rap_h9_bw_summary_matrix_charbel into autosomes and X
#!!!! ask Charbel: 3353 bins annotated as chrY in female H9 !!!!! probably X paralogs ? I remove them but maybe include as X ?
rap_h9_bw_summary_matrix_charbel=bin_tab

rap_h9_bw_summary_matrix_charbel_x=rap_h9_bw_summary_matrix_charbel[grep("chrX",rn(rap_h9_bw_summary_matrix_charbel)),]
rap_h9_bw_summary_matrix_charbel_autosomes=rap_h9_bw_summary_matrix_charbel[-c(grep("chrX",rn(rap_h9_bw_summary_matrix_charbel)),grep("chrY",rn(rap_h9_bw_summary_matrix_charbel)),grep("chrUn",rn(rap_h9_bw_summary_matrix_charbel)),grep("chrM",rn(rap_h9_bw_summary_matrix_charbel))),]

# need to sum every 100 rows to have coverage score for 1 Mb bins


rap_h9_bw_summary_matrix_charbel_x %>% 
  group_by(grp = rep(row_number(), length.out = n(), each = 100)) %>%
  summarise(across(everything(), list(sum))) -> rap_h9_bw_summary_matrix_charbel_1mb_bins_x

bin_gene_density=NULL

for (i in seq(0,156040895,1e6)){
  
  name=as.character(i)

  bin_gene_density[[name]]=sapply(chrX_annot_gene$V4,function(x){ifelse(as.numeric(x)>i & as.numeric(x)<i+1e6,T,F)})

}

bin_gene_density_df=do.call("cbind",bin_gene_density)
bin_gene_density=colSums(bin_gene_density_df)

rap_h9_bw_summary_matrix_charbel_1mb_bins_x$grp=NULL

rap_h9_plot=NULL
rap_h9_plot$sample=c(rep(cn(rap_h9_bw_summary_matrix_charbel_1mb_bins_x),each=nrow(rap_h9_bw_summary_matrix_charbel_1mb_bins_x)))
rap_h9_plot$bin_gene_density=rep(bin_gene_density,ncol(rap_h9_bw_summary_matrix_charbel_1mb_bins_x))
rap_h9_plot$xist_coverage_score=unlist(c(rap_h9_bw_summary_matrix_charbel_1mb_bins_x[,1:ncol(rap_h9_bw_summary_matrix_charbel_1mb_bins_x)]))
rap_h9_plot=as.data.frame(rap_h9_plot)

# pdf(file="rap_h9_XIST_coverage_function_of_gene_density.pdf")
ggplot(data = rap_h9_plot, aes(x=bin_gene_density,y=xist_coverage_score)) +
  geom_point() +
  facet_wrap(~sample) +
  stat_poly_line() +
  stat_poly_eq()
dev.off()

## transposable elements
transposable_elements=read.table("rmsk.txt", header=F, sep="\t")
transposable_elements=transposable_elements[,c(6:8,10:13)]
colnames(transposable_elements)=c("chr","start","end","strand","TE_name","TE_family","TE_subfamily")

transposable_elements_x=transposable_elements[transposable_elements$chr=='chrX',]

te_tab=read.table("rap_h9_bw_summary_transposable_elements_x.tab", header=T, sep="\t", comment.char = '')
colnames(te_tab)=sub(".*X\\.","",sub(".hg38.BPM.all.bw.*","",cn(te_tab))) #probably additional fields for TE name, TE family...

te_tab_x=te_tab[te_tab$chr=='chrX',]


dim(transposable_elements_x)==dim(te_tab_x) #TRUE
tail(transposable_elements_x$start==te_tab_x$start) #TRUE

rap_te_df=merge(transposable_elements_x,te_tab_x, by=c('chr','start','end'))

rap_te_df %>%
  group_by(TE_family) %>%
  group_split() -> rap_te_family

rap_te_family_df_plot=NULL

for (i in seq(length(rap_te_family))){
  chr=rep(rap_te_family[[i]]$chr,4) # 4 samples
  start=rep(rap_te_family[[i]]$start,4)
  end=rep(rap_te_family[[i]]$end,4)
  strand=rep(rap_te_family[[i]]$strand,4)
  TE_name=rep(rap_te_family[[i]]$TE_name,4)
  TE_family=rep(rap_te_family[[i]]$TE_family,4)
  TE_subfamily=rep(rap_te_family[[i]]$TE_subfamily,4)
  sample=rep(c('PXGL_Rep1', 'PXGL_Rep2', 'Primed_Rep1', 'Primed_Rep2'),each=nrow(rap_te_family[[i]]))
  value=c(rap_te_family[[i]]$RAP_XIST_PXGL_Rep1,rap_te_family[[i]]$RAP_XIST_PXGL_Rep2,rap_te_family[[i]]$RAP_XIST_Primed_Rep1,rap_te_family[[i]]$RAP_XIST_Primed_Rep2)
  rap_te_family_df_plot[[i]]=data.frame(chr=chr,start=start,end=end,TE_name=TE_name,TE_family=TE_family,TE_subfamily=TE_subfamily,value=value,sample=sample)
}

# for (i in seq(length(rap_te_family_df_plot))){
#   
#   rap_boxplot=ggplot(rap_te_family_df_plot[[i]], aes(x=sample, y=value, fill=sample)) +
#     geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
#     geom_point(size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
#     facet_wrap(~TE_subfamily, scales="free", nrow=6, ncol=5)+
#     theme_bw()+
#     theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
#     labs(title=paste0("TE_family: ",sub("\\?.*","",unique(rap_te_family_df_plot[[i]]$TE_family))))+
#     scale_x_discrete(expand = c(0, 0.5)) +
#     ylab("XIST Coverage score") +
#     guides(fill=guide_legend(ncol=1)) +
#     theme(axis.line.x = element_blank(),
#           axis.line.y = element_blank(),
#           axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           panel.grid.major.x = element_blank(),
#           panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
#           panel.grid.minor.x = element_blank(),
#           panel.grid.minor.y = element_blank())+
#     stat_boxplot(geom='errorbar', linetype=1, width=0.8, size=0.25)+
#     stat_compare_means(vjust=1,label="p.signif",size=3,ref.group = ".all.", method="wilcox.test", paired=T)  
#   
#   
#   pdf(paste(sep="_",timeNow(),"rap_h9_boxplots",sub("\\?.*","",unique(rap_te_family_df_plot[[i]]$TE_family)),"transposable_elements.pdf"),width = 11.7,height = 12)
#   print(rap_boxplot)
#   dev.off()
#   
# }

#################################

rownames(te_tab)=paste(sep = "_",te_tab$chr,te_tab$start,te_tab$end) #combine with TE_name
te_tab$chr=NULL
te_tab$start=NULL
te_tab$end=NULL

# boxplots

rap_boxplot=ggplot(dataDeSeqForPlot, aes(x=te_family, y=value, fill=cellType)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  geom_point(size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~cell_type, scales="free", nrow=6, ncol=5)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  scale_x_discrete(expand = c(0, 0.5)) +
  ylab("XIST Coverage score") +
  xlab("TE family")+
  scale_fill_manual("TE family :", values = colors)+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  stat_boxplot(geom='errorbar', linetype=1, width=0.8, size=0.25)+
  scale_y_continuous(trans="pseudo_log",breaks =c(0,5,10,50,100), limits=c(0,max(dataDeSeqForPlot$value)+100))+
  stat_compare_means(vjust=1,label="p.signif",size=3,ref.group = "tsc", method="wilcox.test", paired=FALSE)  


pdf(paste(sep="_",timeNow(),"rap_h9_boxplots_transposable_elements.pdf"),width = 11.7,height = 12)
rap_boxplot
dev.off()

# do it on downsampled bam, for 'normalisation'/quantitative analysis

# determine min tot reads bam from stats, take it for downsampling bam
bam_stat_files=list.files(path = 'stats/', pattern = 'stat.txt')
names(bam_stat_files)=bam_stat_files

bam_reads_tot=NULL

for(file in bam_stat_files){

  
  bam_stat=read.table(paste0('stats/',file), header=F, sep="\t")
  colnames(bam_stat)=c('chr','length','reads_mapped','reads_unmapped')
  tot=sum(as.numeric(bam_stat$reads_mapped))
  bam_reads_tot[[file]]=tot
  
}

bam_reads_tot=unlist(bam_reads_tot)

min(bam_reads_tot)

# check bam_stat_downsampled
bam_stat_files_downsampled=list.files(path = 'stats/', pattern = 'stat_downsampled.txt')
names(bam_stat_files_downsampled)=bam_stat_files_downsampled

bam_reads_tot_downsampled=NULL

for(file in bam_stat_files_downsampled){
  
  
  bam_stat_downsampled=read.table(paste0('stats/',file), header=F, sep="\t")
  colnames(bam_stat_downsampled)=c('chr','length','reads_mapped','reads_unmapped')
  tot=sum(as.numeric(bam_stat_downsampled$reads_mapped))
  bam_reads_tot_downsampled[[file]]=tot
  
}

bam_reads_tot_downsampled=unlist(bam_reads_tot_downsampled)
bam_reads_tot_downsampled

################same analysis from downsampled bam -> bw

sample_annot=read.csv("SampleSheet.csv", header=T, comment.char = '')

te_tab_downsampled=read.table("rap_h9_bw_summary_transposable_elements_downsampled_chrX.tab", header=T, sep="\t", comment.char = '')
colnames(te_tab_downsampled)=c("chr","start","end",sample_annot$Sample_Name) #ok col tab same order as sample_annot

te_tab_downsampled_x=te_tab_downsampled[te_tab_downsampled$chr=='chrX',]

# removal of Input, same rationale as for IgG 'forced' sequencing, as discussed with Magali.
te_tab_downsampled_x=te_tab_downsampled_x[,!grepl('Input',cn(te_tab_downsampled_x))]


dim(transposable_elements_x)[1]==dim(te_tab_downsampled_x)[1] #TRUE 
tail(transposable_elements_x$start==te_tab_downsampled_x$start) #TRUE

rap_te_downsampled_df=merge(transposable_elements_x,te_tab_downsampled_x, by=c('chr','start','end'))

rap_te_downsampled_df %>%
  group_by(TE_family) %>%
  group_split() -> rap_te_downsampled_family

rap_te_downsampled_family_df_plot=NULL

for (i in seq(length(rap_te_downsampled_family))){
  chr=rep(rap_te_downsampled_family[[i]]$chr,4) # 4 samples
  start=rep(rap_te_downsampled_family[[i]]$start,4)
  end=rep(rap_te_downsampled_family[[i]]$end,4)
  strand=rep(rap_te_downsampled_family[[i]]$strand,4)
  TE_name=rep(rap_te_downsampled_family[[i]]$TE_name,4)
  TE_family=rep(rap_te_downsampled_family[[i]]$TE_family,4)
  TE_subfamily=rep(rap_te_downsampled_family[[i]]$TE_subfamily,4)
  sample=rep(c('PXGL_Rep1', 'PXGL_Rep2', 'Primed_Rep1', 'Primed_Rep2'),each=nrow(rap_te_downsampled_family[[i]]))
  value=c(rap_te_downsampled_family[[i]]$RAP_XIST_PXGL_Rep1,rap_te_downsampled_family[[i]]$RAP_XIST_PXGL_Rep2,rap_te_downsampled_family[[i]]$RAP_XIST_Primed_Rep1,rap_te_downsampled_family[[i]]$RAP_XIST_Primed_Rep2)
  rap_te_downsampled_family_df_plot[[i]]=data.frame(chr=chr,start=start,end=end,TE_name=TE_name,TE_family=TE_family,TE_subfamily=TE_subfamily,value=value,sample=sample)
}

for (i in seq(length(rap_te_downsampled_family_df_plot))){

  rap_boxplot=ggplot(rap_te_downsampled_family_df_plot[[i]], aes(x=sample, y=value, fill=sample)) +
    geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
    geom_point(size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
    facet_wrap(~TE_subfamily, scales="free", nrow=6, ncol=5)+
    theme_bw()+
    theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
    labs(title=paste0("TE_family: ",sub("\\?.*","",unique(rap_te_downsampled_family_df_plot[[i]]$TE_family))))+
    scale_x_discrete(expand = c(0, 0.5)) +
    ylab("XIST Coverage score") +
    guides(fill=guide_legend(ncol=1)) +
    theme(axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())+
    stat_boxplot(geom='errorbar', linetype=1, width=0.8, size=0.25)+
    stat_compare_means(vjust=1,label="p.signif",size=3,ref.group = ".all.", method="wilcox.test", paired=T)


  pdf(paste(sep="_",timeNow(),"rap_h9_boxplots",sub("\\?.*","",unique(rap_te_downsampled_family_df_plot[[i]]$TE_family)),"transposable_elements_downsampled.pdf"),width = 11.7,height = 12)
  print(rap_boxplot)
  dev.off()

}

## correlation TE_density downsampled
bin_tab_downsampled=read.table("rap_h9_no_igg_norm_scores_per_bin_downsampled.tab", header=T, sep="\t", comment.char = '')
colnames(bin_tab_downsampled)=c("chr","start","end",sample_annot$Sample_Name) #ok col tab same order as sample_annot


bin_tab_downsampled_x=bin_tab_downsampled[bin_tab_downsampled$chr=='chrX',]

rownames(bin_tab_downsampled_x)=paste(sep = "_",bin_tab_downsampled_x$chr,bin_tab_downsampled_x$start,bin_tab_downsampled_x$end)
bin_tab_downsampled_x$chr=NULL
bin_tab_downsampled_x$start=NULL
bin_tab_downsampled_x$end=NULL

# removal of Input, same rationale as for IgG 'forced' sequencing, as discussed with Magali.
bin_tab_downsampled_x=bin_tab_downsampled_x[,!grepl('Input',cn(bin_tab_downsampled_x))]

# bin_tab_downsampled_x=bin_tab_downsampled_x[complete.cases(bin_tab_downsampled_x),]

levels(as.factor(sub("_.*","",rn(bin_tab_downsampled_x))))

rap_h9_downsampled_bw_summary_matrix_charbel_x=bin_tab_downsampled_x
rap_h9_downsampled_bw_summary_matrix_charbel_autosomes=bin_tab_downsampled[-c(grep("chrX",rn(rap_h9_downsampled_bw_summary_matrix_charbel)),grep("chrY",rn(rap_h9_downsampled_bw_summary_matrix_charbel)),grep("chrUn",rn(rap_h9_downsampled_bw_summary_matrix_charbel)),grep("chrM",rn(rap_h9_downsampled_bw_summary_matrix_charbel))),]

# need to sum every 100 rows to have coverage score for 1 Mb bins


rap_h9_downsampled_bw_summary_matrix_charbel_x %>% 
  group_by(grp = rep(row_number(), length.out = n(), each = 100)) %>%
  summarise(across(everything(), list(sum))) -> rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x

bin_te_density=NULL

for (i in seq(0,156040895,1e6)){
  
  name=as.character(i)
  
  bin_te_density[[name]]=sapply(transposable_elements_x$start,function(x){ifelse(as.numeric(x)>i & as.numeric(x)<i+1e6,T,F)})
  
}

bin_te_density_df=do.call("cbind",bin_te_density)
bin_te_density=colSums(bin_te_density_df)

rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x$grp=NULL

rap_h9_downsampled_plot=NULL
rap_h9_downsampled_plot$sample=c(rep(cn(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x),each=nrow(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x)))
rap_h9_downsampled_plot$bin_te_density=rep(bin_te_density,ncol(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x))
rap_h9_downsampled_plot$xist_coverage_score=unlist(c(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x[,1:ncol(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x)]))
rap_h9_downsampled_plot=as.data.frame(rap_h9_downsampled_plot)

hist(rap_h9_downsampled_plot$bin_te_density)

pdf(file="rap_h9_downsampled_XIST_coverage_function_of_te_density.pdf")
ggplot(data = rap_h9_downsampled_plot, aes(x=bin_te_density,y=xist_coverage_score)) +
  geom_point() +
  facet_wrap(~sample) +
  stat_poly_line() +
  stat_poly_eq()
dev.off()


#################te_density_corr_plot_per_family
transposable_elements_x %>%
  group_by(TE_family) %>%
  group_split() -> transposable_elements_x_family

for(j in seq(length(transposable_elements_x_family))){
  bin_te_density=NULL

  for (i in seq(0,156040895,1e6)){
      
      name=as.character(i)
      
      bin_te_density[[name]]=sapply(transposable_elements_x_family[[j]]$start,function(x){ifelse(as.numeric(x)>i & as.numeric(x)<i+1e6,T,F)})
      
    }
    
  bin_te_density_df=do.call("cbind",bin_te_density)
  bin_te_density=colSums(bin_te_density_df)
  
  rap_h9_downsampled_plot=NULL
  rap_h9_downsampled_plot$sample=c(rep(cn(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x),each=nrow(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x)))
  rap_h9_downsampled_plot$bin_te_density=rep(bin_te_density,ncol(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x))
  rap_h9_downsampled_plot$xist_coverage_score=unlist(c(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x[,1:ncol(rap_h9_downsampled_bw_summary_matrix_charbel_1mb_bins_x)]))
  rap_h9_downsampled_plot=as.data.frame(rap_h9_downsampled_plot)
  
  hist(rap_h9_downsampled_plot$bin_te_density)
  
  plot=ggplot(data = rap_h9_downsampled_plot, aes(x=bin_te_density,y=xist_coverage_score)) +
    geom_point() +
    facet_wrap(~sample) +
    stat_poly_line() +
    stat_poly_eq()
  pdf(file=paste0("rap_h9_downsampled_XIST_coverage_function_of_te_density_",sub("\\?.*","",unique(transposable_elements_x_family[[j]]$TE_family)),".pdf"))
  print(plot)
  dev.off()
}
