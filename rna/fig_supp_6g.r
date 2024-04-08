setwd('C:/Users/surfi/Desktop/charbel_paper_2024/2024_04/charbel_paper_github_v1/rna/')

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))


# BiocManager::install('EnhancedVolcano')

library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(Matrix)
library(stringr)
library(BiocManager)
library(devtools)
library(pcaMethods)
library(data.table)
library(parallel)
library(DESeq2)
library(ggpubr)
library(tidyverse)
library(EnhancedVolcano)

gene_chr=read.table("../../hg38_gene_chr_list.txt")
gene_chr$V2=NULL
colnames(gene_chr)=c("gene","chr")

raw_count_files=list.files(path = 'raw_counts/',pattern = "D1507.*_raw_counts.txt")
sample_annot=read.csv("sample_annot_all.txt", header=T, comment.char = '#')
sample_annot=sample_annot[sample_annot$dataset=="rna_seq_d1507",]
rownames(sample_annot)=sample_annot$Sample_ID
sample_annot$sample_id=sample_annot$Sample_ID
sample_annot$group=sub('_[^_]*$', '', sample_annot$Sample_Name)
sample_annot$group=str_replace_all(sample_annot$group,"CTL","1")
sample_annot$group=str_replace_all(sample_annot$group,"Mix","2")
sample_annot$replicate=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$replicate[sample_annot$replicate=="primed"]<-"Rep1"
sample_annot$cell_type=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$cell_type[sample_annot$cell_type!="primed"]<-"pxgl"

mmat=NULL

for (i in seq(length(raw_count_files))){
  mmat[[i]]=read.table(paste0("raw_counts/",raw_count_files[i]))
  colnames(mmat[[i]])=c("gene",sub("_.*","",raw_count_files[i]))
  rownames(mmat[[i]])=mmat[[i]]$gene
  mmat[[i]]$gene=NULL
}

mmat_df=do.call("cbind", mmat)

sample_annot=sample_annot[cn(mmat_df),]

mmat_mtx=as.matrix(mmat_df)

mmat_mtx_cpm=CPM(mmat_mtx)

# X/A median chrX genes/median autosomal genes

# all(cn(mmat_mtx_cpm)==rn(sample_annot))

mmat_mtx_cpm_x_genes=mmat_mtx_cpm[gene_chr$gene[gene_chr$chr=="chrX"],]
# mmat_mtx_cpm_x_genes=mmat_mtx_cpm_x_genes[rowSums(mmat_mtx_cpm_x_genes)>5,] # discard low-expressed genes

mmat_mtx_cpm_autosomal_genes=mmat_mtx_cpm[gene_chr$gene[!(gene_chr$chr%in%c("chrX","chrY"))],]
# mmat_mtx_cpm_autosomal_genes=mmat_mtx_cpm_autosomal_genes[rowMeans(mmat_mtx_cpm_autosomal_genes)>5,] # discard low-expressed genes

stat_x=apply(mmat_mtx_cpm_x_genes,2,function(x){median(as.numeric(x)[as.numeric(x)>0])}) # > 0 keep only expressed genes, otherwise 0 inflation, is absence of measure, not a measure
stat_a=apply(mmat_mtx_cpm_autosomal_genes,2,function(x){median(as.numeric(x)[as.numeric(x)>0])}) # > 0 keep only expressed genes, otherwise 0 inflation, is absence of measure, not a measure
stat_ratio=stat_x/stat_a

stat_4_plot=data.frame(sample=names(stat_x),
                       x_a_ratio=stat_ratio,
                       group=sample_annot$group,
                       cell_type=sample_annot$cell_type,
                       cell_type_group=paste(sep="_",sample_annot$cell_type,sample_annot$group))

stat_4_plot %>%
  # filter(cell_type !="primed") %>%
  group_by(cell_type_group) %>%
    mutate(sd=sd(x_a_ratio), mean=mean(x_a_ratio)) %>%
    ggplot() + 
  geom_bar(aes(x=cell_type_group,y=x_a_ratio,fill=cell_type_group),stat = 'summary',fun="mean",position = 'dodge')+
  geom_point(aes(x=cell_type_group,y=x_a_ratio,group=cell_type_group),position = position_dodge(width = 0.9))+
  # facet_wrap(~category, nrow=1, ncol=2)+
  geom_errorbar(aes(x=cell_type_group,y=x_a_ratio,group=cell_type_group,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
            theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
            ggtitle("X/A ratio")+
            #scale_x_discrete(expand = c(0, 0.5)) +
            ylab("X/A ratio") +
            guides(fill=guide_legend(ncol=1)) +
            theme(axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.ticks.x=element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank())+
    scale_y_continuous(breaks = seq(0,1.1, 0.1)) + ylim(0,1)

# ggsave("dataset_2/spen_x_a_ratio_0_1_all.pdf",width = 6,height = 7)


stat_4_plot %>%
  filter(cell_type !="primed") %>%
  group_by(cell_type_group) %>%
  mutate(sd=sd(x_a_ratio), mean=mean(x_a_ratio)) %>%
  ggplot() + 
  geom_bar(aes(x=cell_type_group,y=x_a_ratio,fill=cell_type_group),stat = 'summary',fun="mean",position = 'dodge')+
  geom_point(aes(x=cell_type_group,y=x_a_ratio,group=cell_type_group),position = position_dodge(width = 0.9))+
  # facet_wrap(~category, nrow=1, ncol=2)+
  geom_errorbar(aes(x=cell_type_group,y=x_a_ratio,group=cell_type_group,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("X/A ratio")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("X/A ratio") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(breaks = seq(0,1.1, 0.1)) + ylim(0,1)+
  stat_compare_means(aes(x=cell_type_group,y=x_a_ratio),vjust=1,label="p.signif",size=3,ref.group = "pxgl_siScr", method="t.test", paired=FALSE)


# ggsave("dataset_2/spen_x_a_ratio_0_1_pxgl.pdf",width = 6,height = 7)



########## all chr median expr

stat_chrom=NULL

for(chrom in unique(gene_chr$chr)[unique(gene_chr$chr)!='chrY']){
  
  
  df=mmat_mtx_cpm[gene_chr$gene[gene_chr$chr==chrom],]
  
  df=df[rowSums(df)>5,] # discard low-expressed genes
  stat_chrom[[chrom]]=colMedians(df)
  names(stat_chrom[[chrom]])=cn(mmat_mtx_cpm)
}

stat_chrom_df=as.data.frame(do.call('rbind',stat_chrom))
stat_chrom_df$chr=rn(stat_chrom_df)

stat_chrom_df %>%
  tidyr::pivot_longer(!chr, names_to = "sample", values_to = "median_expr") -> df_plot

# rn(sample_annot)==cn(stat_chrom_df) # T except for column chr

df_plot$cell_type=rep(sample_annot$cell_type,nrow(stat_chrom_df))

df_plot$group=rep(sample_annot$group,nrow(stat_chrom_df))

df_plot$cell_type_group=paste(sep="_",df_plot$cell_type,df_plot$group)

df_plot$chr=factor(df_plot$chr,levels=c(paste0('chr',seq(22)),'chrX'))

df_plot=as.data.frame(df_plot)

df_plot$cell_type_group=factor(df_plot$cell_type_group,levels=c("pxgl_siScr","pxgl_siSPEN_1","pxgl_siSPEN_2",
                                                                "primed_siScr","primed_siSPEN_1","primed_siSPEN_2"))


df_plot %>%
  filter(cell_type_group%in%c("pxgl_siScr","pxgl_siSPEN_1")) %>%
  list() -> df_plot_list

p_val_df=do.call('rbind',lapply(df_plot_list,function(x){
  
  x %>%
    group_by(chr) %>%
    group_split -> list_chr_pval
  
  names(list_chr_pval)=c(paste0('chr',seq(22)),'chrX') # checked manually, OK
  
  p_val_list=do.call('rbind',lapply(list_chr_pval,function(x){compare_means(x,formula=median_expr~cell_type_group,method='t.test',paired = F)}))
  
  df=cbind(x,rep(p_val_list$p.format,each=6),rep(p_val_list$p.signif,each=6))
  
  colnames(df)=c(cn(x),'p.format','p.signif')
  
  return(df)})
)

# pdf(file="dataset_2/spen_1_pxgl_barplot_median_expr_htseq.pdf",width=30,height=15)
p_val_df %>%
  group_by(cell_type_group,chr) %>%
  mutate(sd=sd(median_expr), mean=mean(median_expr)) %>%
  ggplot(aes(x=chr, fill=cell_type_group, group=cell_type_group)) + 
  geom_bar(aes(y=median_expr),stat = "summary", fun.y = "mean", position = "dodge")+
  geom_point(aes(y=median_expr),position = position_dodge(width = .9))+
  geom_errorbar(aes(y=median_expr,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("SPEN chrom median expression")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("median expression in CPM") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  geom_text(aes(y=25,label=p.signif),position = position_dodge(width = .9))

# dev.off()



df_plot %>%
  filter(cell_type_group%in%c("pxgl_siScr","pxgl_siSPEN_2")) %>%
  list() -> df_plot_list

p_val_df=do.call('rbind',lapply(df_plot_list,function(x){
  
  x %>%
    group_by(chr) %>%
    group_split -> list_chr_pval
  
  names(list_chr_pval)=c(paste0('chr',seq(22)),'chrX') # checked manually, OK
  
  p_val_list=do.call('rbind',lapply(list_chr_pval,function(x){compare_means(x,formula=median_expr~cell_type_group,method='t.test',paired = F)}))
  
  df=cbind(x,rep(p_val_list$p.format,each=6),rep(p_val_list$p.signif,each=6))
  
  colnames(df)=c(cn(x),'p.format','p.signif')
  
  return(df)})
)

# pdf(file="dataset_2/spen_2_pxgl_barplot_median_expr_htseq.pdf",width=30,height=15)
p_val_df %>%
  group_by(cell_type_group,chr) %>%
  mutate(sd=sd(median_expr), mean=mean(median_expr)) %>%
  ggplot(aes(x=chr, fill=cell_type_group, group=cell_type_group)) + 
  geom_bar(aes(y=median_expr),stat = "summary", fun.y = "mean", position = "dodge")+
  geom_point(aes(y=median_expr),position = position_dodge(width = .9))+
  geom_errorbar(aes(y=median_expr,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("SPEN chrom median expression")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("median expression in CPM") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  geom_text(aes(y=25,label=p.signif),position = position_dodge(width = .9))

# dev.off()


# pdf(file="dataset_2/spen_all_barplot_median_expr_htseq.pdf",width=30,height=15)

df_plot %>%
  group_by(cell_type_group,chr) %>%
  mutate(sd=sd(median_expr), mean=mean(median_expr)) %>%
  ggplot(aes(x=chr, fill=cell_type_group, group=cell_type_group)) + 
  geom_bar(aes(y=median_expr),stat = "summary", fun.y = "mean", position = "dodge")+
  geom_point(aes(y=median_expr),position = position_dodge(width = .9))+
  geom_errorbar(aes(y=median_expr,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("SPEN chrom median expression")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("median expression in CPM") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# dev.off()


# pdf(file="dataset_2/spen_pxgl_barplot_median_expr_htseq.pdf",width=30,height=15)

df_plot %>%
  filter(cell_type=="pxgl") %>%
  group_by(cell_type_group,chr) %>%
  mutate(sd=sd(median_expr), mean=mean(median_expr)) %>%
  ggplot(aes(x=chr, fill=cell_type_group, group=cell_type_group)) + 
  geom_bar(aes(y=median_expr),stat = "summary", fun.y = "mean", position = "dodge")+
  geom_point(aes(y=median_expr),position = position_dodge(width = .9))+
  geom_errorbar(aes(y=median_expr,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("SPEN chrom median expression")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("median expression in CPM") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# dev.off()


