setwd('C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/')

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

gene_chr=read.table("C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/hg38_gene_chr_list.txt")
gene_chr$V2=NULL
colnames(gene_chr)=c("gene","chr")

raw_count_files=list.files(path = 'dataset_2/',pattern = "_raw_counts.txt")
sample_annot=read.csv("dataset_2/SampleSheet.csv", header=T, comment.char = '#')
rownames(sample_annot)=paste0("D1507",sample_annot$Sample_ID)
sample_annot$sample_id=paste0("D1507",sample_annot$Sample_ID)
sample_annot$group=sub('_[^_]*$', '', sample_annot$Sample_Name)
sample_annot$group=str_replace_all(sample_annot$group,"CTL","1")
sample_annot$group=str_replace_all(sample_annot$group,"Mix","2")
sample_annot$replicate=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$replicate[sample_annot$replicate=="primed"]<-"Rep1"
sample_annot$cell_type=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$cell_type[sample_annot$cell_type!="primed"]<-"pxgl"

mmat=NULL

for (i in seq(length(raw_count_files))){
  mmat[[i]]=read.table(paste0("dataset_2/",raw_count_files[i]))
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

ggsave("dataset_2/spen_x_a_ratio_0_1_all.pdf",width = 6,height = 7)


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


ggsave("dataset_2/spen_x_a_ratio_0_1_pxgl.pdf",width = 6,height = 7)
