setwd('C:/Users/surfi/Desktop/charbel_paper_2024/2024_04/charbel_paper_github_v1/rna/')

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))


# BiocManager::install('EnhancedVolcano')
# BiocManager::install("sva")
# BiocManager::install("gtools")

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
library(biomaRt)
library(sva)
library(gtools)


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

mmat_mtx_cpm=mmat_mtx_cpm[c("XIST","NANOG","POU5F1","POU5F1B","TFCP2L1","KLF5","KLF4","DNMT3L","SPEN"),]

stat_4_plot=data.frame(value=as.vector(mmat_mtx_cpm),
                       sample=rep(cn(mmat_mtx_cpm),each=nrow(mmat_mtx_cpm)),
                       gene=rep(rn(mmat_mtx_cpm),ncol(mmat_mtx_cpm)),
                       cell_type_group=rep(paste(sep="_",sample_annot$cell_type,sample_annot$group),each=nrow(mmat_mtx_cpm)))


stat_4_plot %>%
  group_by(cell_type_group,gene) %>%
  mutate(sd=sd(value), mean=mean(value)) %>%
  ggplot() + 
  geom_bar(aes(x=cell_type_group,y=value,fill=cell_type_group),stat = 'summary',fun="mean",position = 'dodge')+
  geom_point(aes(x=cell_type_group,y=value,group=cell_type_group),position = position_dodge(width = 0.9))+
  facet_wrap(~gene, nrow=5, ncol=5,scales = "free")+
  geom_errorbar(aes(x=cell_type_group,y=value,group=cell_type_group,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("Markers expression (CPM)")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("Number of transcripts per million of mRNA molecules") +
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
  scale_y_continuous(limits = c(0,NA))

# ggsave("dataset_2/spen_dataset_2_markers_expression_all.pdf",width = 10,height = 10)


stat_4_plot %>%
  filter(cell_type_group %in% c("pxgl_siScr","pxgl_siSPEN_1","pxgl_siSPEN_2")) %>%
  group_by(cell_type_group,gene) %>%
  mutate(sd=sd(value), mean=mean(value)) %>%
  ggplot() + 
  geom_bar(aes(x=cell_type_group,y=value,fill=cell_type_group),stat = 'summary',fun="mean",position = 'dodge')+
  geom_point(aes(x=cell_type_group,y=value,group=cell_type_group),position = position_dodge(width = 0.9))+
  facet_wrap(~gene, nrow=5, ncol=5,scales = "free")+
  geom_errorbar(aes(x=cell_type_group,y=value,group=cell_type_group,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("Markers expression (CPM)")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("Number of transcripts per million of mRNA molecules") +
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
  scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means(aes(x=cell_type_group,y=value),vjust=1,label="p.signif",size=3,ref.group = "pxgl_siScr", method="t.test", paired=FALSE)

# ggsave("dataset_2/spen_dataset_2_markers_expression_pxgl_stat.pdf",width = 10,height = 10)
