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


gene_chr=read.table("hg38_gene_chr_list.txt")
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

# mmat_mtx_cpm[c("DNMT3L","SPEN","CD24","NANOG"),]
# all markers consistent, but SPEN only reduced by factor ~2

# variableGenes<-getMostVariableGenes4(mmat_mtx_cpm,minCount = 1)
# thres<- 0.5
# qplotDensity(variableGenes$residuals)+geom_vline(xintercept = thres)
# variableGenesNames<-rn(variableGenes)[variableGenes$residuals>thres]
# len(variableGenesNames)

mmat_df_cpm=as.data.frame(mmat_mtx_cpm)

# mmat_df=mmat_df[variableGenesNames,]
mmat_df=mmat_df[rowSums(mmat_df)>0,]
mmat_mtx=as.matrix(mmat_df)

# de_seq2

sample_annot %>%
  filter(cell_type!="primed") %>% # only one replicate, not possible DE-Seq2
  filter(group%in%c('siScr','siSPEN_1')) -> sample_annot_dds

dds=NULL
dds <- DESeqDataSetFromMatrix(
  countData=mmat_mtx[,rn(sample_annot_dds)], #integers
  colData=sample_annot_dds,
  design=~group)

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind=FALSE)


dds <- DESeq(dds)
res <- results(dds)

resultsNames(dds)

resOrdered <- res[order(res$pvalue),]

summary(res)

sum(res$padj < 0.05, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)

summary(res05)

print(EnhancedVolcano(res,labSize=2,drawConnectors=T,max.overlaps = 100,pCutoff =0.05,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title=paste(sep="_",unique(sample_annot_dds$group)[1],"vs",unique(sample_annot_dds$group)[2])))

sum(res05$padj < 0.05, na.rm=TRUE)

# ggsave(paste0("",paste(sep="_",unique(sample_annot_dds$group)[1],"vs",unique(sample_annot_dds$group)[2],"volcano_plot_deseq2.pdf")),width = 10,height = 10)
# dev.off()

resSig <- subset(resOrdered, padj < 0.05)

de_seq_signig_cpm=as.data.frame(mmat_mtx_cpm[rn(resSig),])
# de_seq_signig_cpm$chr=unlist(unique(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]})))

gene_chr_res=gene_chr

gene_chr_res %>%
  group_by(gene) %>%
  summarise(chr=paste(unique(chr), collapse=";")) %>%
  as.data.frame() -> gene_chr_res


rownames(gene_chr_res)=gene_chr_res$gene

gene_chr_res=gene_chr_res[rn(res),]

res$chr=gene_chr_res$chr

# write.tsv(res,
#           file=paste0("",paste(sep="_",unique(sample_annot_dds$group)[1],"vs",unique(sample_annot_dds$group)[2],"rna_paper_de_seq2.tsv")))


