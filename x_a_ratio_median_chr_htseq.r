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


gene_chr=read.table("hg38_gene_chr_list.txt")
gene_chr$V2=NULL
colnames(gene_chr)=c("gene","chr")

raw_count_files=list.files(path = '.',pattern = "_raw_counts.txt")
sample_annot=fastRead("sample_annot_all.txt", header=T, sep=",",as.matrix = F)
sample_annot$sample_id=rn(sample_annot)

# sample_annot$condition[sample_annot$condition=='C'] <- "Control"
# sample_annot$condition[sample_annot$condition=='T'] <- "XIST_KD"

mmat=NULL

for (i in seq(length(raw_count_files))){
  mmat[[i]]=read.table(raw_count_files[i])
  colnames(mmat[[i]])=c("gene",sub("_.*","",raw_count_files[i]))
  rownames(mmat[[i]])=mmat[[i]]$gene
  mmat[[i]]$gene=NULL
}

mmat_df=do.call("cbind", mmat)
mmat_df=mmat_df[gene_chr$gene,]

mmat_df=mmat_df[,!grepl('SRR529',cn(mmat_df))] # remove Vallot samples here

sample_annot=sample_annot[cn(mmat_df),]

mmat_mtx=as.matrix(mmat_df)

mmat_mtx_cpm=CPM(mmat_mtx)

variableGenes<-getMostVariableGenes4(mmat_mtx_cpm,minCount = 1)
thres<- 0.5
qplotDensity(variableGenes$residuals)+geom_vline(xintercept = thres)
variableGenesNames<-rn(variableGenes)[variableGenes$residuals>thres]
len(variableGenesNames)

mmat_df_all_genes=mmat_df

mmat_df_to_write=mmat_df_all_genes[gene_chr$gene,]
mmat_df_to_write$chr=unlist(sapply(rn(mmat_df_to_write), function(x){gene_chr$chr[gene_chr$gene==x]}))
# 
# write.tsv(mmat_df_to_write,
#           file="rna_gro_htseq_raw_counts_with_chr.tsv")


mmat_df=mmat_df[variableGenesNames,]
mmat_mtx=as.matrix(mmat_df)

# de_seq2

sample_annot %>%
  group_split(dataset) -> sample_annot_group

names(sample_annot_group)=c("D1271","D1269","D1269_bis","D165")

sample_annot_group=lapply(sample_annot_group, function(x){as.data.frame(x)})

for (i in names(sample_annot_group)){
  
  sample_annot_dds=sample_annot_group[[i]]
  rownames(sample_annot_dds)=sample_annot_dds$sample_id
  
  if (i=="D165"){
    
    dds=NULL
    dds <- DESeqDataSetFromMatrix(
      countData=mmat_mtx[,rn(sample_annot_dds)], #integers
      colData=sample_annot_dds,
      design=~cell_type)
    
    # tail(CPM(assay(dds))==CPM(mmat_mtx)) #TRUE
    
    dds <- estimateSizeFactors(dds)
    
    # vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    vsd <- vst(dds, blind=FALSE)
    
    
    dds <- DESeq(dds)
    res <- results(dds)

    resultsNames(dds)
    
    # resLFC <- lfcShrink(dds, coef="Treatment_UT_vs_T", type="apeglm")
    # resLFC
    
    resOrdered <- res[order(res$pvalue),]
    
    summary(res)
    
    sum(res$padj < 0.05, na.rm=TRUE)
    
    res05 <- results(dds, alpha=0.05)
    
    summary(res05)
    
    EnhancedVolcano(res,pCutoff =0.05,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title='DEseq2 naive vs primed')
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    # # dev.off()
    # print(plotMA(res, ylim=c(-10,10)))
    # ggsave(paste(sep="_",i,"plotMA_rna_paper_deseq2.pdf"))
    # dev.off()
    
    # plotMA(resLFC, ylim=c(-5,5))
    
    resSig <- subset(resOrdered, padj < 0.05)

    
    # write.tsv(mmat_mtx_cpm,
    #           file=paste(sep="_",i,"rna_paper_cpm.tsv"))
    # 
    # write.tsv(assay(vsd),
    #           file=paste(sep="_",i,"rna_paper_vsd.tsv"))

    de_seq_signig_cpm=as.data.frame(mmat_mtx_cpm[rn(resSig),])
    de_seq_signig_cpm$chr=unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]}))
    # names(unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]})))==rn(resSig) # NOT ALL TRUE ! BUT OK, only cause .1 .2 changes ...
    
    # write.tsv(de_seq_signig_cpm,
    #           file=paste(sep="_",i,"naive_primed_de_seq2_signif_gene_expr_cpm.tsv"))

    # write.tsv(assay(vsd)[rn(resSig),],
    #           file=paste(sep="_",i,"rna_paper_de_seq2_signif_gene_expr_vsd.tsv"))

    # write.tsv(res,
    #           file=paste(sep="_",i,"rna_paper_de_seq2.tsv"))

  }
  if (i=="D1269"){
    
    dds=NULL
    
    dds <- DESeqDataSetFromMatrix(
      countData=mmat_mtx[,rn(sample_annot_dds)], #integers
      colData=sample_annot_dds,
      design=~condition)
    
    # tail(CPM(assay(dds))==CPM(mmat_mtx)) #TRUE
    
    dds <- estimateSizeFactors(dds)
    
    # vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    vsd <- vst(dds, blind=FALSE)
    
    
    dds <- DESeq(dds)
    res <- results(dds)

    resultsNames(dds)
    
    # resLFC <- lfcShrink(dds, coef="Treatment_UT_vs_T", type="apeglm")
    # resLFC
    
    resOrdered <- res[order(res$pvalue),]
    
    summary(res)
    
    sum(res$padj < 0.05, na.rm=TRUE)
    
    res05 <- results(dds, alpha=0.05)
    
    summary(res05)
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    
    EnhancedVolcano(res,pCutoff =0.05,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title="DEseq2 PXGL XIST DOX control vs treatment")
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    # # dev.off()
    # print(plotMA(res, ylim=c(-10,10)))
    # ggsave(paste(sep="_",i,"plotMA_rna_paper_deseq2.pdf"))
    # dev.off()
    
    resSig <- subset(resOrdered, padj < 0.05)

    
    # write.tsv(mmat_mtx_cpm,
    #           file=paste(sep="_",i,"rna_paper_cpm.tsv"))
    # 
    # write.tsv(assay(vsd),
    #           file=paste(sep="_",i,"rna_paper_vsd.tsv"))

    de_seq_signig_cpm=as.data.frame(mmat_mtx_cpm[rn(resSig),])
    de_seq_signig_cpm$chr=unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]}))
    # names(unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]})))==rn(resSig) # NOT ALL TRUE ! BUT OK, only cause .1 .2 changes ...
    
    # write.tsv(de_seq_signig_cpm,
    #           file=paste(sep="_",i,"pxgl_dox_de_seq2_signif_gene_expr_cpm.tsv"))

    # write.tsv(assay(vsd)[rn(resSig),],
    #           file=paste(sep="_",i,"rna_paper_de_seq2_signif_gene_expr_vsd.tsv"))

    # write.tsv(res,
    #           file=paste(sep="_",i,"rna_paper_de_seq2.tsv"))
  }
  if (i=="D1269_bis"){
    
    dds <- DESeqDataSetFromMatrix(
      countData=mmat_mtx[,rn(sample_annot_dds)], #integers
      colData=sample_annot_dds,
      design=~condition)
    
    # tail(CPM(assay(dds))==CPM(mmat_mtx)) #TRUE
    
    dds <- estimateSizeFactors(dds)
    
    # vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    vsd <- vst(dds, blind=FALSE)
    
    
    dds <- DESeq(dds)
    res <- results(dds)

    resultsNames(dds)
    
    # resLFC <- lfcShrink(dds, coef="Treatment_UT_vs_T", type="apeglm")
    # resLFC
    
    resOrdered <- res[order(res$pvalue),]
    
    summary(res)
    
    sum(res$padj < 0.05, na.rm=TRUE)
    
    res05 <- results(dds, alpha=0.05)
    
    summary(res05)
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    
    EnhancedVolcano(res,pCutoff =0.05,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title="DEseq2 PXGL XIST KO control vs treatment")
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    # # dev.off()
    # print(plotMA(res, ylim=c(-10,10)))
    # ggsave(paste(sep="_",i,"plotMA_rna_paper_deseq2.pdf"))
    # dev.off()
    
    resSig <- subset(resOrdered, padj < 0.05)

    
    #   write.tsv(mmat_mtx_cpm, 
    #             file=paste(sep="_",i,"rna_paper_cpm.tsv"))
    #   
    #   write.tsv(assay(vsd), 
    #             file=paste(sep="_",i,"rna_paper_vsd.tsv"))
    #   
    de_seq_signig_cpm=as.data.frame(mmat_mtx_cpm[rn(resSig),])
    de_seq_signig_cpm$chr=unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]}))
    # names(unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]})))==rn(resSig) # NOT ALL TRUE ! BUT OK, only cause .1 .2 changes ...
    
    # write.tsv(de_seq_signig_cpm,
    #             file=paste(sep="_",i,"pxgl_wt_ko_de_seq2_signif_gene_expr_cpm.tsv"))
    #   
    #   write.tsv(assay(vsd)[rn(resSig),], 
    #             file=paste(sep="_",i,"rna_paper_de_seq2_signif_gene_expr_vsd.tsv"))
    #   
    #   write.tsv(res, 
    #             file=paste(sep="_",i,"rna_paper_de_seq2.tsv"))    
  }
  if(i=="D1271"){
    
    dds <- DESeqDataSetFromMatrix(
      countData=mmat_mtx[,rn(sample_annot_dds)], #integers
      colData=sample_annot_dds,
      design=~condition)
    
    # tail(CPM(assay(dds))==CPM(mmat_mtx)) #TRUE
    
    dds <- estimateSizeFactors(dds)
    
    # vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    vsd <- vst(dds, blind=FALSE)
    
    
    dds <- DESeq(dds)
    res <- results(dds)

    resultsNames(dds)
    
    # resLFC <- lfcShrink(dds, coef="Treatment_UT_vs_T", type="apeglm")
    # resLFC
    
    resOrdered <- res[order(res$pvalue),]
    
    summary(res)
    
    sum(res$padj < 0.05, na.rm=TRUE)
    
    res05 <- results(dds, alpha=0.05)
    
    summary(res05)
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    
    EnhancedVolcano(res,pCutoff =0.05,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title="GRO-SEQ PXGL XIST DOX control vs treatment")
    
    sum(res05$padj < 0.05, na.rm=TRUE)
    # # dev.off()
    # print(plotMA(res, ylim=c(-10,10)))
    # ggsave(paste(sep="_",i,"plotMA_rna_paper_deseq2.pdf"))
    # dev.off()
    
    resSig <- subset(resOrdered, padj < 0.05)

    
    #   write.tsv(mmat_mtx_cpm, 
    #             file=paste(sep="_",i,"rna_paper_cpm.tsv"))
    #   
    #   write.tsv(assay(vsd), 
    #             file=paste(sep="_",i,"rna_paper_vsd.tsv"))
    #   
    de_seq_signig_cpm=as.data.frame(mmat_mtx_cpm[rn(resSig),])
    de_seq_signig_cpm$chr=unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]}))
    # names(unlist(sapply(rn(resSig), function(x){gene_chr$chr[gene_chr$gene==x]})))==rn(resSig) # NOT ALL TRUE ! BUT OK, only cause .1 .2 changes ...
    
    # write.tsv(de_seq_signig_cpm,
    #             file=paste(sep="_",i,"gro_seq_de_seq2_signif_gene_expr_cpm.tsv"))
    #   
    #   write.tsv(assay(vsd)[rn(resSig),], 
    #             file=paste(sep="_",i,"rna_paper_de_seq2_signif_gene_expr_vsd.tsv"))
    #   
    #   write.tsv(res, 
    #             file=paste(sep="_",i,"rna_paper_de_seq2.tsv"))    
  }
}

# X/A median chrX genes/median autosomal genes

mmat_mtx_cpm_x_genes=mmat_mtx_cpm[gene_chr$gene[gene_chr$chr=="chrX"],]
mmat_mtx_cpm_x_genes=mmat_mtx_cpm_x_genes[rowSums(mmat_mtx_cpm_x_genes)>5,] # discard low-expressed genes

mmat_mtx_cpm_autosomal_genes=mmat_mtx_cpm[gene_chr$gene[!(gene_chr$chr%in%c("chrX","chrY"))],]
mmat_mtx_cpm_autosomal_genes=mmat_mtx_cpm_autosomal_genes[rowMeans(mmat_mtx_cpm_autosomal_genes)>5,] # discard low-expressed genes



stat_x=colMedians(mmat_mtx_cpm_x_genes)
stat_a=colMedians(mmat_mtx_cpm_autosomal_genes)
stat_ratio=stat_x/stat_a

stat_4_plot=NULL
stat_4_plot$sample=names(stat_x)
stat_4_plot$x_a_ratio=stat_ratio
stat_4_plot$cell_type=sample_annot$cell_type
stat_4_plot$condition=sample_annot$condition
stat_4_plot$dataset=sample_annot$dataset
stat_4_plot$name=sample_annot$Sample_Name
stat_4_plot_df=do.call('cbind',stat_4_plot)
stat_4_plot_df=as.data.frame(stat_4_plot_df)

stat_4_plot_df %>%
  group_split(dataset) -> stat_plot_list


for (i in seq(length(stat_plot_list))){
  
  plot_df=stat_plot_list[[i]]
  plot_df=as.data.frame(plot_df)
  plot_df$x_a_ratio=as.numeric(plot_df$x_a_ratio)
  
  plot_df %>%
    group_by(cell_type,condition) %>%
    mutate(sd=sd(x_a_ratio), mean=mean(x_a_ratio)) %>%
    as.data.frame() -> plot_df
  
  if (i==1){
    
    print(ggplot(plot_df, aes(x=condition, y=x_a_ratio, fill=condition, group=condition)) + 
            geom_bar(stat = "summary", fun.y = "mean", position = "dodge")+
            geom_point(size=0.25,alpha=0.8)+
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge") +
            theme_bw()+
            theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
            ggtitle("GRO-seq naive XIST DOX+ vs DOX-")+
            #scale_x_discrete(expand = c(0, 0.5)) +
            ylab("X/A ratio") +
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
            stat_compare_means(vjust=1,label="p.signif",size=3,method="wilcox.test", paired=FALSE))
    
    ggsave("barplot_gro_seq_x_a_ratio_htseq.pdf")
    # dev.off()
    
    
  }
  if (i==2){
    
    print(ggplot(plot_df, aes(x=condition, y=x_a_ratio, fill=condition, group=condition)) +
            geom_bar(stat = "summary", fun.y = "mean", position = "dodge")+
            geom_point(size=0.25,alpha=0.8)+
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge") +
            theme_bw()+
            theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
            ggtitle("RNA-seq PXGL XIST DOX+ vs DOX-")+
            #scale_x_discrete(expand = c(0, 0.5)) +
            ylab("X/A ratio") +
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
            stat_compare_means(vjust=1,label="p.signif",size=3,method="wilcox.test", paired=FALSE))
    
    ggsave("barplot_rna_seq_pxgl_dox_x_a_ratio_htseq.pdf")
    # dev.off()
    
    
  }
  if (i==3){
    
    print(ggplot(plot_df, aes(x=condition, y=x_a_ratio, fill=condition, group=condition)) +
            geom_bar(stat = "summary", fun.y = "mean", position = "dodge")+
            geom_point(size=0.25,alpha=0.8)+
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge") +
            theme_bw()+
            theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
            ggtitle("RNA-seq PXGL XIST KO vs WT")+
            #scale_x_discrete(expand = c(0, 0.5)) +
            ylab("X/A ratio") +
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
            stat_compare_means(vjust=1,label="p.signif",size=3,method="wilcox.test", paired=FALSE))
    
    ggsave("barplot_rna_seq_pxgl_ko_vs_wt_x_a_ratio_htseq.pdf")
    # dev.off()
    
    
  }
  if (i==4){
    
    print(ggplot(plot_df, aes(x=cell_type, y=x_a_ratio, fill=cell_type, group=cell_type)) +
            geom_bar(stat = "summary", fun.y = "mean", position = "dodge")+
            geom_point(size=0.25,alpha=0.8)+
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge") +
            theme_bw()+
            theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
            ggtitle("RNA-seq Naive vs Primed")+
            #scale_x_discrete(expand = c(0, 0.5)) +
            ylab("X/A ratio") +
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
            stat_compare_means(vjust=1,label="p.signif",size=3,method="wilcox.test", paired=FALSE))
    
    ggsave("barplot_rna_seq_naive_vs_primed_x_a_ratio_htseq.pdf")
    # dev.off()
    
  }
}


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

df_plot$dataset=rep(sample_annot$dataset,nrow(stat_chrom_df))
df_plot$cell_type=rep(sample_annot$cell_type,nrow(stat_chrom_df))
df_plot$condition=rep(sample_annot$condition,nrow(stat_chrom_df))
df_plot$dox_ko=rep(sample_annot$dox_ko,nrow(stat_chrom_df))
df_plot$dataset_simple=df_plot$dataset

df_plot$dataset_simple[grep('rna',df_plot$dataset_simple)]<-'rna'
df_plot$dataset_simple[grep('gro',df_plot$dataset_simple)]<-'gro'

df_plot$dataset_cellType_condition=paste(sep="_",df_plot$dataset_simple,df_plot$cell_type,df_plot$dox_ko)

df_plot$chr=factor(df_plot$chr,levels=c(paste0('chr',seq(22)),'chrX'))

df_plot=as.data.frame(df_plot)
df_plot=df_plot[!grepl('D165',df_plot$sample),]

df_plot$dataset_cellType_condition=factor(df_plot$dataset_cellType_condition,levels=c('gro_naive_C',
                                                                                      'gro_naive_T',
                                                                                      'rna_naive_C',
                                                                                      'rna_naive_T',
                                                                                      'rna_naive_WT',
                                                                                      'rna_naive_KO'))




my_comparisons=list(c("gro_naive_C","gro_naive_T"),c('rna_naive_C','rna_naive_T'),
                    c("rna_naive_WT","rna_naive_KO"))


df_plot %>%
  group_by(dataset)%>%
  group_split()->df_plot_list

p_val_df=do.call('rbind',lapply(df_plot_list,function(x){

  x %>%
  group_by(chr) %>%
  group_split -> list_chr_pval

  names(list_chr_pval)=c(paste0('chr',seq(22)),'chrX') # checked manually, OK

  p_val_list=do.call('rbind',lapply(list_chr_pval,function(x){compare_means(x,formula=median_expr~dataset_cellType_condition,method='t.test',paired = F)}))

  df=cbind(x,rep(p_val_list$p.format,each=6),rep(p_val_list$p.signif,each=6))

  colnames(df)=c(cn(x),'p.format','p.signif')

  return(df)})
)

pdf(file="barplot_median_expr_htseq.pdf",width=30,height=15)
p_val_df %>%
  group_by(dataset_cellType_condition,chr) %>%
  mutate(sd=sd(median_expr), mean=mean(median_expr)) %>%
  ggplot(aes(x=chr, fill=dataset_cellType_condition, group=dataset_cellType_condition)) + 
            geom_bar(aes(y=median_expr),stat = "summary", fun.y = "mean", position = "dodge")+
            geom_point(aes(y=median_expr),position = position_dodge(width = .9))+
            geom_errorbar(aes(y=median_expr,ymin=mean-sd, ymax=mean+sd), position="dodge") +
            theme_bw()+
            theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
            ggtitle("chrom median expression")+
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

dev.off()
  
# ggsave("barplot_gro_seq_median_expr_htseq.pdf")
# dev.off()
    
