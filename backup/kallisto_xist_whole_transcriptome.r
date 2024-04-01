setwd('C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/')

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))


# BiocManager::install('EnhancedVolcano')
# BiocManager::install("sva")
# BiocManager::install("gtools")

library(biomaRt)
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
library(sva)
library(gtools)


gene_chr=read.table("hg38_gene_chr_list.txt")
gene_chr$V2=NULL
colnames(gene_chr)=c("gene","chr")

sample_annot=fastRead("sample_annot_all.txt", header=T, sep=",",as.matrix = F)
sample_annot$sample_id=rn(sample_annot)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)

results <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "chromosome_name"),mart = mart)

# use output from kallisto

raw_count_files=list.files(path = 'kallisto_quant_xist_rf_stranded_whole_transcriptome/',pattern = "_abundance.tsv")

# sample_annot$condition[sample_annot$condition=='C'] <- "Control"
# sample_annot$condition[sample_annot$condition=='T'] <- "XIST_KD"

mmat=NULL

for (i in raw_count_files){
  mmat[[i]]=read.table(paste0('kallisto_quant_xist_rf_stranded_whole_transcriptome/',i),header=T)
  colnames(mmat[[i]])=c("target_id","length","eff_length","est_counts","tpm")
  rownames(mmat[[i]])=mmat[[i]]$target_id
  mmat[[i]]$target_id=NULL
  mmat[[i]]$length=NULL
  mmat[[i]]$eff_length=NULL
  mmat[[i]]$est_counts=NULL
  colnames(mmat[[i]])=sub("_.*","",i)
}

mmat_df=do.call("cbind", mmat)

sample_annot=sample_annot[cn(mmat_df),]

mmat_mtx=as.matrix(mmat_df)

gene_list=sapply(rn(mmat_mtx),function(x){
  
  list_gene=results$hgnc_symbol[results$ensembl_transcript_id==sub(":.*","",x)]
  print(list_gene)
  return(list_gene)
  
}) #takes time ~10 min

gene_list_corrected=NULL

for (name in names(gene_list)){
  
  gene_list_corrected[[name]]=ifelse(length(gene_list[[name]])==0,name,gene_list[[name]]) #character(0) when gene symbol not found
  print(gene_list_corrected[[name]])
} #takes time ~10 min

names(gene_list_corrected)=names(gene_list)

genes_corrected=unlist(gene_list_corrected)

mmat_plot=as.data.frame(mmat_mtx)

mmat_plot$gene=genes_corrected

# mmat_plot[grepl('_XIST|NANOG|DNMT3L|CD24',mmat_plot$gene),c("gene",rn(sample_annot)[sample_annot$dox_ko=="WT"&sample_annot$cell_type%in%c("primed","naive")])]
# mmat_plot_subset=mmat_plot[grepl('_XIST|ENST00000229307|ENST00000606017|ENST00000628202',rn(mmat_plot)),c("gene",rn(sample_annot)[sample_annot$dox_ko=="WT"&sample_annot$cell_type%in%c("primed","naive")])]

# ENST00000229307 : NANOG
# ENST00000606017 : CD24
# ENST00000628202 : DNMT3L

mmat_plot_subset=mmat_plot[grepl('_XIST',rn(mmat_plot)),c("gene",rn(sample_annot)[sample_annot$dox_ko=="WT"&sample_annot$cell_type%in%c("primed","naive")])]


mmat_plot_subset %>%
  tidyr::pivot_longer(!gene, names_to = "sample", values_to = "tpm") -> mmat_plot_subset
  
mmat_plot_subset$cell_type=mmat_plot_subset$sample
mmat_plot_subset$cell_type[grep('D165',mmat_plot_subset$cell_type)]<-"primed"
mmat_plot_subset$cell_type[grep('D1269',mmat_plot_subset$cell_type)]<-"naive"

mmat_plot_subset %>%
  ggplot() +
  geom_boxplot(aes(x=cell_type, y=tpm, fill=cell_type),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  geom_point(aes(x=cell_type, y=tpm, fill=cell_type),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~gene, nrow=4, ncol=4,scales = 'free')+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
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
  # stat_boxplot(aes(x=cell_type, y=tpm),geom='errorbar', linetype=1, width=0.8, size=0.25)+
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))+
  stat_compare_means(aes(x=cell_type, y=tpm),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("XIST isoforms in WT naive vs primed\nKallisto")

ggsave("xist_isoforms_kallisto_whole_transcriptome.pdf",width = 15,height = 10)
