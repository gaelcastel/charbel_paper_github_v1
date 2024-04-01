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
library(rstatix)

bi_bi_up=read.table("C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/embryo_pseudobulk/bi_bi_up.gtf", header=F, sep="\t")
colnames(bi_bi_up)=c('chrom','source','feature','start',
                     'end',
                     'score',
                     'strand',
                     'frame',
                     'attribute')

mono_bi_up=read.table("C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/embryo_pseudobulk/mono_bi_up.gtf", header=F, sep="\t")
colnames(mono_bi_up)=c('chrom','source','feature','start',
                       'end',
                       'score',
                       'strand',
                       'frame',
                       'attribute')

bi_bi=read.table("C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/embryo_pseudobulk/bi_bi.gtf", header=F, sep="\t")
colnames(bi_bi)=c('chrom','source','feature','start',
                  'end',
                  'score',
                  'strand',
                  'frame',
                  'attribute')

up_gtf=rbind(bi_bi_up,
             mono_bi_up)

constant_gtf=bi_bi # I checked, these genes are common to wt-ko/ctl-ttt, cf file: C:\Users\gael\charbel_paper\rna\charbel_paper_ifb\embryo_pseudobulk\charbel_petro_zhou_v2.r


up_gtf%>%
  filter(feature=="gene")->up_gtf

up_gtf$chrom[up_gtf$chrom=="X"]="chrX"
up_gtf$category=rep('up',nrow(up_gtf))
up_gtf$gene=sub('.*gene_name ','',unlist(lapply(str_split(pattern=";",up_gtf$attribute),function(x){x[grep('gene_name',x)]})))

constant_gtf%>%
  filter(feature=="gene")->constant_gtf

constant_gtf$chrom[constant_gtf$chrom=="X"]="chrX"
constant_gtf$category=rep('constant',nrow(constant_gtf))
constant_gtf$gene=sub('.*gene_name ','',unlist(lapply(str_split(pattern=";",constant_gtf$attribute),function(x){x[grep('gene_name',x)]})))

gene_cat=c(up_gtf$gene,constant_gtf$gene)
gene_category=NULL

for(gene in gene_cat){
  
  if(gene %in% up_gtf$gene){gene_category[[gene]]="up"}
  if(gene %in% constant_gtf$gene){gene_category[[gene]]="constant"}
  
}

gene_category=unlist(gene_category)
# all(names(gene_category)==gene_cat)

gene_cat=data.frame(gene=gene_cat,category=gene_category)

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

sample_annot=sample_annot[cn(mmat_df),]

mmat_mtx=as.matrix(mmat_df)

mmat_mtx_cpm=CPM(mmat_mtx)

samples=rn(sample_annot)[sample_annot$cell_type%in%c("naive","primed")
                         & sample_annot$dox_ko=="WT"
                         & sample_annot$dataset!="rna_seq_dvallot"]

dataForViolinBU=reshape2::melt(data.frame(t(log2(mmat_mtx_cpm[gene_cat$gene,samples]+1)),sample=sample_annot[samples,"Sample_Name"],
                                          cell_type=sample_annot[samples,"cell_type"],check.names = F),
                               value.name="log_cpm",variable.name="gene",id.vars=c("sample","cell_type"))

dataForViolinBU$category=sapply(dataForViolinBU$gene,function(x){
  gene_cat[x,"category"]
  
})


pdf("gene_categories_expression_in_stem_cells_violin_plots_wt.pdf",height = 8,width = 5)
print(ggplot(dataForViolinBU,aes(x=category,y=log_cpm,fill=cell_type))+
        geom_violin(scale = "width")+
        geom_boxplot(outlier.shape = NA,width=0.3,position=position_dodge(.9))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,size=8),
              panel.background = element_rect(fill = NA,colour="black"),
              panel.grid.major = element_line(colour = "black"),
              panel.grid.minor = element_line(colour = NA),
              panel.grid.major.x = element_line(colour = NA),
              strip.text = element_text(face = "bold.italic"),
              strip.background = element_rect(fill=NA,color="black")
        )+
        ggtitle("Gene categories expression in stem cells")+
        stat_compare_means(label.y = 9.8,label="p",size=3, method="t.test", paired=FALSE)+
        stat_compare_means(label.y = 10,label="p.signif",size=3, method="t.test", paired=FALSE)   
)
dev.off()



samples=rn(sample_annot)[sample_annot$dataset!="rna_seq_dvallot"]

dataForViolinBU=reshape2::melt(data.frame(t(log2(mmat_mtx_cpm[gene_cat$gene,samples]+1)),sample=sample_annot[samples,"Sample_Name"],
                                          dataset=sample_annot[samples,"dataset"],
                                          cell_type_condition_dataset=paste(sep="_",sample_annot[samples,"cell_type"],
                                                                            sample_annot[samples,"dox_ko"],
                                                                            sample_annot[samples,"dataset"]),check.names = F),
                               value.name="log_cpm",variable.name="gene",id.vars=c("sample","dataset","cell_type_condition_dataset"))

dataForViolinBU$category=sapply(dataForViolinBU$gene,function(x){
  gene_cat[x,"category"]
  
})

dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_C_rna_seq_d1269"]<-"naive_C"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_T_rna_seq_d1269"]<-"naive_T"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_WT_rna_seq_d1269_bis"]<-"naive_WT"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_KO_rna_seq_d1269_bis"]<-"naive_KO"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_C_gro_seq_d1271"]<-"naive_C_gro"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_T_gro_seq_d1271"]<-"naive_T_gro"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="primed_WT_rna_seq_d165"]<-"primed_WT"
dataForViolinBU$cell_type_condition_dataset[dataForViolinBU$cell_type_condition_dataset=="naive_cult_WT_rna_seq_d165"]<-"naive_cult_WT"

dataForViolinBU$cell_type_condition_dataset=factor(dataForViolinBU$cell_type_condition_dataset,levels=c("naive_WT",
                                                                                                        "naive_KO",
                                                                                                        "naive_C",
                                                                                                        "naive_T",
                                                                                                        "naive_C_gro",
                                                                                                        "naive_T_gro",
                                                                                                        "naive_cult_WT",
                                                                                                        "primed_WT"))

# my_comparisons=list(c("naive_WT","naive_KO"),
#                     c("naive_C","naive_T"),
#                     c("naive_C_gro","naive_T_gro"),
#                     c("naive_cult_WT","primed_WT"))

p<-ggplot(dataForViolinBU)+
        geom_violin(aes(x=category,y=log_cpm,fill=cell_type_condition_dataset),scale = "width")+
        geom_boxplot(aes(x=category,y=log_cpm,fill=cell_type_condition_dataset),outlier.shape = NA,width=0.3,position=position_dodge(.9))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,size=8),
              panel.background = element_rect(fill = NA,colour="black"),
              panel.grid.major = element_line(colour = "black"),
              panel.grid.minor = element_line(colour = NA),
              panel.grid.major.x = element_line(colour = NA),
              strip.text = element_text(face = "bold.italic"),
              strip.background = element_rect(fill=NA,color="black")
        )+
        ggtitle("Gene categories expression in stem cells")
        # stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE)

res.stat <- dataForViolinBU %>% 
  group_by(category,dataset)%>%
  t_test(data = ., log_cpm ~ cell_type_condition_dataset,paired = F) %>% 
  # adjust_pvalue(method = "bonferroni") %>%
  # add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "category", fun = "max",formula = log_cpm ~ cell_type_condition_dataset,dodge = 0.9) %>%
  ungroup()


pdf("gene_categories_expression_in_stem_cells_violin_plots_all_cell_types_pval.pdf",height = 6,width = 8)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif", y.position=res.stat$y.position+2
        )+
        stat_pvalue_manual(
          res.stat, label = "p", tip.length = 0.01
        ))

dev.off()
