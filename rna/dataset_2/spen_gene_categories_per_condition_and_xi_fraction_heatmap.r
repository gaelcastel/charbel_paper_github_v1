# vcf from bcftools already on SNP in H9 CHR X !!!!! BUT CAUTION !!! :not only heterozygous SNP !!! so I applied detect_bi on informative_H9_wgs.vcf
# Could be relevant to re-run bcftools intersect on this vcf -> bed only informative SNP heterozygous chrX to re-generate samples.vcf

setwd('C:/Users/surfi/Desktop/charbel_paper_2024/2024_04/charbel_paper_github_v1/rna/')

# ex. biallelic: DP4=0,3,2,0 : 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles
# ex. monoallelic: ; DP4=0,0,1,1: 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles. BUT also DP4=1,1,10,1 (proba/stat)

# NB: each dataset distinct read length, and some are single-end, so I guess need to determine heterozygous SNPs
# in primed of each dataset

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

# BiocManager::install("statebins")
# BiocManager::install("wesanderson")
library(statebins)
library(wesanderson)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
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
library(ComplexHeatmap)

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


spen_ctl_ttt_2_snp=fastRead('dataset_2/pxgl_siScr_pxgl_siSPEN_2_xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv',
                           header=T, sep="\t",as.matrix = F)
spen_ctl_ttt_2_snp=spen_ctl_ttt_2_snp[,cn(spen_ctl_ttt_2_snp)[grep('gene|xi_mean_fraction',cn(spen_ctl_ttt_2_snp))]]
spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siScr[spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siScr==0]<-0.01 
spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siSPEN_2[spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siSPEN_2==0]<-0.01
spen_ctl_ttt_2_snp$fc=spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siSPEN_2/spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siScr

rownames(spen_ctl_ttt_2_snp)=spen_ctl_ttt_2_snp$gene_snp

gene_cat_spen_2=gene_cat[gene_cat$gene%in%rn(spen_ctl_ttt_2_snp),]

spen_ctl_ttt_2_snp=spen_ctl_ttt_2_snp[gene_cat_spen_2$gene,]

spen_ctl_ttt_2_snp$xist_sensitive_category=gene_cat_spen_2$category

# spen_ctl_ttt_2_snp %>%
#   write.tsv(.,file="dataset_2/spen_ctl_ttt_2_xist_sensitive_gene_categories.tsv",row.names = F)

# levels(as.factor(spen_ctl_ttt_2_snp$category))

spen_ctl_ttt_1_snp=fastRead('dataset_2/pxgl_siScr_pxgl_siSPEN_1_xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv',
                            header=T, sep="\t",as.matrix = F)
spen_ctl_ttt_1_snp=spen_ctl_ttt_1_snp[,cn(spen_ctl_ttt_1_snp)[grep('gene|xi_mean_fraction',cn(spen_ctl_ttt_1_snp))]]
spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siScr[spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siScr==0]<-0.01 
spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siSPEN_1[spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siSPEN_1==0]<-0.01
spen_ctl_ttt_1_snp$fc=spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siSPEN_1/spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siScr

rownames(spen_ctl_ttt_1_snp)=spen_ctl_ttt_1_snp$gene_snp

gene_cat_spen_1=gene_cat[gene_cat$gene%in%rn(spen_ctl_ttt_1_snp),]

spen_ctl_ttt_1_snp=spen_ctl_ttt_1_snp[gene_cat_spen_1$gene,]

spen_ctl_ttt_1_snp$xist_sensitive_category=gene_cat_spen_1$category

# spen_ctl_ttt_1_snp %>%
#   write.tsv(.,file="dataset_2/spen_ctl_ttt_1_xist_sensitive_gene_categories.tsv",row.names = F)


# heatmap Xi fraction

naive_ctl_ttt_snp=fastRead('rna_seq_naive_control_rna_seq_naive_treatment_xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv',
                           header=T, sep="\t",as.matrix = F)
naive_ctl_ttt_snp=naive_ctl_ttt_snp[,cn(naive_ctl_ttt_snp)[grep('gene|xi_mean_fraction',cn(naive_ctl_ttt_snp))]]
naive_ctl_ttt_snp$xi_mean_fraction_rna_seq_naive_control[naive_ctl_ttt_snp$xi_mean_fraction_rna_seq_naive_control==0]<-0.01 
naive_ctl_ttt_snp$xi_mean_fraction_rna_seq_naive_treatment[naive_ctl_ttt_snp$xi_mean_fraction_rna_seq_naive_treatment==0]<-0.01
naive_ctl_ttt_snp$fc=naive_ctl_ttt_snp$xi_mean_fraction_rna_seq_naive_treatment/naive_ctl_ttt_snp$xi_mean_fraction_rna_seq_naive_control

rownames(naive_ctl_ttt_snp)=naive_ctl_ttt_snp$gene_snp
gene_cat_spen_1=gene_cat[gene_cat$gene%in%rn(naive_ctl_ttt_snp),]
naive_ctl_ttt_snp=naive_ctl_ttt_snp[gene_cat_spen_1$gene,]
naive_ctl_ttt_snp$xist_sensitive_category=gene_cat_spen_1$category


naive_wt_ko_snp=fastRead('rna_seq_naive_wt_rna_seq_naive_ko_xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv',
                         header=T, sep="\t",as.matrix = F)
naive_wt_ko_snp=naive_wt_ko_snp[,cn(naive_wt_ko_snp)[grep('gene|xi_mean_fraction',cn(naive_wt_ko_snp))]]
naive_wt_ko_snp$xi_mean_fraction_rna_seq_naive_wt[naive_wt_ko_snp$xi_mean_fraction_rna_seq_naive_wt==0]<-0.01 
naive_wt_ko_snp$xi_mean_fraction_rna_seq_naive_ko[naive_wt_ko_snp$xi_mean_fraction_rna_seq_naive_ko==0]<-0.01
naive_wt_ko_snp$fc=naive_wt_ko_snp$xi_mean_fraction_rna_seq_naive_ko/naive_wt_ko_snp$xi_mean_fraction_rna_seq_naive_wt

rownames(naive_wt_ko_snp)=naive_wt_ko_snp$gene_snp
gene_cat_spen_1=gene_cat[gene_cat$gene%in%rn(naive_wt_ko_snp),]
naive_wt_ko_snp=naive_wt_ko_snp[gene_cat_spen_1$gene,]
naive_wt_ko_snp$xist_sensitive_category=gene_cat_spen_1$category


list(naive_ctl_ttt_snp,naive_wt_ko_snp,spen_ctl_ttt_1_snp,spen_ctl_ttt_2_snp) %>%
  reduce(merge, by = "gene_snp",all=F) -> xi_fraction_merge

xi_fraction_merge=xi_fraction_merge[,grepl("xi_mean_fraction|xist_sensitive_category.x|gene_snp",cn(xi_fraction_merge))]

rownames(xi_fraction_merge)=xi_fraction_merge$gene_snp
xi_fraction_merge$gene_snp=NULL

xi_fraction_merge$xi_mean_fraction_pxgl_siScr.y=NULL
xi_fraction_merge$xist_sensitive_category.x=NULL
xi_fraction_merge$xist_sensitive_category.x.1=NULL


colnames(xi_fraction_merge)=c("naive_ctl","naive_xist_kd","naive_wt","naive_xist_ko","pxgl_siScr","pxgl_siSPEN_1","pxgl_siSPEN_2")

gene_cat_htmp=gene_cat[gene_cat$gene%in%rn(xi_fraction_merge),]

xi_fraction_merge=xi_fraction_merge[gene_cat_htmp$gene,]

pdf(file = "dataset_2/spen_xi_fraction_heatmap_gene_categories.pdf",width = 8,height = 10)

Heatmap(xi_fraction_merge, cluster_rows = F,cluster_columns = T, split = gene_cat_htmp$category, border = T, name="Xi fraction")

dev.off()



