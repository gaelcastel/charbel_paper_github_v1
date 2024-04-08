# vcf from bcftools already on SNP in H9 CHR X !!!!! BUT CAUTION !!! :not only heterozygous SNP !!! so I applied detect_bi on informative_H9_wgs.vcf
# Could be relevant to re-run bcftools intersect on this vcf -> bed only informative SNP heterozygous chrX to re-generate samples.vcf

setwd('C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/')

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

spen_ctl_ttt_2_snp=fastRead('dataset_2/pxgl_siScr_pxgl_siSPEN_2_xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv',
                           header=T, sep="\t",as.matrix = F)
spen_ctl_ttt_2_snp=spen_ctl_ttt_2_snp[,cn(spen_ctl_ttt_2_snp)[grep('gene|xi_mean_fraction',cn(spen_ctl_ttt_2_snp))]]
spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siScr[spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siScr==0]<-0.01 
spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siSPEN_2[spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siSPEN_2==0]<-0.01
spen_ctl_ttt_2_snp$fc=spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siSPEN_2/spen_ctl_ttt_2_snp$xi_mean_fraction_pxgl_siScr


# !!!! Ex: bi_bi_change -> bi_bi_up VS bi_bi_down

spen_ctl_ttt_2_snp$category=spen_ctl_ttt_2_snp$gene_snp # just for initialization
spen_ctl_ttt_2_snp$category<-apply(spen_ctl_ttt_2_snp,1,function(x){
  
  if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])<=0.25 & data.table::between(as.numeric(x['fc']),1/1.2,1.2)){return('mono_mono')}# I checked, no one > 0.75
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])<=0.25 & as.numeric(x['fc'])>=1.2){return('mono_mono_up')}
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])<=0.25 & as.numeric(x['fc'])<=1/1.2){return('mono_mono_down')}
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])>=0.25 & data.table::between(as.numeric(x['fc']),1/1.2,1.2)){return('bi_bi')}
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])>=0.25 & as.numeric(x['fc'])>=1.2){return('bi_bi_up')}
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])>=0.25 & as.numeric(x['fc'])<=1/1.2){return('bi_bi_down')}
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])>=0.25){return('mono_bi_up')}
  else if(as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siSPEN_2'])<=0.25){return('bi_mono_down')}
  
})

spen_ctl_ttt_2_snp$category[spen_ctl_ttt_2_snp$category%in%c("bi_bi")]<-"constant"
spen_ctl_ttt_2_snp$category[spen_ctl_ttt_2_snp$category%in%c("bi_bi_up","mono_bi_up")]<-"up"

spen_ctl_ttt_2_snp %>%
  filter(category%in%c("constant","up")) %>%
  summarise(gene=gene_snp,category=category) %>%
  write.tsv(.,file="spen_ctl_ttt_2_gene_categories.tsv",row.names = F)

# levels(as.factor(spen_ctl_ttt_2_snp$category))

spen_ctl_ttt_1_snp=fastRead('dataset_2/pxgl_siScr_pxgl_siSPEN_1_xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv',
                         header=T, sep="\t",as.matrix = F)
spen_ctl_ttt_1_snp=spen_ctl_ttt_1_snp[,cn(spen_ctl_ttt_1_snp)[grep('gene|xi_mean_fraction',cn(spen_ctl_ttt_1_snp))]]
spen_ctl_ttt_1_snp$xi_mean_fraction_rna_seq_naive_wt[spen_ctl_ttt_1_snp$xi_mean_fraction_rna_seq_naive_wt==0]<-0.01 
spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siScr[spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siScr==0]<-0.01
spen_ctl_ttt_1_snp$fc=spen_ctl_ttt_1_snp$xi_mean_fraction_pxgl_siScr/spen_ctl_ttt_1_snp$xi_mean_fraction_rna_seq_naive_wt

spen_ctl_ttt_1_snp$category=spen_ctl_ttt_1_snp$gene_snp # just for initialization
spen_ctl_ttt_1_snp$category<-apply(spen_ctl_ttt_1_snp,1,function(x){
  
  if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & data.table::between(as.numeric(x['fc']),1/1.2,1.2)){return('mono_mono')}# I checked, no one > 0.75
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & as.numeric(x['fc'])>=1.2){return('mono_mono_up')}
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25 & as.numeric(x['fc'])<=1/1.2){return('mono_mono_down')}
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & data.table::between(as.numeric(x['fc']),1/1.2,1.2)){return('bi_bi')}
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & as.numeric(x['fc'])>=1.2){return('bi_bi_up')}
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25 & as.numeric(x['fc'])<=1/1.2){return('bi_bi_down')}
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])<=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])>=0.25){return('mono_bi_up')}
  else if(as.numeric(x['xi_mean_fraction_rna_seq_naive_wt'])>=0.25 & as.numeric(x['xi_mean_fraction_pxgl_siScr'])<=0.25){return('bi_mono_down')}
  
})

spen_ctl_ttt_1_snp$category[spen_ctl_ttt_1_snp$category%in%c("bi_bi")]<-"constant"
spen_ctl_ttt_1_snp$category[spen_ctl_ttt_1_snp$category%in%c("bi_bi_up","mono_bi_up")]<-"up"

spen_ctl_ttt_1_snp %>%
  filter(category%in%c("constant","up")) %>%
  summarise(gene=gene_snp,category=category) %>%
  write.tsv(.,file="spen_ctl_ttt_1_gene_categories.tsv",row.names = F)
