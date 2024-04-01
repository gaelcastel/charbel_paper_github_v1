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
library(rstatix)

gene_chr=read.table("hg38_gene_chr_list.txt")
gene_chr$V2=NULL
colnames(gene_chr)=c("gene","chr")

raw_count_files=list.files(path = 'raw_counts/',pattern = "_raw_counts.txt")
sample_annot=fastRead("./sample_annot_all.txt", header=T, sep=",",as.matrix = F)
sample_annot$sample_id=rn(sample_annot)
sample_annot$cell_type_condition=paste(sep="_",sample_annot$cell_type,sample_annot$condition)

sample_annot$sample_category=paste(sep="_",sample_annot$cell_type_condition,
                                   sample_annot$dox_ko,
                                   sample_annot$dataset)

sample_annot$sample_category=str_replace(sample_annot$sample_category,"_control","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_rna","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_treatment","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_seq","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d1269","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_bis","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d165","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d1271","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_dvallot","_vallot")

sample_annot$sample_category=factor(sample_annot$sample_category,
                                    levels = c("primed_WT",
                                               "primed_non_eroded_WT_vallot",
                                               "primed_eroded_WT_vallot",
                                               "primed_KO_XACT",
                                               "naive_cult_WT",
                                               "naive_cult_KO_XACT",
                                               "naive_WT",
                                               "naive_KO",
                                               "naive_C",
                                               "naive_T",
                                               "naive_C_gro",
                                               "naive_T_gro"))      

sample_annot$group=as.vector(sample_annot$sample_category)
sample_annot$group[sample_annot$group%in%c("primed_WT","primed_KO_XACT")]<-"group1"
sample_annot$group[sample_annot$group%in%c("naive_cult_WT","naive_cult_KO_XACT")]<-"group2"
sample_annot$group[sample_annot$group%in%c("primed_non_eroded_WT_vallot","primed_eroded_WT_vallot")]<-"group3"
sample_annot$group[sample_annot$group%in%c("naive_WT","naive_KO")]<-"group4"
sample_annot$group[sample_annot$group%in%c("naive_C","naive_T")]<-"group5"
sample_annot$group[sample_annot$group%in%c("naive_C_gro","naive_T_gro")]<-"group6"


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

# NB Céline : XACT KO dans les naive cult n'a pas aussi bien marché que dans les primed, donc logique "fuite" d'expression ?

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
                       sample_category=sample_annot$sample_category,
                       condition=sample_annot$condition,
                       cell_type_condition=sample_annot$cell_type_condition,
                       dataset=sample_annot$dataset,
                       dox_ko=sample_annot$dox_ko,
                       name=sample_annot$Sample_Name,
                       dataset_cell_type_condition_dox_ko=paste(sep="_",sample_annot$cell_type_condition,sample_annot$dox_ko,sample_annot$dataset))

stat_4_plot %>%
  filter(dataset !="rna_seq_dvallot") %>%
  # filter(sample_category %in% c("naive_cult_WT","naive_cult_KO_XACT")) %>%
  group_by(sample_category) %>%
    mutate(sd=sd(x_a_ratio), mean=mean(x_a_ratio)) %>%
    ggplot() + 
  geom_bar(aes(x=group,y=x_a_ratio,fill=sample_category),stat = 'summary',fun="mean",position = 'dodge')+
  geom_point(aes(x=group,y=x_a_ratio,group=sample_category),position = position_dodge(width = 0.9))+
  # facet_wrap(~category, nrow=1, ncol=2)+
  geom_errorbar(aes(x=group,y=x_a_ratio,group=sample_category,ymin=mean-sd, ymax=mean+sd), position="dodge") +
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
    scale_y_continuous(breaks = seq(0,1.1, 0.1)) + ylim(0,1.1)-> p

res.stat <- stat_4_plot %>%
  filter(dataset !="rna_seq_dvallot") %>%
  # filter(sample_category %in% c("naive_cult_WT","naive_cult_KO_XACT")) %>%
  group_by(group) %>%
  t_test(data = ., x_a_ratio ~ sample_category,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "group", fun = "max",formula = x_a_ratio ~ sample_category,step.increase = 0.05,dodge = 0.9) %>%
  ungroup()

# pdf("xist_and_markers_htseq_gene_expression.pdf",width = 9,height = 6)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif", y.position = res.stat$y.position
        ) + 
        stat_pvalue_manual(
          res.stat, label = "p", y.position = res.stat$y.position+0.05
        ) 
)

# ggsave("naive_cult_wt_xact_ko_x_a_ratio_0_1.pdf",width = 6,height = 7)

    
