setwd('D:/rougeulle_lab_ordi_sept_2023/charbel_paper/rna/charbel_paper_ifb/')

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
# library(EnhancedVolcano)


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

mmat_df=mmat_df[,!grepl('1271|SRR',cn(mmat_df))] #remove GRO-Seq and Vallot samples

sample_annot=sample_annot[cn(mmat_df),]

colnames(mmat_df)=sample_annot$Sample_Name

mmat_mtx=as.matrix(mmat_df)

mmat_mtx_cpm=as.data.frame(CPM(mmat_mtx))

df_plot=mmat_mtx_cpm[c("CD24","DNMT3L","NANOG","XIST","SPEN","XACT"),]
df_plot$gene=rn(df_plot)

df_plot %>%
  tidyr::pivot_longer(!c(gene), names_to = "sample", values_to = "CPM") -> df_plot

df_plot$cell_type=unlist(lapply(str_split(pattern = "_",df_plot$sample),function(x){
  
  x[grep("PXGL|naive|primed",x)]
  
}))

df_plot$cell_type[df_plot$cell_type=="naive"]<-"naive_cult"
df_plot$cell_type[df_plot$cell_type=="PXGL"]<-"naive"

df_plot$condition=unlist(lapply(str_split(pattern = "_",df_plot$sample),function(x){
  
  x[x%in%c("WT","KO","C","T","H9")]
  
}))

df_plot$condition[df_plot$condition%in%c("C","WT")]<-"ctl"
df_plot$condition[df_plot$condition%in%c("T","KO")]<-"ttt"
df_plot$condition[df_plot$condition%in%c("H9")]<-"ctl"

df_plot$cell_type_condition=paste(sep="_",df_plot$cell_type,df_plot$condition)

df_plot %>%
  filter(cell_type%in%c("primed","naive","naive_cult") & condition=="ctl") %>%
  ggplot() +
  geom_boxplot(aes(x=cell_type, y=CPM, fill=cell_type),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~gene, nrow=3, ncol=3,scales = 'free')+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  # stat_boxplot(aes(x=category, y=tpm),geom='errorbar', linetype=1, width=0.8, size=0.25)+
  scale_y_continuous(limits = c(0, NA))+
  # stat_compare_means(aes(x=cell_type, y=CPM),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("Marker gene expression Charbel revision") -> p

print(p)

ggsave("D:/charbel_december_2023/stem_cell_expression.pdf",width=8,height = 7)

res.stat <- df_plot %>%
  filter(cell_type%in%c("primed","naive","naive_cult") & condition=="ctl") %>%
  group_by(gene)%>%
  t_test(data = ., CPM ~ cell_type_condition,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "cell_type_condition", fun = "max",formula = CPM ~ cell_type_condition,step.increase = 0.05,dodge = 0.9) %>%
  ungroup()

# pdf("xist_and_markers_htseq_gene_expression.pdf",width = 9,height = 6)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif"
        )
      
)


ggsave("D:/charbel_december_2023/stem_cell_expression_stat.pdf",width=8,height = 7)







##############

df_plot %>%
  filter(cell_type_condition%in%c("primed_ctl","naive_ctl","naive_ttt")) %>%
  ggplot() +
  geom_boxplot(aes(x=cell_type_condition, y=CPM, fill=cell_type_condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~gene, nrow=2, ncol=3,scales = 'free')+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  # stat_boxplot(aes(x=category, y=tpm),geom='errorbar', linetype=1, width=0.8, size=0.25)+
  scale_y_continuous(limits = c(0, NA))+
  # stat_compare_means(aes(x=cell_type_condition, y=CPM, group=cell_type_condition),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("XIST gene expression - htseq-count all condtions") -> p

res.stat <- df_plot %>%
  filter(cell_type_condition%in%c("primed_ctl","naive_ctl","naive_ttt")) %>%
  group_by(gene)%>%
  t_test(data = ., CPM ~ cell_type_condition,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "cell_type_condition", fun = "max",formula = CPM ~ cell_type_condition,step.increase = 0.05,dodge = 0.9) %>%
  ungroup()

# pdf("xist_and_markers_htseq_gene_expression.pdf",width = 9,height = 6)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif"
        )
      
)

# dev.off()
