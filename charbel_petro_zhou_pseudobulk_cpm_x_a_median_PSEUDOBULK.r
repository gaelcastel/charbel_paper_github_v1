#!/shared/#!/sh#!/sh#!/shre/minred/ifbstor1/software/miniconda/envs/r-4.2.1/bin/Rscript

setwd("D:/rougeulle_lab_ordi_sept_2023/charbel_paper/rna/charbel_paper_ifb/embryo_pseudobulk/")

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

library(dplyr)
library(ggplot2)
library(stringr)
library(sva)
library(ggpubr)
library(rstatix)
library(SingleCellExperiment)
library(preprocessCore)
library(lemon)
library(devtools)
library("Matrix.utils")

timeNow<-function(x=NULL, y=NULL, m=NULL, d=NULL , h=NULL, min=NULL,ymd_h_m=NULL){x=date();
date=unlist(str_split(pattern = "-", Sys.Date()))
time=unlist(str_split(pattern = ":",unlist(str_split(pattern=" ",Sys.time()))[2]))
return(paste(c(date,time), collapse = "_"))
}

gtf=read.table("D:/rougeulle_lab_ordi_sept_2023/charbel_paper/rna/charbel_paper_ifb/Homo_sapiens.GRCh38.90.gtf", header=F, sep="\t")
colnames(gtf)=c('chrom','source','feature','start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute')

gtf=gtf[!grepl('_dup',gtf$attribute),]


bi_bi_up=read.table("bi_bi_up.gtf", header=F, sep="\t")
colnames(bi_bi_up)=c('chrom','source','feature','start',
                     'end',
                     'score',
                     'strand',
                     'frame',
                     'attribute')

mono_bi_up=read.table("mono_bi_up.gtf", header=F, sep="\t")
colnames(mono_bi_up)=c('chrom','source','feature','start',
                       'end',
                       'score',
                       'strand',
                       'frame',
                       'attribute')

bi_bi=read.table("bi_bi.gtf", header=F, sep="\t")
colnames(bi_bi)=c('chrom','source','feature','start',
                  'end',
                  'score',
                  'strand',
                  'frame',
                  'attribute')

up_gtf=rbind(bi_bi_up,
             mono_bi_up) # I checked, these genes are common to wt-ko/ctl-ttt, cf file: C:\Users\gael\charbel_paper\rna\charbel_paper_ifb\embryo_pseudobulk\charbel_petro_zhou_v2.r
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

gene_chr=read.table("../hg38_gene_chr_list.txt")
gene_chr$V2=NULL
colnames(gene_chr)=c("gene","chr")

petro_zhou_sample_annot=fastRead("petro_zhou_sample_annot_cell_type_sex.tsv", header=T, sep="\t",as.matrix = F)

# PETRO ZHOU RAW COUNTS
load("2022_06_20_14_03_petro_total_raw_counts_not_corr.Rdata") #petro_total_raw_counts
load("2022_06_20_14_22_zhou_total_raw_counts_not_corr.Rdata") #zhou_total_raw_counts


# all(rn(assay(petro_total_raw_counts))==rn(assay(zhou_total_raw_counts)))
petro_total_raw=assay(petro_total_raw_counts)
zhou_total_raw=assay(zhou_total_raw_counts)

petro_zhou_raw=cbind(petro_total_raw,zhou_total_raw)

petro_total_raw=NULL
zhou_total_raw=NULL

colnames(petro_zhou_raw)=cn(petro_zhou_raw)
rownames(petro_zhou_raw)=rn(petro_zhou_raw)

petro_zhou_raw=as.data.frame(petro_zhou_raw)


petro_zhou_sample_annot=petro_zhou_sample_annot[rn(petro_zhou_sample_annot)%in%cn(petro_zhou_raw) & petro_zhou_sample_annot$sex=="female",]

petro_zhou_raw=petro_zhou_raw[,cn(petro_zhou_raw)%in%rn(petro_zhou_sample_annot)]

petro_zhou_raw$gene=sub("_ENSG.*","",rownames(petro_zhou_raw))

petro_zhou_raw %>%
  group_by(gene) %>%
  filter(n() > 1) %>%
  summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>%
  ungroup() %>%
  as.data.frame() -> petro_zhou_raw_duplicates

petro_zhou_raw %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  as.data.frame() -> petro_zhou_raw_unique

petro_zhou_raw=rbind(petro_zhou_raw_unique,petro_zhou_raw_duplicates)

rownames(petro_zhou_raw)=petro_zhou_raw$gene

petro_zhou_raw$gene=NULL

# all(cn(petro_zhou_raw)==rn(petro_zhou_sample_annot))

petro_zhou_sample_annot$cell_type_day_dataset=paste(sep="_",petro_zhou_sample_annot$finalClusters,petro_zhou_sample_annot$day,petro_zhou_sample_annot$dataset)

petro_zhou_sample_annot$cell_type_day_dataset=sub('_petro.*','',petro_zhou_sample_annot$cell_type_day_dataset)
petro_zhou_sample_annot$cell_type_day_dataset=sub('_zhou.*','',petro_zhou_sample_annot$cell_type_day_dataset)
petro_zhou_sample_annot$cell_type_day_dataset[petro_zhou_sample_annot$cell_type_day_dataset=="Morula_5arly"]<-"Morula_5"
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,".early", "")
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,".late|.medium|.medium1|.medium3", "")
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,"medium_|late_|Pre.", "")
petro_zhou_sample_annot$cell_type_day_dataset=sub('.*early_','',petro_zhou_sample_annot$cell_type_day_dataset)
petro_zhou_sample_annot$cell_type_day_dataset=sub('arly.*','',petro_zhou_sample_annot$cell_type_day_dataset)
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,"TB", "TE")
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,"TE1|TE2|TE3", "TE")
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,"EightCells", "8C")
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,"B1_B2", "BLASTO")
petro_zhou_sample_annot$cell_type_day_dataset=str_replace(petro_zhou_sample_annot$cell_type_day_dataset,"Morula", "MORULA")



petro_zhou_sample_annot_female_filtered=petro_zhou_sample_annot[petro_zhou_sample_annot$sex=="female" & !grepl("apoptosis|ys|BLASTO_6",petro_zhou_sample_annot$cell_type_day_dataset),] # subset on female embryos for chrX related questions

petro_zhou_raw_female=petro_zhou_raw[,rn(petro_zhou_sample_annot_female_filtered)]

petro_zhou_sample_annot_female_filtered$day[petro_zhou_sample_annot_female_filtered$day=="5arly"]<-"5"
petro_zhou_sample_annot_female_filtered$day=factor(petro_zhou_sample_annot_female_filtered$day,levels=c("3","4","5","6","7","8","10","12"))

petro_zhou_sample_annot_female_filtered$cell_type_day_dataset=factor(petro_zhou_sample_annot_female_filtered$cell_type_day_dataset,
                                                                     levels=c("8C_3","MORULA_4","MORULA_5",
                                                                              "BLASTO_4","BLASTO_5",
                                                                              "EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                                                              "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                                                              "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12",
                                                                              "ST_8","ST_10","ST_12","EVT_6","EVT_10","EVT_12"
                                                                     ))


petro_zhou_raw_female_pseudobulk <- aggregate.Matrix(t(petro_zhou_raw_female), 
                                                     groupings = petro_zhou_sample_annot_female_filtered$cell_type_day_dataset, fun = "sum") 


petro_zhou_raw_female_pseudobulk=as.data.frame(as.matrix(petro_zhou_raw_female_pseudobulk))

petro_zhou_raw_female_pseudobulk=t(petro_zhou_raw_female_pseudobulk)
petro_zhou_raw_female_pseudobulk=CPM(petro_zhou_raw_female_pseudobulk)
petro_zhou_raw_female_pseudobulk=as.data.frame(petro_zhou_raw_female_pseudobulk)

petro_zhou_raw_female_pseudobulk$gene=sub("_ENSG.*","",rownames(petro_zhou_raw_female_pseudobulk))

petro_zhou_raw_female_pseudobulk_x_genes=petro_zhou_raw_female_pseudobulk[petro_zhou_raw_female_pseudobulk$gene%in%gene_chr$gene[gene_chr$chr=="chrX"],]
# petro_zhou_raw_female_pseudobulk_x_genes=petro_zhou_raw_female_pseudobulk_x_genes[rowSums(petro_zhou_raw_female_pseudobulk_x_genes)>5,] # discard low-expressed genes

petro_zhou_raw_female_pseudobulk_a_genes=petro_zhou_raw_female_pseudobulk[petro_zhou_raw_female_pseudobulk$gene%in%gene_chr$gene[!(gene_chr$chr%in%c("chrX","chrY"))],]
# petro_zhou_raw_female_pseudobulk_a_genes=petro_zhou_raw_female_pseudobulk_a_genes[rowMeans(petro_zhou_raw_female_pseudobulk_a_genes)>5,] # discard low-expressed genes

petro_zhou_raw_female_pseudobulk_x_genes$gene=NULL
petro_zhou_raw_female_pseudobulk_a_genes$gene=NULL

stat_x=apply(petro_zhou_raw_female_pseudobulk_x_genes,2,function(x){median(as.numeric(x)[as.numeric(x)>0])}) # > 0 keep only expressed genes, otherwise 0 inflation, is absence of measure, not a measure
stat_a=apply(petro_zhou_raw_female_pseudobulk_a_genes,2,function(x){median(as.numeric(x)[as.numeric(x)>0])}) # > 0 keep only expressed genes, otherwise 0 inflation, is absence of measure, not a measure

petro_zhou_raw_female_pseudobulk %>%
  tidyr::pivot_longer(!c(gene), names_to = "cell_type_day_dataset", values_to = "count") -> petro_zhou_raw_female_tibble

# petro_zhou_raw_female_tibble %>%
#   group_by(cell_type_day_dataset,gene) %>%
#   summarise(count=sum(count)) -> petro_zhou_raw_female_tibble

# petro_zhou_raw_female_pseudobulk=NULL

petro_zhou_raw_female_tibble$cell_type=sub('_[^_]*$', '',petro_zhou_raw_female_tibble$cell_type_day_dataset)

petro_zhou_raw_female_tibble %>%
  filter(gene%in%constant_gtf$gene) %>%
  mutate(category="constant") -> petro_zhou_raw_female_constant

petro_zhou_raw_female_tibble %>%
  filter(gene%in%up_gtf$gene) %>%
  mutate(category="up") -> petro_zhou_raw_female_up

petro_zhou_raw_female_tibble=NULL

petro_zhou_raw_female_tibble_plot=rbind(petro_zhou_raw_female_constant,
                                   petro_zhou_raw_female_up)



petro_zhou_raw_female_tibble_plot %>%
  mutate(log_expr=count) -> petro_zhou_raw_female_tibble_plot

petro_zhou_raw_female_tibble_plot$cell_type=sub("_.*","",petro_zhou_raw_female_tibble_plot$cell_type_day_dataset)


petro_zhou_raw_female_tibble_plot%>%
  filter(cell_type_day_dataset%in%c("8C_3","MORULA_4","BLASTO_4","BLASTO_5","EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                    "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                    "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12")) -> petro_zhou_raw_female_tibble_plot


petro_zhou_raw_female_tibble_plot$cell_type_day_dataset=factor(petro_zhou_raw_female_tibble_plot$cell_type_day_dataset,
                                                                 levels=c("8C_3","MORULA_4","BLASTO_4","BLASTO_5","EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                                                          "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                                                          "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12"))

petro_zhou_raw_female_tibble_plot %>%  
  group_by(category)%>%
  mutate(median_cat=median(count)) %>%
  group_by(cell_type_day_dataset,category)%>%
  mutate(sd=sd(count), median=median(count)) -> petro_zhou_raw_female_tibble_plot_df

stat_x=stat_x[names(stat_x)%in%unique(petro_zhou_raw_female_tibble_plot_df$cell_type_day_dataset)]
stat_x=stat_x[levels(petro_zhou_raw_female_tibble_plot_df$cell_type_day_dataset)]

stat_a=stat_a[names(stat_a)%in%unique(petro_zhou_raw_female_tibble_plot_df$cell_type_day_dataset)]
stat_a=stat_a[levels(petro_zhou_raw_female_tibble_plot_df$cell_type_day_dataset)]

petro_zhou_raw_female_tibble_plot_df$x_median=apply(petro_zhou_raw_female_tibble_plot_df,1,function(x){stat_x[x["cell_type_day_dataset"]]})
petro_zhou_raw_female_tibble_plot_df$a_median=apply(petro_zhou_raw_female_tibble_plot_df,1,function(x){stat_a[x["cell_type_day_dataset"]]})

petro_zhou_raw_female_tibble_plot_df$group=rep("A",nrow(petro_zhou_raw_female_tibble_plot_df))

petro_zhou_raw_female_tibble_plot_df%>%
  ggplot(aes(x=cell_type_day_dataset)) +
  geom_boxplot(aes(fill=category, y=count),lwd=0.15,position=position_dodge(width = 0.8),outlier.shape = NA)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  geom_pointline(aes(y=x_median,group=group), colour="green")+
  geom_pointline(aes(y=a_median,group=group), colour="red")+
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA),breaks=c(0,10,100,1000))+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("CPM gene expression") +
  ggtitle("Gene categories upon XIST depletion in Naive hESCs\nExpression in female embryos") +
  # xlab("Cell Type")+
  geom_blank(aes(y = 0)) + #y_axis start at 0
  # scale_fill_manual("Cell Type :", values = colors)+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# ggsave(paste(sep="_","gene_categories_in_embryo_female_quantile_cpm_norw_petro_zhou_x_a_medians_PSEUDOBULK_v2.pdf"),width=8,height=6)


petro_zhou_raw_female_tibble_plot_df%>%  
  ggplot(aes(x=cell_type_day_dataset)) +
  geom_bar(aes(y=count,fill=category),stat = 'summary',fun="median",position = 'dodge')+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  geom_pointline(aes(y=x_median,group=group), colour="green")+
  geom_pointline(aes(y=a_median,group=group), colour="red")+
  # facet_wrap(~branch, nrow=1, ncol=4, scales = "fixed")+
  # geom_point(aes(y=count,fill=category))+
  # geom_text(aes(y=count,label=gene))+
  # scale_y_continuous(trans="pseudo_log",limits = c(0, NA),breaks=c(0,10,100,1000))+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("CPM gene expression") +
  ggtitle("Gene categories upon XIST depletion in Naive hESCs\nExpression in female embryos (median)") +
  # xlab("Cell Type")+
  geom_blank(aes(y = 0)) + #y_axis start at 0
  # scale_fill_manual("Cell Type :", values = colors)+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# res.stat <- petro_zhou_raw_female_tibble_plot_df_new_categories %>% 
#   group_by(cell_type_day_dataset)%>%
#   t_test(data = ., count ~ category,paired = F) %>% 
#   adjust_pvalue(method = "bonferroni") %>%
#   # add_significance("p.adj") %>%
#   add_significance("p") %>%
#   ungroup()%>%
#   add_xy_position(x = "cell_type_day_dataset", fun = "max",formula = count ~ category,dodge = 0.9) %>%
#   ungroup()
# 
# print(p + 
#         stat_pvalue_manual(
#           res.stat, label = "p.adj.signif", y.position=c(250,245,240), bracket.size = 0.5, tip.length = 0.001
#         )
#       # +
#       #   stat_pvalue_manual(
#       #     res.stat, label = "p", tip.length = 0.01
#       #   )
#       )

ggsave(paste(sep="_","BARPLOTS_gene_categories_in_embryo_female_quantile_cpm_norw_petro_zhou_x_a_medians_old_categories_PSEUDOBULK.pdf"),width=8*1.33,height=6)






petro_zhou_raw_female_tibble_plot_df_new_categories=petro_zhou_raw_female_tibble_plot_df

up_new=c("SH3BGRL","DMD","NBDY","IGSF1","INTS6L","RTL8B","DANT1","DANT2","RTL8A","HDAC8","GSPT2")

petro_zhou_raw_female_tibble_plot_df_new_categories$category[petro_zhou_raw_female_tibble_plot_df_new_categories$gene%in%up_new]<-"up_new"
petro_zhou_raw_female_tibble_plot_df_new_categories$category[petro_zhou_raw_female_tibble_plot_df_new_categories$category=="up"]<-"up_old"

petro_zhou_raw_female_tibble_plot_df_new_categories%>%
  group_by(category)%>%
  mutate(median_cat=median(count)) %>%
  group_by(cell_type_day_dataset,category)%>%
  mutate(sd=sd(count), median=median(count)) -> petro_zhou_raw_female_tibble_plot_df_new_categories

petro_zhou_raw_female_tibble_plot_df_new_categories%>%  
  ggplot(aes(x=cell_type_day_dataset)) +
  geom_boxplot(aes(fill=category, y=count),lwd=0.15,position=position_dodge(width = 0.8),outlier.shape = NA)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df_new_categories$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df_new_categories$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df_new_categories$median_cat)[3], color='darkgreen',linetype='dashed',size=0.75)+
  geom_pointline(aes(y=x_median,group=group), colour="green")+
  geom_pointline(aes(y=a_median,group=group), colour="red")+
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA),breaks=c(0,10,100,1000))+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("CPM gene expression") +
  ggtitle("Gene categories upon XIST depletion in Naive hESCs\nExpression in female embryos") +
  # xlab("Cell Type")+
  geom_blank(aes(y = 0)) + #y_axis start at 0
  # scale_fill_manual("Cell Type :", values = colors)+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# ggsave(paste(sep="_","gene_categories_in_embryo_female_quantile_cpm_norw_petro_zhou_x_a_medians_new_categories_PSEUDOBULK_v2.pdf"),width=8*1.33,height=6)

petro_zhou_raw_female_tibble_plot_df_new_categories$branch=petro_zhou_raw_female_tibble_plot_df_new_categories$cell_type
petro_zhou_raw_female_tibble_plot_df_new_categories$branch[petro_zhou_raw_female_tibble_plot_df_new_categories$branch%in%c("8C",
                                                                                                                           "MORULA",
                                                                                                                           "BLASTO")]<-"Unspecified"

petro_zhou_raw_female_tibble_plot_df_new_categories$branch=factor(petro_zhou_raw_female_tibble_plot_df_new_categories$branch,
                                                                  levels=c("Unspecified","EPI","PrE","TE"))

petro_zhou_raw_female_tibble_plot_df_new_categories%>%  
  ggplot(aes(x=cell_type_day_dataset)) +
  geom_bar(aes(y=count,fill=category),stat = 'summary',fun="median",position = 'dodge')+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df_new_categories$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df_new_categories$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df_new_categories$median_cat)[3], color='darkgreen',linetype='dashed',size=0.75)+
  geom_pointline(aes(y=x_median,group=group), colour="green")+
  geom_pointline(aes(y=a_median,group=group), colour="red")+
  # facet_wrap(~branch, nrow=1, ncol=4, scales = "fixed")+
  # geom_point(aes(y=count,fill=category))+
  # geom_text(aes(y=count,label=gene))+
  # scale_y_continuous(trans="pseudo_log",limits = c(0, NA),breaks=c(0,10,100,1000))+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("CPM gene expression") +
  ggtitle("Gene categories upon XIST depletion in Naive hESCs\nExpression in female embryos (median)") +
  # xlab("Cell Type")+
  geom_blank(aes(y = 0)) + #y_axis start at 0
  # scale_fill_manual("Cell Type :", values = colors)+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# res.stat <- petro_zhou_raw_female_tibble_plot_df_new_categories %>% 
#   group_by(cell_type_day_dataset)%>%
#   t_test(data = ., count ~ category,paired = F) %>% 
#   adjust_pvalue(method = "bonferroni") %>%
#   # add_significance("p.adj") %>%
#   add_significance("p") %>%
#   ungroup()%>%
#   add_xy_position(x = "cell_type_day_dataset", fun = "max",formula = count ~ category,dodge = 0.9) %>%
#   ungroup()
# 
# print(p + 
#         stat_pvalue_manual(
#           res.stat, label = "p.adj.signif", y.position=c(250,245,240), bracket.size = 0.5, tip.length = 0.001
#         )
#       # +
#       #   stat_pvalue_manual(
#       #     res.stat, label = "p", tip.length = 0.01
#       #   )
#       )

# ggsave(paste(sep="_","BARPLOTS_gene_categories_in_embryo_female_quantile_cpm_norw_petro_zhou_x_a_medians_new_categories_PSEUDOBULK.pdf"),width=8*1.33,height=6)







# NO overall higher global expression in Zhou compared to Petropoulos at day 6

petro_zhou_sample_annot_female_filtered %>%
  filter(day==6) %>% # Petro and Zhou cells at day 6, so apparently NO "overall" higher expression due to Zhou vs Petro (=NO technical bias)
  View()
