#!/shared/#!/sh#!/sh#!/shre/minred/ifbstor1/software/miniconda/envs/r-4.2.1/bin/Rscript


setwd("D:/rougeulle_lab_ordi_sept_2023/charbel_paper/rna/charbel_paper_ifb/embryo_pseudobulk/")

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

library(dplyr)
library(ggplot2)
library(stringr)
library(sva)
library(ggpubr)
library(rstatix)
library(remotes)
# install_github("https://github.com/bmbolstad/preprocessCore")
library(preprocessCore)

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

petro_zhou_sample_annot=fastRead("petro_zhou_sample_annot_cell_type_sex.tsv", header=T, sep="\t",as.matrix = F)

# PETRO ZHOU RAW COUNTS
load("2022_06_20_14_03_petro_total_raw_counts_not_corr.Rdata") #petro_total_raw_counts
load("2022_06_20_14_22_zhou_total_raw_counts_not_corr.Rdata") #zhou_total_raw_counts

# all(rn(assay(petro_total_raw_counts))==rn(assay(zhou_total_raw_counts)))
petro_total_raw_counts=assay(petro_total_raw_counts)
zhou_total_raw_counts=assay(zhou_total_raw_counts)

petro_zhou_total_raw_counts=cbind(petro_total_raw_counts,
                                  zhou_total_raw_counts)

petro_zhou_norm=preprocessCore::normalize.quantiles(log2(CPM(petro_zhou_total_raw_counts)+1))

colnames(petro_zhou_norm)=cn(petro_zhou_total_raw_counts)
rownames(petro_zhou_norm)=rn(petro_zhou_total_raw_counts)

petro_zhou_norm=as.data.frame(petro_zhou_norm)

petro_zhou_norm$gene=sub("_ENSG.*","",rownames(petro_zhou_norm))

selected_genes=c("XIST", "XACT", "SPEN", "NANOG", "GATA3")

petro_zhou_norm=petro_zhou_norm[petro_zhou_norm$gene%in%selected_genes,]

petro_zhou_norm$gene=NULL

petro_zhou_sample_annot=petro_zhou_sample_annot[rn(petro_zhou_sample_annot)%in%cn(petro_zhou_norm),]
petro_zhou_norm=petro_zhou_norm[,cn(petro_zhou_norm)%in%rn(petro_zhou_sample_annot)]

petro_zhou_norm=petro_zhou_norm[,rn(petro_zhou_sample_annot)]

# all(cn(petro_zhou_norm)==rn(petro_zhou_sample_annot))

t_petro_zhou_norm=as.data.frame(t(petro_zhou_norm))
# all(rn(t_petro_zhou_norm)==rn(petro_zhou_sample_annot))

t_petro_zhou_norm$cell_type=petro_zhou_sample_annot$finalClusters
t_petro_zhou_norm$dataset=petro_zhou_sample_annot$dataset
t_petro_zhou_norm$day=petro_zhou_sample_annot$day
t_petro_zhou_norm$sex=petro_zhou_sample_annot$sex
t_petro_zhou_norm$cell_type_day_dataset=paste(sep="_",t_petro_zhou_norm$cell_type,t_petro_zhou_norm$day,t_petro_zhou_norm$dataset)

t_petro_zhou_norm$cell_type_day_dataset=sub('_petro.*','',t_petro_zhou_norm$cell_type_day_dataset)
t_petro_zhou_norm$cell_type_day_dataset=sub('_zhou.*','',t_petro_zhou_norm$cell_type_day_dataset)
t_petro_zhou_norm$cell_type_day_dataset[t_petro_zhou_norm$cell_type_day_dataset=="Morula_5arly"]<-"Morula_5"
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,".early", "")
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,".late|.medium|.medium1|.medium3", "")
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,"medium_|late_|Pre.", "")
t_petro_zhou_norm$cell_type_day_dataset=sub('.*early_','',t_petro_zhou_norm$cell_type_day_dataset)
t_petro_zhou_norm$cell_type_day_dataset=sub('arly.*','',t_petro_zhou_norm$cell_type_day_dataset)
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,"TB", "TE")
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,"TE1|TE2|TE3", "TE")
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,"EightCells", "8C")
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,"B1_B2", "BLASTO")
t_petro_zhou_norm$cell_type_day_dataset=str_replace(t_petro_zhou_norm$cell_type_day_dataset,"Morula", "MORULA")

# subset on female embryos for chrX related questions

t_petro_zhou_norm %>%
  filter(sex=="female") -> t_petro_zhou_norm_female

t_petro_zhou_norm_female$sample=rn(t_petro_zhou_norm_female)

t_petro_zhou_norm_female %>%
  tidyr::pivot_longer(!c(sample,
                         cell_type,
                         dataset,
                         day,
                         sex,
                         cell_type_day_dataset), names_to = "gene_id", values_to = "count") -> petro_zhou_raw_female_tibble

petro_zhou_raw_female_tibble$gene=sub("_ENSG.*","",petro_zhou_raw_female_tibble$gene_id)

petro_zhou_raw_female_tibble %>%
  filter(!grepl("apoptosis|ys|BLASTO_6",cell_type_day_dataset)) -> petro_zhou_raw_female_tibble

petro_zhou_raw_female_tibble$day[petro_zhou_raw_female_tibble$day=="5arly"]<-"5"
petro_zhou_raw_female_tibble$day=factor(petro_zhou_raw_female_tibble$day,levels=c("3","4","5","6","7","8","10","12"))

petro_zhou_raw_female_tibble$cell_type_day_dataset=factor(petro_zhou_raw_female_tibble$cell_type_day_dataset,
                                                          levels=c("8C_3","MORULA_4","MORULA_5",
                                                                   "BLASTO_4","BLASTO_5",
                                                                   "EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                                                   "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                                                   "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12",
                                                                   "ST_8","ST_10","ST_12","EVT_6","EVT_10","EVT_12"
                                                          ))

petro_zhou_raw_female_tibble$cell_type=sub("_.*","",petro_zhou_raw_female_tibble$cell_type_day_dataset)

petro_zhou_raw_female_tibble %>%
  group_by(cell_type_day_dataset,gene) %>%
  filter(n()>=2) %>% # at least two cells per condition
  ungroup() %>%
  split(f = as.factor(.$gene)) -> petro_zhou_raw_female_tibble_category

violin_plots_callithrix_embryo_lineage_markers_constant=lapply(petro_zhou_raw_female_tibble_category,function(x){ggplot(x,aes(x=day,y=2^(count)-1,fill=cell_type_day_dataset))+
    geom_violin(scale = "width")+
    geom_boxplot(outlier.shape = NA,width=0.3,position=position_dodge(.9))+
    facet_wrap(~gene+cell_type, scales = "free_x",nrow = 2,ncol = 4) +
    guides(fill=guide_legend(ncol=1)) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
    # scale_y_continuous(limits = c(0, NA),trans = "pseudo_log")
  # stat_compare_means(label.y = 9.8,label="p",size=3, method="t.test", paired=FALSE)+
  # stat_compare_means(label.y = 10,label="p.signif",size=3, method="t.test", paired=FALSE)
})










# petro_zhou_raw_female_tibble %>%
#   group_by(cell_type_day_dataset)%>%
#   wilcox_test(count ~ category) %>% 
#   adjust_pvalue(method = 'bonferroni')

# pdf(paste(sep="_","gene_categories_in_embryo_female_norm_expr.pdf"),width=15,height=10)
petro_zhou_raw_female_tibble%>%
  # filter(cell_type_day_dataset!="MORULA_5")%>%
  ggplot(aes(x=cell_type_day_dataset, y=count,fill=category)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.shape = NA)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("Log norm gene expression") +
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
# scale_y_continuous(trans="pseudo_log")+

# stat_compare_means(vjust=1,label="p.signif",size=3, method="t.test", paired=F)

dev.off()



# petro_zhou_raw_female_tibble %>%
#   group_by(cell_type_day_dataset)%>%
#   wilcox_test(count ~ category) %>% 
#   adjust_pvalue(method = 'bonferroni')

petro_zhou_raw_female_tibble%>%
  filter(cell_type_day_dataset!="MORULA_5")%>%
  filter(cell_type_day_dataset!="8C_3")%>%
  group_by(category)%>%
  mutate(median_cat=median(count)) %>%
  group_by(cell_type_day_dataset,category)%>%
  mutate(sd=sd(count), median=median(count)) -> petro_zhou_raw_female_tibble_df

ggplot(petro_zhou_raw_female_tibble_df,aes(x=cell_type_day_dataset, y=count, fill=category)) + 
  # geom_errorbar(aes(ymin=median , ymax=median+sd),position = "dodge")+
  geom_bar(stat = "summary",fun='median', position = "dodge")+
  # geom_point(aes(x=cell_type_day_dataset, y=count, fill=category, group=category),size=0.25,alpha=0.8)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_df$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_df$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("Gene categories upon XIST depletion in Naive hESCs\nExpression in female embryos")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("Median norm gene expression") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
# scale_y_continuous(trans="pseudo_log",limits = c(0, NA))

# ggsave(paste(sep="_","gene_categories_in_embryo_female_median_norm_expr_barplots.pdf"),width=8,height=6)



petro_zhou_raw_female_tibble%>%
  filter(cell_type_day_dataset!="MORULA_5")%>%
  # filter(cell_type_day_dataset!="8C_3")%>%
  group_by(category)%>%
  mutate(median_cat=median(count)) %>%
  group_by(cell_type_day_dataset,category)%>%
  mutate(sd=sd(count), median=median(count)) -> petro_zhou_raw_female_tibble_df

petro_zhou_raw_female_tibble_df%>%
  ggplot(aes(x=cell_type_day_dataset, y=count,fill=category)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.shape = NA)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_df$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_df$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("Log norm gene expression") +
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
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))

# ggsave(paste(sep="_","gene_categories_in_embryo_female_norm_expr_boxplots_with_8c.pdf"),width=8,height=6)


# petro_zhou_raw_female_tibble %>%
#   filter(cell_type_day_dataset=="8C_3") %>%
#   write.tsv(file = "gene_categories.tsv")


###################################SUM_GENE_CATEGORY###################################

petro_zhou_sample_annot_sum_gene=petro_zhou_sample_annot

t_petro_zhou_total_raw_counts=as.data.frame(t(petro_zhou_total_raw_counts))

t_petro_zhou_total_raw_counts=t_petro_zhou_total_raw_counts[rn(t_petro_zhou_total_raw_counts)%in%rn(petro_zhou_sample_annot_sum_gene),]
petro_zhou_sample_annot_sum_gene=petro_zhou_sample_annot_sum_gene[rn(t_petro_zhou_total_raw_counts),]

# all(rn(t_petro_zhou_total_raw_counts)==rn(petro_zhou_sample_annot_sum_gene))

t_petro_zhou_total_raw_counts$cell_type=petro_zhou_sample_annot_sum_gene$finalClusters
t_petro_zhou_total_raw_counts$dataset=petro_zhou_sample_annot_sum_gene$dataset
t_petro_zhou_total_raw_counts$day=petro_zhou_sample_annot_sum_gene$day
t_petro_zhou_total_raw_counts$sex=petro_zhou_sample_annot_sum_gene$sex
t_petro_zhou_total_raw_counts$cell_type_day_dataset=paste(sep="_",t_petro_zhou_total_raw_counts$cell_type,t_petro_zhou_total_raw_counts$day,t_petro_zhou_total_raw_counts$dataset)

t_petro_zhou_total_raw_counts$cell_type_day_dataset=sub('_petro.*','',t_petro_zhou_total_raw_counts$cell_type_day_dataset)
t_petro_zhou_total_raw_counts$cell_type_day_dataset=sub('_zhou.*','',t_petro_zhou_total_raw_counts$cell_type_day_dataset)
t_petro_zhou_total_raw_counts$cell_type_day_dataset[t_petro_zhou_total_raw_counts$cell_type_day_dataset=="Morula_5arly"]<-"Morula_5"
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,".early", "")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,".late|.medium|.medium1|.medium3", "")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,"medium_|late_|Pre.", "")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=sub('.*early_','',t_petro_zhou_total_raw_counts$cell_type_day_dataset)
t_petro_zhou_total_raw_counts$cell_type_day_dataset=sub('arly.*','',t_petro_zhou_total_raw_counts$cell_type_day_dataset)
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,"TB", "TE")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,"TE1|TE2|TE3", "TE")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,"EightCells", "8C")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,"B1_B2", "BLASTO")
t_petro_zhou_total_raw_counts$cell_type_day_dataset=str_replace(t_petro_zhou_total_raw_counts$cell_type_day_dataset,"Morula", "MORULA")

t_petro_zhou_total_raw_counts$cell_type=sub('_[^_]*$', '',t_petro_zhou_total_raw_counts$cell_type_day_dataset)

t_petro_zhou_total_raw_counts %>%
  filter(sex=="female") %>%
  tidyr::pivot_longer(!c(cell_type,
                         dataset,
                         day,
                         sex,
                         cell_type_day_dataset), names_to = "gene_id", values_to = "count") -> petro_zhou_raw_female_tibble

petro_zhou_raw_female_tibble$gene=sub("_ENSG.*","",petro_zhou_raw_female_tibble$gene_id)

petro_zhou_raw_female_tibble %>%
  group_by(cell_type_day_dataset) %>%
  mutate(sum_cell_type=sum(count)) %>%
  ungroup() %>%
  group_by(cell_type_day_dataset,gene)%>%
  summarise(norm_gene_count=unique(sum(count)/sum_cell_type*1e6),gene=unique(gene)) -> petro_zhou_raw_female_tibble


petro_zhou_raw_female_tibble_wider=petro_zhou_raw_female_tibble
petro_zhou_raw_female_tibble_wider$cell_type=NULL

petro_zhou_raw_female_tibble_wider %>%
  tidyr::pivot_wider(names_from = cell_type_day_dataset, values_from = norm_gene_count) %>% 
  as.data.frame() -> petro_zhou_raw_female_tibble_wider

rownames(petro_zhou_raw_female_tibble_wider)=petro_zhou_raw_female_tibble_wider$gene
petro_zhou_raw_female_tibble_wider$gene=NULL


petro_zhou_raw_female_tibble_wider_quantile_norm=preprocessCore::normalize.quantiles(as.matrix(log2(petro_zhou_raw_female_tibble_wider+1)))

rownames(petro_zhou_raw_female_tibble_wider_quantile_norm)=rn(petro_zhou_raw_female_tibble_wider)
colnames(petro_zhou_raw_female_tibble_wider_quantile_norm)=cn(petro_zhou_raw_female_tibble_wider)

petro_zhou_raw_female_tibble_wider_quantile_norm=as.data.frame(petro_zhou_raw_female_tibble_wider_quantile_norm)

petro_zhou_raw_female_tibble_wider_quantile_norm$gene=rn(petro_zhou_raw_female_tibble_wider_quantile_norm)

petro_zhou_raw_female_tibble_wider_quantile_norm %>%
  tidyr::pivot_longer(!c(gene), names_to = "cell_type_day_dataset", values_to = "count") -> petro_zhou_raw_female_tibble

petro_zhou_raw_female_tibble$cell_type=sub('_[^_]*$', '',petro_zhou_raw_female_tibble$cell_type_day_dataset)

petro_zhou_raw_female_tibble %>%
  filter(gene%in%constant_gtf$gene) %>%
  mutate(category="constant") -> petro_zhou_raw_female_constant

petro_zhou_raw_female_tibble %>%
  filter(gene%in%up_gtf$gene) %>%
  mutate(category="up") -> petro_zhou_raw_female_up

petro_zhou_raw_female_tibble_plot=rbind(petro_zhou_raw_female_constant,
                                   petro_zhou_raw_female_up)



petro_zhou_raw_female_tibble_plot %>%
  mutate(log_expr=count) -> petro_zhou_raw_female_tibble_plot

petro_zhou_raw_female_tibble_plot$cell_type_day_dataset=factor(petro_zhou_raw_female_tibble_plot$cell_type_day_dataset,
                                                          levels=c("8C_3","MORULA_4","MORULA_5",
                                                                   "BLASTO_4","BLASTO_5","BLASTO_6",
                                                                   "EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                                                   "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                                                   "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12",
                                                                   "ST_8","ST_10","ST_12","EVT_6","EVT_10","EVT_12",
                                                                   "ysTE_10","ysTE_12","TE.apoptosis_8"
                                                          ))

petro_zhou_raw_female_tibble_plot$cell_type=sub("_.*","",petro_zhou_raw_female_tibble_plot$cell_type_day_dataset)


petro_zhou_raw_female_tibble_plot%>%
  filter(cell_type_day_dataset%in%c("8C_3","MORULA_4","BLASTO_4","BLASTO_5","EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                    "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                    "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12"))%>%
  group_by(category)%>%
  mutate(median_cat=median(2^(count)-1)) %>%
  group_by(cell_type_day_dataset,category)%>%
  mutate(sd=sd(2^(count)-1), median=median(2^(count)-1)) -> petro_zhou_raw_female_tibble_plot_df

petro_zhou_raw_female_tibble_plot_df%>%
  ggplot(aes(x=cell_type_day_dataset, y=2^(count)-1,fill=category)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.shape = NA)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df$median_cat)[1], color='#F8766D',linetype='dashed',size=0.75)+
  geom_hline(yintercept = unique(petro_zhou_raw_female_tibble_plot_df$median_cat)[2], color='#00BFC4',linetype='dashed',size=0.75)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("quantile cpm norm gene expression") +
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
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))

ggsave(paste(sep="_","gene_categories_in_embryo_female_quantile_cpm_norw_petro_zhou.pdf"),width=8,height=6)


# pdf(paste(sep="_","gene_categories_in_embryo_female_norm_expr.pdf"),width=15,height=10)
petro_zhou_raw_female_tibble_plot%>%
  filter(cell_type_day_dataset%in%c("8C_3","MORULA_4","BLASTO_4","BLASTO_5","EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                    "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                    "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12"))%>%
  ggplot(aes(x=cell_type_day_dataset, y=log_expr,fill=category)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.shape = NA)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ylab("Log norm gene expression") +
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
# scale_y_continuous(trans="pseudo_log")+

# dev.off()
