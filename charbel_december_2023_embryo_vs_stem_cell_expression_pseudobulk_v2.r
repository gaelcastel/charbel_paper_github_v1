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
library(SingleCellExperiment)
library(Matrix.utils)

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

# Embryo

petro_zhou_sample_annot=fastRead("petro_zhou_sample_annot_cell_type_sex.tsv", header=T, sep="\t",as.matrix = F)

# PETRO ZHOU RAW COUNTS
load("2022_06_20_14_03_petro_total_raw_counts_not_corr.Rdata") #petro_total_raw_counts
load("2022_06_20_14_22_zhou_total_raw_counts_not_corr.Rdata") #zhou_total_raw_counts


# all(rn(assay(petro_total_raw_counts))==rn(assay(zhou_total_raw_counts)))
petro_total_raw=assay(petro_total_raw_counts)
zhou_total_raw=assay(zhou_total_raw_counts)

petro_zhou_raw=cbind(petro_total_raw,zhou_total_raw)

# petro_zhou_norm=preprocessCore::normalize.quantiles(log2(CPM(petro_zhou_raw_counts)+1))

colnames(petro_zhou_raw)=cn(petro_zhou_raw)
rownames(petro_zhou_raw)=rn(petro_zhou_raw)

petro_zhou_raw=as.data.frame(petro_zhou_raw)


petro_zhou_sample_annot=petro_zhou_sample_annot[rn(petro_zhou_sample_annot)%in%cn(petro_zhou_raw),]
petro_zhou_raw=petro_zhou_raw[,cn(petro_zhou_raw)%in%rn(petro_zhou_sample_annot)]

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

selected_genes=c("XIST", "XACT", "SPEN", "NANOG", "GATA3")

petro_zhou_raw_female_pseudobulk=petro_zhou_raw_female_pseudobulk[petro_zhou_raw_female_pseudobulk$gene%in%selected_genes,]

petro_zhou_raw_female_pseudobulk$gene=NULL

t_petro_zhou_cpm_female=as.data.frame(t(petro_zhou_raw_female_pseudobulk))

t_petro_zhou_cpm_female$cell_type_day=rn(t_petro_zhou_cpm_female)

t_petro_zhou_cpm_female %>%
  tidyr::pivot_longer(!c(cell_type_day), names_to = "gene_id", values_to = "count") -> petro_zhou_cpm_female_tibble

petro_zhou_cpm_female_tibble$gene=sub("_ENSG.*","",petro_zhou_cpm_female_tibble$gene_id)
petro_zhou_cpm_female_tibble$day=sub(".*_","",petro_zhou_cpm_female_tibble$cell_type_day)
petro_zhou_cpm_female_tibble$day=factor(petro_zhou_cpm_female_tibble$day,levels=c("3","4","5","6","7","8","10","12"))
petro_zhou_cpm_female_tibble$cell_type_day=factor(petro_zhou_cpm_female_tibble$cell_type_day,
                                                          levels=c("8C_3","MORULA_4","MORULA_5",
                                                                   "BLASTO_4","BLASTO_5",
                                                                   "EPI_5","EPI_6","EPI_7","EPI_8","EPI_10","EPI_12",
                                                                   "PrE_6","PrE_7","PrE_8","PrE_10","PrE_12",
                                                                   "TE_5", "TE_6","TE_7","TE_8","TE_10","TE_12",
                                                                   "ST_8","ST_10","ST_12","EVT_6","EVT_10","EVT_12"
                                                          ))

petro_zhou_cpm_female_tibble$cell_type=sub("_.*","",petro_zhou_cpm_female_tibble$cell_type_day)


# Stem cells


setwd('D:/rougeulle_lab_ordi_sept_2023/charbel_paper/rna/charbel_paper_ifb/')


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

df_plot=mmat_mtx_cpm[c("XIST", "XACT", "SPEN", "NANOG", "GATA3"),]
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


# Embryo vs Stem cell

df_plot$cell_type_day=df_plot$cell_type_condition
df_plot$count=df_plot$CPM

df_plot_stem_cells=df_plot[,c("cell_type_day","gene","count")]
df_plot_embryo=petro_zhou_cpm_female_tibble[,c("cell_type_day","gene","count")]

df_plot_embryo_vs_stem_cells=rbind(df_plot_embryo,df_plot_stem_cells)

df_plot_embryo_vs_stem_cells %>%
  group_by(gene,cell_type_day) %>%
  mutate(sd=sd(count), mean=mean(count)) %>%
  ggplot() + 
  geom_bar(aes(x=cell_type_day,y=count,fill=cell_type_day),stat = 'summary',fun="mean",position = 'dodge')+
  geom_point(aes(x=cell_type_day,y=count,group=cell_type_day),position = position_dodge(width = 0.9))+
  facet_wrap(~gene, nrow=2, ncol=4, scales = "free")+
  geom_errorbar(aes(x=cell_type_day,y=count,group=cell_type_day,ymin=mean-sd, ymax=mean+sd), position="dodge") +
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  ggtitle("Embryo (pseudobulk) vs Stem cell (bulk)\nGene expression")+
  #scale_x_discrete(expand = c(0, 0.5)) +
  ylab("CPM") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
  # scale_y_continuous(breaks = seq(0,1.1, 0.1)) + ylim(0,1)
  # stat_compare_means(aes(x=cell_type_day,y=count),vjust=1,label="p.signif",size=3,ref.group = "pxgl_siScr", method="t.test", paired=FALSE)
