setwd("C:/Users/gael/charbel_paper/rna/")

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

.libPaths("C:/Program Files/R/R-4.2.3/library")

library(doParallel)
# library(BiocParallel)
library(dplyr)
library(ggplot2)
library(stringr)
library(sva)
library(ggpubr)
library(rstatix)
library(data.table)
library(tidyr)

timeNow<-function(x=NULL, y=NULL, m=NULL, d=NULL , h=NULL, min=NULL,ymd_h_m=NULL){x=date();
date=unlist(str_split(pattern = "-", Sys.Date()))
time=unlist(str_split(pattern = ":",unlist(str_split(pattern=" ",Sys.time()))[2]))
return(paste(c(date,time), collapse = "_"))
}

sample_annot=read.csv("sample_annot_gro_seq.csv", header=T)
rownames(sample_annot)=sample_annot$Sample_ID
sample_annot$condition=sample_annot$Sample_Name
sample_annot$condition[grep("_C",sample_annot$condition)]<-"control"
sample_annot$condition[grep("_T",sample_annot$condition)]<-"treatment"
sample_annot$replicate=sample_annot$Sample_Name
sample_annot$replicate[grep("A11",sample_annot$replicate)]<-"rep1"
sample_annot$replicate[grep("B3",sample_annot$replicate)]<-"rep2"
sample_annot$replicate[grep("B4",sample_annot$replicate)]<-"rep3"

matrices=list.files(path = 'matrices/',pattern='.matrix.gz')

mtx=sapply(matrices,function(x){
  fread(file = paste0('matrices/',x), header = F)
  })

rn_x=lapply(mtx,function(x){
  
  return(paste(sep="_",x$V1,x$V2,x$V3,x$V4))
})

df_list=lapply(mtx,function(x){x[,-c(1:6)]})

df_genebody=df_list[[1]]
df_tss=df_list[[2]]

rownames(df_genebody)=rn_x[[1]]
rownames(df_tss)=rn_x[[2]]
# all(rn(df_tss)==rn(df_genebody))

# split by samples

colnames(df_genebody)=c(paste(sep="_","D1271T117",as.character(seq(1,400))),paste(sep="_","D1271T118",as.character(seq(1,400))),paste(sep="_","D1271T119",as.character(seq(1,400))),paste(sep="_","D1271T120",as.character(seq(1,400))),paste(sep="_","D1271T121",as.character(seq(1,400))),paste(sep="_","D1271T122",as.character(seq(1,400))))
colnames(df_tss)=c(paste(sep="_","D1271T117",as.character(seq(1,20))),paste(sep="_","D1271T118",as.character(seq(1,20))),paste(sep="_","D1271T119",as.character(seq(1,20))),paste(sep="_","D1271T120",as.character(seq(1,20))),paste(sep="_","D1271T121",as.character(seq(1,20))),paste(sep="_","D1271T122",as.character(seq(1,20))))

ratio_tss_genebody_list=NULL

for(sample in unique(sub("_.*","",cn(df_genebody)))){
  
  ratio_tss_genebody_list[[sample]]=(rowSums(select(df_tss, cn(df_tss)[grep(sample,cn(df_tss))]))+0.001)/(rowSums(select(df_genebody, cn(df_genebody)[grep(sample,cn(df_genebody))]))+0.001)

}

# split by gene categories

# "group_boundaries":[0,117,505] in zcat gro_seq_gene_body.matrix.gz |head
# "constant.gtf","up.gtf" ; constant genes = from 1:117 ; up genes

ratio_constant=lapply(ratio_tss_genebody_list,function(x){x[1:117]})
ratio_up=lapply(ratio_tss_genebody_list,function(x){x[118:505]})

ratio_constant_df=as.data.frame(do.call('cbind',ratio_constant))
ratio_up_df=as.data.frame(do.call('cbind',ratio_up))

ratio_constant_df$gene=rn(df_genebody)[1:117]
ratio_up_df$gene=rn(df_genebody)[118:505]

ratio_constant_df %>%
  tidyr::pivot_longer(!gene, names_to = "sample", values_to = "log2ratio") -> ratio_constant_df

ratio_up_df %>%
  tidyr::pivot_longer(!gene, names_to = "sample", values_to = "log2ratio") -> ratio_up_df

ratio_constant_df$condition=ratio_constant_df$sample
ratio_up_df$condition=ratio_up_df$sample

ratio_constant_df$condition=unlist(sapply(ratio_constant_df$condition,function(x){sample_annot[x,'condition']}))
ratio_up_df$condition=unlist(sapply(ratio_up_df$condition,function(x){sample_annot[x,'condition']}))

ratio_constant_df$replicate=ratio_constant_df$sample
ratio_up_df$replicate=ratio_up_df$sample

ratio_constant_df$replicate=unlist(sapply(ratio_constant_df$replicate,function(x){sample_annot[x,'replicate']}))
ratio_up_df$replicate=unlist(sapply(ratio_up_df$replicate,function(x){sample_annot[x,'replicate']}))

ratio_constant_df$category=rep("constant",nrow(ratio_constant_df))
ratio_up_df$category=rep("up",nrow(ratio_up_df))

ratio_constant_df$sample_name=ratio_constant_df$sample
ratio_up_df$sample_name=ratio_up_df$sample

ratio_constant_df$sample_name=unlist(sapply(ratio_constant_df$sample_name,function(x){sample_annot[x,'Sample_Name']}))
ratio_up_df$sample_name=unlist(sapply(ratio_up_df$sample_name,function(x){sample_annot[x,'Sample_Name']}))


ratio_df=rbind(ratio_constant_df,ratio_up_df)

pdf(file = "gro_seq_tss_gene_body_ratio_by_gene_categories_upon_xist_depletion.pdf",width=5,height = 7)

ratio_df %>%
  group_by(category,condition,replicate)%>%
  summarise(median_ratio=median(log2ratio),sample=unique(sample),sample_name=unique(sample_name)) %>%
  ungroup()%>%
  group_by(category,condition)%>%
  mutate(mean=mean(median_ratio),sd=sd(median_ratio))%>%
  ggplot(aes(x=condition,y=median_ratio,fill=condition))+
    geom_bar(stat = 'summary')+
    geom_point()+
    facet_wrap(~category, nrow=1, ncol=2)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge") +
    stat_compare_means(ref.group = "control", vjust=0.2, label="p.signif",size=3, method="t.test", paired=F)  +
    stat_compare_means(ref.group = "control", vjust=1, label="p.format",size=3, method="t.test", paired=F)  + 
    ggtitle('GRO-Seq TSS/gene_body ratio\nby gene categories upon\nXIST depletion') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
