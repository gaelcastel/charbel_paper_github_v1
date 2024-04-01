setwd("../atac/")

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


constant=read.table("compute_matrices_for_boxplot/constant.gtf", header=F, sep="\t")
colnames(constant)=c('chrom','source','feature','start',
                     'end',
                     'score',
                     'strand',
                     'frame',
                     'attribute')
constant$category=rep("constant",nrow(constant))

up=read.table("compute_matrices_for_boxplot/up.gtf", header=F, sep="\t")
colnames(up)=c('chrom','source','feature','start',
               'end',
               'score',
               'strand',
               'frame',
               'attribute')
up$category=rep("up",nrow(up))

transcripts_gtf=rbind(constant,up)

transcripts_gtf %>%
  filter(feature=="transcript") -> transcripts_gtf

transcripts_gtf$transcript=sub(".*transcript_id ","",unlist(lapply(str_split(pattern=";",transcripts_gtf$attribute),function(x){x[grep("transcript_id",x)]})))

matrices=list.files(path = 'matrices_boxplot/',pattern='.matrix.gz')
names(matrices)=sub("_up.*","",matrices)
matrices=as.list(matrices)

df_list=lapply(matrices,function(x){
  df=fread(file = paste0('matrices_boxplot/',x),data.table = F)
  rownames(df)=df$V4
  df=df[,7:ncol(df)]
  naive_ctl_1=rowSums(df[,1:200])
  naive_ctl_2=rowSums(df[,201:400])
  naive_ctl_3=rowSums(df[,401:600])
  
  naive_ttt_1=rowSums(df[,601:800])
  naive_ttt_2=rowSums(df[,801:1000])
  naive_ttt_3=rowSums(df[,1001:1200])
  
  ctl_df=data.frame(naive_ctl_1,naive_ctl_2,naive_ctl_3)
  ttt_df=data.frame(naive_ttt_1,naive_ttt_2,naive_ttt_3)
  

  
  return(data.frame(ctl=rowMeans(ctl_df),ttt=rowMeans(ttt_df)))
  
})

# all(rn(df)==transcripts_gtf$transcript)

df=do.call('cbind',df_list)
df$gene=rn(df)
df$category=transcripts_gtf$category

df%>%
  tidyr::pivot_longer(!c(gene,category), names_to = "sample", values_to = "norm_count") -> df


df$condition=df$sample
df$condition=sub(".*\\.","",df$condition)

df %>%
  filter(condition=="ctl") %>%
  ggplot() +
  geom_boxplot(aes(x=category, y=norm_count, fill=category),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  # stat_boxplot(aes(x=category, y=tpm),geom='errorbar', linetype=1, width=0.8, size=0.25)+
  stat_compare_means(aes(x=category, y=norm_count, group=category),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("ATAC-Seq XIST sensitive/unsensitive gene categories\nin WT naive - TSS 5kb")

# ggsave("atac_seq_gene_categories_in_wt_naive_tss_5kb.pdf",width=5,height=7)




df %>%
  ggplot() +
  geom_boxplot(aes(x=category, y=norm_count, fill=condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  # stat_boxplot(aes(x=category, y=tpm),geom='errorbar', linetype=1, width=0.8, size=0.25)+
  stat_compare_means(aes(x=category, y=norm_count, group=condition),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("ATAC-Seq XIST sensitive/unsensitive gene categories\nNaive ctl vs ttt - TSS 5kb")

# ggsave("atac_seq_gene_categories_naive_ctl_ttt_tss_5kb.pdf",width=5,height=7)
