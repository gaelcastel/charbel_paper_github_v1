setwd("../cutrun/")

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
  primed=rowSums(df[,1:200])
  naive_ctl=rowSums(df[,201:400])
  naive_ttt=rowSums(df[,401:600])
  
  return(data.frame(primed=primed,naive_ctl=naive_ctl,naive_ttt=naive_ttt))
  
})

# all(rn(df)==transcripts_gtf$transcript)

df=do.call('cbind',df_list)
df$gene=rn(df)
df$category=transcripts_gtf$category

df%>%
  tidyr::pivot_longer(!c(gene,category), names_to = "sample", values_to = "norm_count") -> df


df$mark=df$sample
df$mark=sub("\\..*","",df$mark)

df$condition=df$sample
df$condition=sub(".*\\.","",df$condition)

df %>%
  filter(condition%in%c("primed","naive_ctl"),mark!="IgG") %>%
  ggplot() +
  geom_boxplot(aes(x=category, y=norm_count, fill=condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~mark, nrow=2, ncol=3,scales = 'free')+
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
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))+
  # stat_compare_means(aes(x=category, y=norm_count, group=condition),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("Chromatin marks CUTRUN gene categories 5kb around TSS\nin WT naive vs Primed") -> p


res.stat <- df %>%
  filter(condition%in%c("primed","naive_ctl"),mark!="IgG") %>%
  group_by(mark,category)%>%
  t_test(data = ., norm_count ~ condition,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "category", fun = "max",formula = norm_count ~ condition,dodge = 0.9) %>%
  ungroup()

#pdf("cutrun_boxplot_wt_naive_primed_tss_5kb.pdf",width = 7,height = 5)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif", y.position=log(res.stat$y.position,2.7)
        )

)

#dev.off()



df %>%
  filter(condition%in%c("naive_ctl"),mark!="IgG") %>%
  ggplot() +
  geom_boxplot(aes(x=category, y=norm_count, fill=category),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~mark, nrow=2, ncol=3,scales = 'free')+
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
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))+
  # stat_compare_means(aes(x=category, y=norm_count, group=condition),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("Chromatin marks CUTRUN gene categories 5kb around TSS\nin WT naive") -> p


res.stat <- df %>%
  filter(condition%in%c("naive_ctl"),mark!="IgG") %>%
  group_by(mark)%>%
  t_test(data = ., norm_count ~ category,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "category", fun = "max",formula = norm_count ~ category,dodge = 0.9) %>%
  ungroup()

#pdf("cutrun_boxplot_wt_naive_tss_5kb.pdf",width = 7,height = 5)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif", y.position=log(res.stat$y.position,2.7)
        )
      
)

#dev.off()




df %>%
  filter(condition%in%c("naive_ctl","naive_ttt"),mark!="IgG") %>%
  ggplot() +
  geom_boxplot(aes(x=category, y=norm_count, fill=condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~mark, nrow=2, ncol=3,scales = 'free')+
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
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))+
  # stat_compare_means(aes(x=category, y=norm_count, group=condition),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("Chromatin marks CUTRUN gene categories 5kb around TSS\nin WT naive ctl vs ttt") -> p


res.stat <- df %>%
  filter(condition%in%c("naive_ctl","naive_ttt"),mark!="IgG") %>%
  group_by(mark,category)%>%
  t_test(data = ., norm_count ~ condition,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "category", fun = "max",formula = norm_count ~ condition,dodge = 0.9) %>%
  ungroup()

#pdf("cutrun_boxplot_naive_ctl_ttt_tss_5kb.pdf",width = 7,height = 5)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif", y.position=log(res.stat$y.position,2.7)
        )
      
) # fig 5e

#dev.off()




df %>%
  filter(mark!="IgG") %>%
  ggplot() +
  geom_boxplot(aes(x=category, y=norm_count, fill=condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  # geom_point(aes(x=category, y=norm_count, fill=condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~mark, nrow=2, ncol=3,scales = 'free')+
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
  scale_y_continuous(trans="pseudo_log",limits = c(0, NA))+
  # stat_compare_means(aes(x=category, y=norm_count, group=condition),vjust=0.1,label="p.signif",size=3, method="t.test", paired=FALSE) +
  ggtitle("Chromatin marks CUTRUN gene categories 5kb around TSS\nall cell types") -> p


res.stat <- df %>%
  filter(mark!="IgG") %>%
  group_by(mark,category)%>%
  t_test(data = ., norm_count ~ condition,paired = F) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_significance("p") %>%
  ungroup()%>%
  add_xy_position(x = "category", fun = "max",formula = norm_count ~ condition,dodge = 0.9) %>%
  ungroup()

#pdf("cutrun_boxplot_all_cell_types_tss_5kb.pdf",width = 7,height = 5)

print(p + 
        stat_pvalue_manual(
          res.stat, label = "p.signif", y.position=log(res.stat$y.position,2.7)
        )
      
)

#dev.off()
