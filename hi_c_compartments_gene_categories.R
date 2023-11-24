setwd("C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/")

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

library(dplyr)
library(ggplot2)
library(stringr)
library(sva)
library(ggpubr)
library(rstatix)
library(intervals)

timeNow<-function(x=NULL, y=NULL, m=NULL, d=NULL , h=NULL, min=NULL,ymd_h_m=NULL){x=date();
date=unlist(str_split(pattern = "-", Sys.Date()))
time=unlist(str_split(pattern = ":",unlist(str_split(pattern=" ",Sys.time()))[2]))
return(paste(c(date,time), collapse = "_"))
}


a_compartment=read.table("A_compartment.bed", header=F, sep="\t")
colnames(a_compartment)=c("chrom","start","end")

b_compartment=read.table("B_compartment.bed", header=F, sep="\t")
colnames(b_compartment)=c("chrom","start","end")

bi_bi_up=read.table("embryo_pseudobulk/bi_bi_up.gtf", header=F, sep="\t")
colnames(bi_bi_up)=c('chrom','source','feature','start',
                     'end',
                     'score',
                     'strand',
                     'frame',
                     'attribute')

mono_bi_up=read.table("embryo_pseudobulk/mono_bi_up.gtf", header=F, sep="\t")
colnames(mono_bi_up)=c('chrom','source','feature','start',
                       'end',
                       'score',
                       'strand',
                       'frame',
                       'attribute')

bi_bi=read.table("embryo_pseudobulk/bi_bi.gtf", header=F, sep="\t")
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

unlist(sapply(a_compartment$start,function(x){ifelse(data.table::between(as.numeric(x),as.numeric(up_gtf$start),as.numeric(up_gtf$end)),"A",NA)}))

int_a <- Intervals(a_compartment[,c("start","end")])
int_b <- Intervals(b_compartment[,c("start","end")])

int_up <- Intervals(up_gtf[,c("start","end")])
rownames(int_up)=up_gtf$gene

int_constant <- Intervals(constant_gtf[,c("start","end")])
rownames(int_constant)=constant_gtf$gene




int_up_int_a=NULL
int_up_int_b=NULL

for(interval in seq(nrow(int_up))){
  
  int_up_int_a[interval]<-length(interval_intersection(int_a,int_up[interval,]))!=0 # determines if there is an intersection between the two
  int_up_int_b[interval]<-length(interval_intersection(int_b,int_up[interval,]))!=0 # determines if there is an intersection between the two
  
}

up_compartment=data.frame(row.names = rn(int_up),
           A=int_up_int_a,
           B=int_up_int_b)

up_compartment$compartment=unlist(apply(up_compartment,1,function(x){
  
  if(x["A"]&!x["B"]){return("A")}
  if(!x["A"]&x["B"]){return("B")}
  if(x["A"]&x["B"]){return("AB")}
  if(!x["A"]&!x["B"]){return("out")}
  
}))

up_compartment$category=rep("up",nrow(up_compartment))


int_constant_int_a=NULL
int_constant_int_b=NULL

for(interval in seq(nrow(int_constant))){
  
  int_constant_int_a[interval]<-length(interval_intersection(int_a,int_constant[interval,]))!=0 # determines if there is an intersection between the two
  int_constant_int_b[interval]<-length(interval_intersection(int_b,int_constant[interval,]))!=0 # determines if there is an intersection between the two
  
}

constant_compartment=data.frame(row.names = rn(int_constant),
                          A=int_constant_int_a,
                          B=int_constant_int_b)

constant_compartment$compartment=unlist(apply(constant_compartment,1,function(x){
  
  if(x["A"]&!x["B"]){return("A")}
  if(!x["A"]&x["B"]){return("B")}
  if(x["A"]&x["B"]){return("AB")}
  if(!x["A"]&!x["B"]){return("out")}
  
}))

constant_compartment$category=rep("constant",nrow(constant_compartment))


gene_category_compartment=rbind(up_compartment,constant_compartment)

pval <- stats::chisq.test(gene_category_compartment$category,gene_category_compartment$compartment,simulate.p.value = TRUE)$p.value
# CAUTION !!! PVAL CHANGES (BUT WITHIN GIVEN RANGE), AS IT IS COMPUTED FROM SIMULATION TO "INCREASE" N OBS


gene_category_compartment%>%
  ggplot(aes(x = category, y = frequency(compartment), fill = compartment)) +
  geom_bar(stat = "identity",position = "fill")+
  ggtitle('Hi-C compartment location\nof gene categories') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, by = 0.1))
  # annotate("text", x=1.5, y=1.1, label=round(pval,5))

ggsave("hi_c_compartment_gene_categories.pdf")
