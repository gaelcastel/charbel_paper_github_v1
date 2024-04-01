setwd("../rap/")

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

cl <- makeCluster(10)  # Use 10 cores
registerDoParallel(cl)

sample_annot_rap=read.csv("sample_annot_rap.csv", header=T)
rownames(sample_annot_rap)=sample_annot_rap$Sample_ID
sample_annot_rap$cellType=unlist(lapply(str_split(sample_annot_rap$Sample_Name,pattern="_"),function(x){x[grep('PXGL|Primed',x)]}))

sample_annot_rap$mark=unlist(lapply(str_split(sample_annot_rap$Sample_Name,pattern="_"),function(x){x[grep('Input|RAP',x)]}))
sample_annot_rap$replicate=unlist(lapply(str_split(sample_annot_rap$Sample_Name,pattern="_"),function(x){x[grep('Rep',x)]}))

sample_annot_rap$target=sample_annot_rap$mark
sample_annot_rap$target[sample_annot_rap$target!='Input']='target'


bed_files=list.files(path = 'bed/',pattern='regions')

bed_list_chrX=NULL

for(bed in bed_files){

  sample=sub('\\.regions.*','',bed)

  bed_list_chrX[[sample]]=as.data.frame(fread(file = paste0('bed/',bed), header = F))
  colnames(bed_list_chrX[[sample]])=c('chr','start','end',sample)
  bed_list_chrX[[sample]]=bed_list_chrX[[sample]][!grepl('_',bed_list_chrX[[sample]]$chr),]
  bed_list_chrX[[sample]]=bed_list_chrX[[sample]][bed_list_chrX[[sample]]$chr=="chrX",]
  bed_list_chrX[[sample]]$feature=paste(sep="_",bed_list_chrX[[sample]]$chr,bed_list_chrX[[sample]]$start,bed_list_chrX[[sample]]$end)
  bed_list_chrX[[sample]]=bed_list_chrX[[sample]][,c(5,4)]
  rownames(bed_list_chrX[[sample]])=bed_list_chrX[[sample]]$feature
  bed_list_chrX[[sample]]$feature=NULL

}


bed_list_chrX=bed_list_chrX[rn(sample_annot_rap)]
bed_list_chrX_df=do.call('cbind',bed_list_chrX)

bed_list_chrX_df=bed_list_chrX_df[,rn(sample_annot_rap)]
colnames(bed_list_chrX_df)=sample_annot_rap$Sample_Name

colorScale<-colorRamp2(quantile(cor(log2(bed_list_chrX_df+1)),probs=c(0.05,0.5,0.95)),c("#00008B","white","red"))

htPearson=Heatmap(cor(log2(bed_list_chrX_df+1)),
                  col = colorScale,clustering_method_columns = "ward.D2",clustering_method_rows = "ward.D2",
                  show_row_names = T,show_column_names = T,row_names_gp = autoGparFontSizeMatrix(nrow(cor(log2(bed_list_chrX_df+1)))),
                  name = "Pearson\ncorrelation", column_title = "Pearson correlation heatmap\nRAP-XIST\nchrX")

#pdf("pearson_correlation_heatmap_rap_xist_chrX_all_bins.pdf",width = 20,height = 20)
htPearson #fig 
# dev.off() 
