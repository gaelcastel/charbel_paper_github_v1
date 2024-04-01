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
library(circlize)
library(ComplexHeatmap)

timeNow<-function(x=NULL, y=NULL, m=NULL, d=NULL , h=NULL, min=NULL,ymd_h_m=NULL){x=date();
date=unlist(str_split(pattern = "-", Sys.Date()))
time=unlist(str_split(pattern = ":",unlist(str_split(pattern=" ",Sys.time()))[2]))
return(paste(c(date,time), collapse = "_"))
}

cl <- makeCluster(10)  # Use 10 cores
registerDoParallel(cl)

sample_annot_cutrun=read.csv("sample_annot_cutrun.csv", header=T)
rownames(sample_annot_cutrun)=sample_annot_cutrun$Sample_ID
sample_annot_cutrun$cellType=unlist(lapply(str_split(sample_annot_cutrun$Sample_Name,pattern="_"),function(x){x[grep('PXGL|Primed',x)]}))
sample_annot_cutrun$condition=unlist(lapply(lapply(str_split(sample_annot_cutrun$Sample_Name,pattern="_"),function(x){x[grep('WT|KO|UT|DOX',x)]}),
                                            function(x) if(identical(x, character(0))) NA_character_ else x))

sample_annot_cutrun$condition[is.na(sample_annot_cutrun$condition)]<-"WT"
sample_annot_cutrun$mark=unlist(lapply(str_split(sample_annot_cutrun$Sample_Name,pattern="_"),function(x){x[grep('H3K|H2A|CTCF|IgG',x)]}))
sample_annot_cutrun$replicate=sub('.*CRISPRi','',unlist(lapply(str_split(sample_annot_cutrun$Sample_Name,pattern="_"),function(x){x[grep('Rep|A11|B3|B4',x)]})))
sample_annot_cutrun$replicate[sample_annot_cutrun$replicate=='A11']='Rep1'
sample_annot_cutrun$replicate[sample_annot_cutrun$replicate=='B3']='Rep2'
sample_annot_cutrun$replicate[sample_annot_cutrun$replicate=='B4']='Rep3'
sample_annot_cutrun$target=sample_annot_cutrun$mark
sample_annot_cutrun$target[sample_annot_cutrun$target!='IgG']='target'


bed_files=list.files(path = 'bed/',pattern='regions')

bed_list=NULL

for(bed in bed_files){

  sample=sub('\\.regions.*','',bed)

  bed_list[[sample]]=as.data.frame(fread(file = paste0('bed/',bed), header = F))
  colnames(bed_list[[sample]])=c('chr','start','end',sample)
  bed_list[[sample]]=bed_list[[sample]][!grepl('_',bed_list[[sample]]$chr),]
  bed_list[[sample]][,sample]=(bed_list[[sample]][,sample])/sum((bed_list[[sample]][,sample]))*1e6

}


bed_list=bed_list[rn(sample_annot_cutrun)]


bed_list_chrX=lapply(bed_list,function(x){
  x=x[x$chr=='chrX',]
  x$feature=paste(sep="_",x$chr,x$start,x$end)
  x=x[,c(5,4)]
  return(x)
})

bed_chrX_df=do.call('cbind',bed_list_chrX)
rownames(bed_chrX_df)=bed_chrX_df[,1]
bed_chrX_df=bed_chrX_df[,!grepl("feature",cn(bed_chrX_df))]
colnames(bed_chrX_df)=sub('\\..*','',cn(bed_chrX_df))

sample_annot_cutrun_subset=sample_annot_cutrun[!grepl("H3K4|ac|CTCF|IgG",sample_annot_cutrun$mark),]

bed_chrX_df_expressed=bed_chrX_df[rowSums(bed_chrX_df)>0,rn(sample_annot_cutrun_subset)]


# select most variable bins in naive

sample_annot_cutrun_naive=sample_annot_cutrun_subset[sample_annot_cutrun_subset$cellType=='PXGL',]

bed_chrX_naive=bed_chrX_df_expressed[,rn(sample_annot_cutrun_naive)]
bed_chrX_naive$rownames=rn(bed_chrX_naive)

bed_chrX_naive %>%
  tidyr::pivot_longer(!rownames, names_to = "sample", values_to = "coverage") -> bed_chrX_naive

bed_chrX_naive$mark=bed_chrX_naive$sample
bed_chrX_naive$mark=unlist(sapply(bed_chrX_naive$mark,function(x){sample_annot_cutrun_naive[x,'mark']}))


bed_chrX_naive %>%
  group_by(mark)%>%
  group_split()->bed_chrX_naive_mark_list

names(bed_chrX_naive_mark_list)=unlist(lapply(bed_chrX_naive_mark_list,function(x){unique(x$mark)}))

bed_chrX_naive_mark_list=lapply(bed_chrX_naive_mark_list,function(x){x[,c('rownames','sample','coverage')]})

bed_chrX_naive_mark_list=lapply(bed_chrX_naive_mark_list,function(x){
  
  x %>%
    tidyr::pivot_wider(names_from = sample, values_from = coverage) %>% 
    as.data.frame() -> x
  
})

var_bins_mark=lapply(bed_chrX_naive_mark_list,function(x){
  
  rownames(x)=x$rownames
  x$rownames=NULL
  
  x=x[rowSums(x)>20,]
  
  variableBins<-getMostVariableGenes4(x)
  thres<- 0.25
  qplotDensity(variableBins$residuals)+geom_vline(xintercept = thres)
  variableBinsNames<-head(rn(variableBins)[order(variableBins$residuals,decreasing = T)],500)
  len(variableBinsNames)
  
  return(variableBinsNames)
  
})


var_bins=unlist(var_bins_mark)


bed_chrX_df_expressed_var_bins=bed_chrX_df_expressed[var_bins,]

# all(cn(bed_chrX_df_expressed_var_bins)==rn(sample_annot_cutrun_subset))

colnames(bed_chrX_df_expressed_var_bins)=sample_annot_cutrun_subset$Sample_Name
colnames(bed_chrX_df_expressed_var_bins)=sub("\\..*","",cn(bed_chrX_df_expressed_var_bins))

# CAUTION !!! CHANGE WT -> KO that was inverted

new_colnames=colnames(bed_chrX_df_expressed_var_bins)
new_colnames=str_replace(new_colnames, "WT", "ko")
new_colnames=str_replace(new_colnames, "KO", "wt")
new_colnames=str_replace(new_colnames, "wt", "WT")
new_colnames=str_replace(new_colnames, "ko", "KO")

colnames(bed_chrX_df_expressed_var_bins)=new_colnames

colorScale<-colorRamp2(quantile(cor(log2(bed_chrX_df_expressed_var_bins+1)),probs=c(0.05,0.5,0.95)),c("#00008B","white","red"))

htPearson=Heatmap(cor(log2(bed_chrX_df_expressed_var_bins+1)),
                  col = colorScale,clustering_method_columns = "ward.D2",clustering_method_rows = "ward.D2",
                  show_row_names = T,show_column_names = T,row_names_gp = autoGparFontSizeMatrix(nrow(cor(log2(bed_chrX_df_expressed_var_bins+1)))),
                  name = "Pearson\ncorrelation", column_title = "Pearson correlation heatmap\nCUTRUN repressive marks\nchrX")

#pdf("pearson_correlation_heatmap_cutrun_repressive_marks_chrX_top_var_bins_on_naive_top_500_var_and_expr_thres_20_rowsums_by_mark.pdf",width = 20,height = 20)
htPearson
#dev.off() 

# CAUTION!!! IF WE TAKE (ALMOST) ALL POSITIONS, THEN CLUSTER ACCORDING TO CLONE BECAUSE NOISE DOMINATES, RATHER THAN CONDITION


bed_chrX_df_repressive_marks=bed_chrX_df
colnames(bed_chrX_df_repressive_marks)=sub("\\..*","",cn(bed_chrX_df_repressive_marks))

sample_annot_wt=sample_annot_cutrun_subset[(sample_annot_cutrun_subset$cellType=='Primed'&sample_annot_cutrun_subset$condition=='WT')|(sample_annot_cutrun_subset$cellType=='PXGL'&sample_annot_cutrun_subset$condition=='KO'),]

bed_chrX_df_repressive_marks=bed_chrX_df_repressive_marks[,rn(sample_annot_wt)]

colnames(bed_chrX_df_repressive_marks)=sample_annot_wt$Sample_Name

# CAUTION !!! CHANGE WT -> KO that was inverted

new_colnames=colnames(bed_chrX_df_repressive_marks)
new_colnames=str_replace(new_colnames, "WT", "ko")
new_colnames=str_replace(new_colnames, "KO", "wt")
new_colnames=str_replace(new_colnames, "wt", "WT")
new_colnames=str_replace(new_colnames, "ko", "KO")

colnames(bed_chrX_df_repressive_marks)=new_colnames

colorScale<-colorRamp2(quantile(cor(log2(bed_chrX_df_repressive_marks+1)),probs=c(0.05,0.5,0.95)),c("#00008B","white","red"))

htPearson=Heatmap(cor(log2(bed_chrX_df_repressive_marks+1)),
                  col = colorScale,clustering_method_columns = "ward.D2",clustering_method_rows = "ward.D2",
                  show_row_names = T,show_column_names = T,row_names_gp = autoGparFontSizeMatrix(nrow(cor(log2(bed_chrX_df_repressive_marks+1)))),
                  name = "Pearson\ncorrelation", column_title = "Pearson correlation heatmap\nCUTRUN repressive marks\nchrX\nall bins WT")

#pdf("pearson_correlation_heatmap_cutrun_repressive_marks_chrX_all_bins_wt.pdf",width = 20,height = 20)
htPearson #supp fig 4a
#dev.off() 


