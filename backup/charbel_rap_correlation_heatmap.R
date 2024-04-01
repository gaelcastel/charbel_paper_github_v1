setwd("C:/Users/gael/charbel_paper/rap/")

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

bed_list_chrx=NULL

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

pdf("pearson_correlation_heatmap_rap_xist_chrX_all_bins.pdf",width = 20,height = 20)
htPearson
dev.off() 





bed_chrX_naive_mark_list=lapply(bed_chrX_naive_mark_list,function(x){
  
  bed_chrX %>%
    tidyr::pivot_wider(names_from = sample, values_from = coverage) %>% 
    as.data.frame() -> x
  
})







# all(unique(bed_all$sample)==rn(sample_annot_rap))
bed_chrX$cellType=bed_chrX$sample
bed_chrX$mark=bed_chrX$sample
bed_chrX$replicate=bed_chrX$sample
bed_chrX$target=bed_chrX$sample


foreach(x=rn(sample_annot_rap)) %do% {

  bed_chrX$cellType[bed_chrX$cellType==x]<-sample_annot_rap[x,'cellType']
  bed_chrX$mark[bed_chrX$mark==x]<-sample_annot_rap[x,'mark']
  bed_chrX$replicate[bed_chrX$replicate==x]<-sample_annot_rap[x,'replicate']
  bed_chrX$target[bed_chrX$target==x]<-sample_annot_rap[x,'target']

}


bed_chrX %>%
  filter(mark=="RAP") 

# save(bed_chrX,file="bed_chrX.Rdata") # takes ~ 10 min

# load("bed_chrX.Rdata")

bed_all_fc=bed_all

bed_all_fc %>%
  group_by(mark) %>%
  group_split() ->bed_all_fc

names(bed_all_fc)=unlist(lapply(bed_all_fc,function(x){unique(x$mark)}))

marks=names(bed_all_fc)[names(bed_all_fc)!='Input']

bed_all_fc_thres=NULL

for(i in c(10)){
  bed_all_fc_thres<-foreach(x=names(bed_all_fc)[names(bed_all_fc)!='Input']) %do% { # takes ~15 min
    
    # bed_all_fc[[x]] %>%
    #   group_by(cellType,replicate) %>%
    #   filter(count>exp(mean(log(count+0.01)))) -> bed_all_fc[[x]] #geom_mean
    
    bed_all_fc_Input=rbind(bed_all_fc[[x]],bed_all_fc[['Input']])
    
    bed_all_fc_Input %>%
      group_by(cellType,replicate,chr,start) %>%
      filter(n()>1) -> bed_all_fc_Input # select common bins to target and Inputs
  
    bed_all_fc_Input %>%
      group_by(cellType,replicate) %>%
      summarise(geom_mean=exp(mean(log(count[mark != "Input"]+0.01))))
    
    bed_all_fc_Input %>%
      group_by(cellType,replicate) %>%
      summarise(chr = chr[mark != "Input"],start = start[mark != "Input"],end = end[mark != "Input"], fc = (count[mark != "Input"]+0.01) / (count[mark == "Input"]+0.01)) %>%
      filter(fc>=i) -> x_fc # avoid /0 ratio ; thres 5 OK
    
    x_fc$mark=rep(x,nrow(x_fc))
    
    bed_all_fc_Input=NULL
    
    return(x_fc)
    
  }
  
  # save(bed_all_fc,file="bed_all_fc.RData")
  # load("bed_all_fc.RData")
  bed_mark=do.call('rbind',bed_all_fc_thres)
  
  chrom_sizes=read.table("hg38.chrom.sizes",header = F)
  colnames(chrom_sizes)=c("chr","size")
  rownames(chrom_sizes)=chrom_sizes$chr
  
  bed_mark%>%
    group_by(chr)%>%
    group_split() -> bed_mark
  
  bed_mark=lapply(bed_mark,function(x){
    x$size=rep(chrom_sizes[unique(x$chr),'size'],nrow(x))
    return(x)})
  
  bed_mark_df=do.call('rbind',bed_mark)
  
  bed_mark_df %>%
    group_by(mark) %>%
    group_split() -> bed_mark_list
  
  names(bed_mark_list)=unlist(lapply(bed_mark_list,function(x){unique(x$mark)}))
  
  # compute % coverage per chr and plot
  for(mark in names(bed_mark_list)){
    
    x=bed_mark_list[[mark]]
    
  
    naive=x[x$cellType=="PXGL",]
    primed=x[x$cellType=="Primed",]
    
    naive%>%
      group_by(chr,start) %>%
      filter(n()>1) %>% 
      filter(row_number()==1) %>% 
      summarise(chr=unique(chr),start=unique(start),end=unique(end)) %>%
      arrange(factor(chr, levels=c(paste0(rep("chr",22),seq(22)),'chrX','chrY','chrM'))) -> naive # select common peaks across replicates
    
    primed%>%
      group_by(chr,start) %>%
      filter(n()>1) %>% 
      filter(row_number()==1) %>% 
      summarise(chr=unique(chr),start=unique(start),end=unique(end)) %>%
      arrange(factor(chr, levels=c(paste0(rep("chr",22),seq(22)),'chrX','chrY','chrM'))) -> primed 
    
    # write.tsv(naive, file=paste(sep="_",mark,'fc',i,"PXGL_RAP.bed"),row.names = F,col.names = F,  eol = "\n")
    # write.tsv(primed, file=paste(sep="_",mark,'fc',i,"Primed_RAP.bed"),row.names = F,col.names = F,  eol = "\n")
    
    x %>%
      group_by(cellType,replicate,chr) %>%
      summarise(coverage_percent=unique(n()*10000/size*100)) %>%
      filter(!(chr%in%c('chrM','chrY'))) -> x # NB: chrY neglictable, likely alignment errors
    
    x %>%
      group_by(cellType,chr) %>%
      mutate(sd=sd(coverage_percent), mean=mean(coverage_percent)) -> x
    
    x %>%
      group_by(cellType) %>%
      mutate(median=median(coverage_percent)) -> x
    
    x$chr=factor(x$chr,levels=c(paste0(rep("chr",22),seq(22)),'chrX'))
    
    print(ggplot(x,aes(x=chr,y=coverage_percent,fill=cellType))+
            geom_bar(stat = 'summary',position = "dodge")+
            geom_point(position = position_dodge(width = 0.9)) +
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge")+
            ggtitle(mark) +
            ylab('% occupancy')+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            geom_hline(yintercept = unique(x$median)[1], color='#00BFC4',linetype='dashed',size=0.75)+
            geom_hline(yintercept = unique(x$median)[2], color='#F8766D',linetype='dashed',size=0.75)+
            theme(axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank())
          
    )
    ggsave(paste0(mark,'_fc_',i,'_coverage_occupancy_rap_flip_back.pdf'))
    
    dev.off()
    
  }

}
