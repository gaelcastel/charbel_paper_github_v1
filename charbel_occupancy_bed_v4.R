setwd("C:/Users/gael/charbel_paper/cutrun/")

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

bed_list=lapply(bed_list,function(x){
  x %>%
    tidyr::pivot_longer(!c(chr,start,end), names_to = "sample", values_to = "count") -> x
  return(x)

})


bed_all=do.call('rbind',bed_list)

# all(unique(bed_all$sample)==rn(sample_annot_cutrun))
bed_all$cellType=bed_all$sample
bed_all$condition=bed_all$sample
bed_all$mark=bed_all$sample
bed_all$replicate=bed_all$sample
bed_all$target=bed_all$sample


foreach(x=rn(sample_annot_cutrun)) %do% {

  bed_all$cellType[bed_all$cellType==x]<-sample_annot_cutrun[x,'cellType']
  bed_all$condition[bed_all$condition==x]<-sample_annot_cutrun[x,'condition']
  bed_all$mark[bed_all$mark==x]<-sample_annot_cutrun[x,'mark']
  bed_all$replicate[bed_all$replicate==x]<-sample_annot_cutrun[x,'replicate']
  bed_all$target[bed_all$target==x]<-sample_annot_cutrun[x,'target']

}

# save(bed_all,file="bed_all.Rdata") # takes ~ 10 min

# load("bed_all.Rdata")

# for(x in c('D725C14')){ # just a fucking bug during foreach
#   bed_all$cellType[bed_all$cellType==x]<-sample_annot_cutrun[x,'cellType']
#   bed_all$condition[bed_all$condition==x]<-sample_annot_cutrun[x,'condition']
#   bed_all$mark[bed_all$mark==x]<-sample_annot_cutrun[x,'mark']
#   bed_all$replicate[bed_all$replicate==x]<-sample_annot_cutrun[x,'replicate']
#   bed_all$target[bed_all$target==x]<-sample_annot_cutrun[x,'target']
# }

for (i in c(1.2)){
  bed_all_fc=bed_all
  
  bed_all_fc %>%
    group_by(mark) %>%
    group_split() ->bed_all_fc
  
  names(bed_all_fc)=unlist(lapply(bed_all_fc,function(x){unique(x$mark)}))
  
  marks=names(bed_all_fc)[names(bed_all_fc)!='IgG']
  
  bed_all_fc<-foreach(x=names(bed_all_fc)[names(bed_all_fc)!='IgG']) %do% { # takes ~15 min
    
    bed_all_fc_igg=rbind(bed_all_fc[[x]],bed_all_fc[['IgG']])
    
    bed_all_fc_igg %>%
      group_by(cellType,condition,replicate,chr,start) %>%
      filter(n()>1) -> bed_all_fc_igg # select common bins to target and IgGs
    
    bed_all_fc_igg %>%
      group_by(cellType,condition,replicate) %>%
      summarise(chr = chr[mark != "IgG"],start = start[mark != "IgG"],end = end[mark != "IgG"], fc = (count[mark != "IgG"]+0.01) / (count[mark == "IgG"]+0.01)) %>%
      filter(fc>=i) -> x_fc # avoid /0 ratio ; thres 5 OK
    
    x_fc$mark=rep(x,nrow(x_fc))
    
    bed_all_fc_igg=NULL
    
    return(x_fc)
    
  }
  
  # save(bed_all_fc,file="bed_all_fc.RData")
  # load("bed_all_fc.RData")
  bed_mark=do.call('rbind',bed_all_fc)
  
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
    
    x$condition[x$cellType=="PXGL" & x$condition=="WT"]<-'ko' # !!!! CHANGE PXGL KO TO WT !!!!
    x$condition[x$cellType=="PXGL" & x$condition=="KO"]<-'wt'
    x$condition<-toupper(x$condition)
    
    naive=x[x$cellType=="PXGL" & x$condition=="WT",]
    primed=x[x$cellType=="Primed" & x$condition=="WT",]
    
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
    
    write.tsv(naive, file=paste(sep="_",mark,'fc',i,"PXGL_WT.bed"),row.names = F,col.names = F,  eol = "\n")
    write.tsv(primed, file=paste(sep="_",mark,'fc',i,"Primed_WT.bed"),row.names = F,col.names = F,  eol = "\n")
    
    x %>%
      group_by(cellType,condition,replicate,chr) %>%
      summarise(coverage_percent=unique(n()*10000/size*100)) %>%
      filter(condition=='WT',!(chr%in%c('chrM','chrY'))) -> x # NB: chrY neglictable, likely alignment errors
    
    x %>%
      group_by(cellType,condition,chr) %>%
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
            ylim(0,75)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            geom_hline(yintercept = unique(x$median)[1], color='#00BFC4',linetype='dashed',size=0.75)+
            geom_hline(yintercept = unique(x$median)[2], color='#F8766D',linetype='dashed',size=0.75)+
            theme(axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank()) +
            coord_flip()        
          
    )
    ggsave(paste0(mark,'_fc_',i,'_coverage_occupancy_wt.pdf'))
    
    dev.off()
    
  }
}

