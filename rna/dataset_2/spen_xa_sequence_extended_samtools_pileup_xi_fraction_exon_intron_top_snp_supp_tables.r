# vcf from bcftools already on SNP in H9 CHR X !!!!! BUT CAUTION !!! :not only heterozygous SNP !!! so I applied detect_bi on informative_H9_wgs.vcf
# Could be relevant to re-run bcftools intersect on this vcf -> bed only informative SNP heterozygous chrX to re-generate samples.vcf

setwd('C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/')

# ex. biallelic: DP4=0,3,2,0 : 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles
# ex. monoallelic: ; DP4=0,0,1,1: 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles. BUT also DP4=1,1,10,1 (proba/stat)

# NB: each dataset distinct read length, and some are single-end, so I guess need to determine heterozygous SNPs
# in primed of each dataset

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

# BiocManager::install("statebins")
# BiocManager::install("wesanderson")
library(statebins)
library(wesanderson)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(Matrix)
library(stringr)
library(BiocManager)
library(devtools)
library(pcaMethods)
library(data.table)
library(parallel)
library(DESeq2)
library(ggpubr)
library(tidyverse)
library(EnhancedVolcano)
library(biomaRt)
library(sva)
library(gtools)

open_file <- function(sample){
    snp <- read.table(paste0("dataset_2/vcf/", sample, "_samtools_mpileup_informative_from_wgs_h9_chrX.pileup"), header = FALSE,
                   sep = "\t")
    colnames(snp) <- c("chrom", "pos", "ref", "dp", "alts", "bq", "mq")
    return(snp)
}

detect_bi <- function(snp){
    snp=snp[snp$allele_min_num/snp$allele_maj_num >= 1/4,]
    return(snp)
}


qual_conv <- function(ascii){
  
  return(as.numeric(ascii_qual[ascii,"ASCII"]))

}


ascii_qual=read.tsv("ascii_qual.tsv",comment.char = "")
ascii_qual$symbol=rn(ascii_qual)

gtf=read.table("hg38.gtf", header=F, sep="\t")
colnames(gtf)=c('chrom','source','feature','start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute')

gtf=gtf[-grep('_',gtf$chrom),]

gtf %>%
  group_by(chrom,attribute) %>%
  summarise(chrom=unique(chrom),strand=unique(strand),start=min(start), end=max(end), feature=unique(attribute)) %>%
  as.data.frame() -> gtf_unique


gtf_unique %>%
  group_by(chrom) %>%
  arrange(start, .by_group=TRUE) %>%
  as.data.frame() -> gtf_unique

gtf_chrX=gtf_unique[gtf_unique$chrom=='chrX',]

chrom_sizes <- structure(list(chromosome = c("chrX"), size = c(156040895L)), .Names = c("chromosome", 
                                                                                        "size"), class = "data.frame", row.names = c(NA, -1L))
chrom_sizes

# hg38 centromere locations
centromeres <- structure(list(chromosome = c("chrX"), start = c(58100000L), end = c(63800000L)), 
                         .Names = c("chromosome", "start", "end"), class = "data.frame", row.names = c(NA, -1L
                         ))
centromeres

# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chrX")
chrom_key <- setNames(object = as.character(c(1)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                      levels = chrom_order)

centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                      levels = chrom_order)
# create a color key for the plot
#group.colors <- c(gain = "red", loss = "blue")

informative_pos_vcf_chrX=read.table('C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/H9_WGS_hg38_filtered_chrX.vcf')
colnames(informative_pos_vcf_chrX)=c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample")
informative_pos_vcf_chrX=informative_pos_vcf_chrX[grep("0/0|0/1|1/1", informative_pos_vcf_chrX$sample),] #remove alt_2 because very low frequency/probability
# and annoying for computation
informative_pos_vcf_chrX=informative_pos_vcf_chrX[informative_pos_vcf_chrX$chrom=="chrX" &
                                                    informative_pos_vcf_chrX$qual >=50,]
informative_pos_vcf_chrX$GT=sapply(str_split(pattern=":", informative_pos_vcf_chrX$sample),function(x) x[1])
informative_pos_vcf_chrX$AD=sapply(str_split(pattern=":", informative_pos_vcf_chrX$sample),function(x) x[2])
informative_pos_vcf_chrX$DP=as.numeric(sapply(str_split(pattern=":", informative_pos_vcf_chrX$sample),function(x) x[3]))
informative_pos_vcf_chrX$GQ=as.numeric(sapply(str_split(pattern=":", informative_pos_vcf_chrX$sample),function(x) x[4]))
informative_pos_vcf_chrX$PL=sapply(str_split(pattern=":", informative_pos_vcf_chrX$sample),function(x) x[5])
informative_pos_vcf_chrX=informative_pos_vcf_chrX[informative_pos_vcf_chrX$DP>=10,]
alleles=lapply(lapply(str_split(pattern=",", informative_pos_vcf_chrX$AD),function(x) x), as.numeric)
alleles=do.call('rbind',alleles)
informative_pos_vcf_chrX$allele_0=alleles[,1]
informative_pos_vcf_chrX$allele_1=alleles[,2]
informative_pos_vcf_chrX$allele_maj_num=apply(informative_pos_vcf_chrX,1,function(x){max(as.numeric(x['allele_0']),as.numeric(x['allele_1']))})
informative_pos_vcf_chrX$allele_min_num=apply(informative_pos_vcf_chrX,1,function(x){min(as.numeric(x['allele_0']),as.numeric(x['allele_1']))})
informative_pos_vcf_chrX=detect_bi(informative_pos_vcf_chrX)

sample_annot=read.csv("dataset_2/SampleSheet.csv", header=T, comment.char = '#')
rownames(sample_annot)=paste0("D1507",sample_annot$Sample_ID)
sample_annot$sample_id=paste0("D1507",sample_annot$Sample_ID)
sample_annot$group=sub('_[^_]*$', '', sample_annot$Sample_Name)
sample_annot$group=str_replace_all(sample_annot$group,"CTL","1")
sample_annot$group=str_replace_all(sample_annot$group,"Mix","2")
sample_annot$replicate=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$replicate[sample_annot$replicate=="primed"]<-"Rep1"
sample_annot$cell_type=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$cell_type[sample_annot$cell_type!="primed"]<-"pxgl"

sample_list=NULL
sample_list_bi=NULL
sample_list_mono=NULL
percent_biallelism=NULL
alleles=NULL
alts_list=NULL
strand_list=NULL
bq_list=NULL
mq_list=NULL

# ! different VCF format from bcftools. e.g:Number of 1) forward ref alleles; 2) reverse ref; 3) forward alt; 4) reverse alt alleles

for (sample in rn(sample_annot)){
  
    sample_list[[sample]]=open_file(sample)
    # nrow(sample_list[[sample]])==nrow(sample_list[[sample]][sample_list[[sample]]$pos %in% informative_pos_vcf_chrX$pos,]) # TRUE ; no need, because samtools pileup performed on these pos
    sample_list[[sample]]=sample_list[[sample]][sample_list[[sample]]$dp>=10,]
    
    alts_list=strsplit(sample_list[[sample]]$alts, "")
    
    alts_list=lapply(alts_list,function(x){x[!(x%in%c('*','#'))]}) #remove deletions
    
    strand_list=lapply(alts_list, function(x){
      
      x[str_detect(x,"[[:upper:]]")]<-"+"
      x[str_detect(x,"[[:lower:]]")]<-"-"
      x[x=="."]<-"+"
      x[x==","]<-"-"
    
    })

    for (i in seq(length(alts_list))){
      
      alts_list[[i]][alts_list[[i]]%in%c(".",",")]<-sample_list[[sample]][i,'ref']
      
    }  
    
    alts_list=lapply(alts_list,function(x){toupper(x)})
    
    alts_list=lapply(alts_list, function(x){
      
      if(dim(table(x))>2){x[x%in%(names(table(x)[order(table(x),decreasing = T)][1:2]))]}else{x}
      
      }) # remove putative sequencing errors (more than 2 base types detected)
        
    sample_list[[sample]]$dp=unlist(lapply(alts_list,function(x){length(x)}))
    
    sample_list[[sample]]$allele_maj_num=unlist(lapply(alts_list,function(x){table(x)[order(table(x),decreasing = T)][1]}))
    sample_list[[sample]]$allele_maj=names(unlist(lapply(alts_list,function(x){table(x)[order(table(x),decreasing = T)][1]})))
    sample_list[[sample]]$allele_min_num=unlist(lapply(alts_list,function(x){tail(table(x)[table(x)==min(table(x))],n=1)}))
    sample_list[[sample]]$allele_min=names(unlist(lapply(alts_list,function(x){tail(table(x)[table(x)==min(table(x))],n=1)})))
    
    sample_list[[sample]]$ref=toupper(sample_list[[sample]]$ref)
    
    informative_pos_vcf_chrX_subset=informative_pos_vcf_chrX[informative_pos_vcf_chrX$pos%in%sample_list[[sample]]$pos,]
    # informative_pos_vcf_chrX_subset$pos==sample_list[[sample]]$pos
    # all(informative_pos_vcf_chrX_subset$pos==sample_list[[sample]]$pos) # because ordered
    
    # CAUTION: avoid this situation, from putative sequencing error:
    # "chrom" "pos" "ref" "dp" "alts" "bq" "mq" "allele_maj_num" "allele_maj" "allele_min_num" "allele_min" "allelism"      
    # chrX 15702622 G 13 ,,,....,,,,., FFFFFkFFFFFFk 51 13 G 0 A mono
    # chrX 15702622 G 12 ,.,.,t,,,,., FkFFFFFFFFFF 52 11 G 1 T mono
    # NOT same allele_min. Take the one from the info_chrX that is different from allele_maj
    
    for (i in seq(nrow(sample_list[[sample]]))){
      
      if(sample_list[[sample]]$allele_min[i]==sample_list[[sample]]$allele_maj[i] |
         !(sample_list[[sample]]$allele_min[i]%in%c(informative_pos_vcf_chrX_subset$ref[i],
                                                    informative_pos_vcf_chrX_subset$alt[i]))){
        
        sample_list[[sample]]$allele_min[i]=c(informative_pos_vcf_chrX_subset$ref[i],informative_pos_vcf_chrX_subset$alt[i])[
        c(informative_pos_vcf_chrX_subset$ref[i],informative_pos_vcf_chrX_subset$alt[i])!=sample_list[[sample]]$allele_maj[i]][1]
        
      
        sample_list[[sample]]$allele_min_num[i]=0
      
      }else{
          
        sample_list[[sample]]$allele_min[i]=sample_list[[sample]]$allele_min[i]
        sample_list[[sample]]$allele_min_num[i]=sample_list[[sample]]$allele_min_num[i]
        
        }
      
    }

    mq_list=strsplit(sample_list[[sample]]$mq, "") # convert ascii bq & mq to numeric values, and compute median
    sample_list[[sample]]$mq=unlist(lapply(mq_list,function(x){median(qual_conv(x),na.rm=T)}))
    sample_list[[sample]]=sample_list[[sample]][sample_list[[sample]]$mq>=50,]
    
    sample_list[[sample]]$allelism=rep(NA,nrow(sample_list[[sample]]))
    sample_list[[sample]]$allelism[sample_list[[sample]]$pos%in%detect_bi(sample_list[[sample]])$pos]<-"bi"
    sample_list[[sample]]$allelism[!(sample_list[[sample]]$pos%in%detect_bi(sample_list[[sample]])$pos)]<-"mono"
    
    sample_list[[sample]]=sample_list[[sample]][nchar(sample_list[[sample]]$allele_min)<2 & nchar(sample_list[[sample]]$allele_maj)<2,] # remove indels positions that were not printed for those reads but yet reported for remaining reads with single base. Induces bias mono
    # nchar(sample_list[[sample]]$allele_maj)<2 : I think not needed but not a problem.
    
    
    sample_list_bi[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='bi',]
    
    sample_list_bi[[sample]][["chrom"]] <- factor(x = sample_list_bi[[sample]][["chrom"]], 
                                                  levels = chrom_order)
    
    
    sample_list_mono[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='mono',]
    # length(sample_list[[sample]]$pos)==length(unique(sample_list[[sample]]$pos)) #TRUE
    
    sample_list_mono[[sample]][["chrom"]] <- factor(x = sample_list_mono[[sample]][["chrom"]], 
                                                    levels = chrom_order)
    
    df_plot=rbind(sample_list_bi[[sample]],sample_list_mono[[sample]])

    df_plot_mono=sample_list_mono[[sample]]
    # df_plot_mono$allelism=rep("mono",nrow(sample_list_mono[[sample]]))
    
    df_plot_bi=sample_list_bi[[sample]]
    # df_plot_bi$allelism=rep("bi",nrow(sample_list_bi[[sample]]))
    
    ggplot() + 
      # base rectangles for the chroms, with numeric value for each chrom on the x-axis
      statebins:::geom_rrect(data = chrom_sizes, aes(xmin = as.numeric(chromosome) - 0.1, 
                                                     xmax = as.numeric(chromosome) + 0.1, 
                                                     ymax = size, ymin = 0), 
                             colour="black", fill = "white") + 
      # rotate the plot 90 degrees
      coord_flip() +
      # black & white color theme 
      theme(axis.text.x = element_text(colour = "black"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank()) + 
      # give the appearance of a discrete axis with chrom labels
      scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
      # add bands for centromeres
      geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.1, 
                                        xmax = as.numeric(chromosome) + 0.1, 
                                        ymax = end, ymin = start)) +
      # add bands for CNA value
      
      geom_rect(data = df_plot, aes(xmin = as.numeric(chrom) - 0.097,
                                    xmax = as.numeric(chrom) + 0.097,
                                    ymax = pos+150000, ymin = pos-150000, fill=allelism)) +
      
      geom_rect(data = df_plot_mono, aes(xmin = as.numeric(chrom) - 0.097,
                                         xmax = as.numeric(chrom) + 0.097,
                                         ymax = pos+150000, ymin = pos-150000), fill='#3beae4') +
      
      geom_rect(data = df_plot_bi, aes(xmin = as.numeric(chrom) - 0.097, 
                                       xmax = as.numeric(chrom) + 0.097, 
                                       ymax = pos+150000, ymin = pos-150000), fill='#e75b4e') + 
      ylab("region (bp)") + 
      # scale_color_gradient(low="gold", high="red", limits = c(-10,10)) +
      # scale_fill_gradient(low="gold", high="red", limits = c(-10,10)) +
      ggtitle(paste(sep="_",sample_annot[sample, c("Sample_Name","cell_type")], collapse ="_"))
    
    # ggsave(paste0("dataset_2/",paste(sep="_",paste(sep="_",sample_annot[sample, c("Sample_Name","cell_type")], collapse ="_"),"mono_bi_snp_heterozygous_h9_samtools_pileup_chrX.pdf")))

    
    # % biallelism chrX primed vs naive
    percent_biallelism[[sample]]=nrow(df_plot_bi)/(nrow(df_plot_bi)+nrow(df_plot_mono))*100
    print(sample)
}

sample=NULL

for(sample in names(sample_list)){
  df_to_be_written=sample_list[[sample]][,c("chrom","pos","dp","allele_maj","allele_maj_num","allele_min","allele_min_num","allelism")]
  write.tsv(df_to_be_written,file=paste0("dataset_2/",paste(sep = "_",sample,"_chrX_snp_table.tsv")),row.names = F)
}

# CHR7

open_file <- function(sample){
  snp <- read.table(paste0("dataset_2/vcf/", sample, "_samtools_mpileup_informative_from_wgs_h9_chr7.pileup"), header = FALSE,
                    sep = "\t")
  colnames(snp) <- c("chrom", "pos", "ref", "dp", "alts", "bq", "mq")
  return(snp)
}

gtf_chr7=gtf_unique[gtf_unique$chrom=='chr7',]

chrom_sizes <- structure(list(chromosome = c("chr7"), size = c(159345973L)), .Names = c("chromosome", 
                                                                                        "size"), class = "data.frame", row.names = c(NA, -1L))
chrom_sizes

# hg38 centromere locations
centromeres <- structure(list(chromosome = c("chr7"), start = c(58054331L), end = c(61054331L)), 
                         .Names = c("chromosome", "start", "end"), class = "data.frame", row.names = c(NA, -1L
                         ))
centromeres

# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chr7")
chrom_key <- setNames(object = as.character(c(1)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                      levels = chrom_order)

centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                      levels = chrom_order)
# create a color key for the plot
#group.colors <- c(gain = "red", loss = "blue")

informative_pos_vcf_chr7=read.table('C:/Users/gael/charbel_paper/rna/charbel_paper_ifb/H9_WGS_hg38_filtered_chr7.vcf')
colnames(informative_pos_vcf_chr7)=c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample")
informative_pos_vcf_chr7=informative_pos_vcf_chr7[grep("0/0|0/1|1/1", informative_pos_vcf_chr7$sample),] #remove alt_2 because very low frequency/probability
# and annoying for computation
informative_pos_vcf_chr7=informative_pos_vcf_chr7[informative_pos_vcf_chr7$chrom=="chr7" &
                                                    informative_pos_vcf_chr7$qual >=50,]
informative_pos_vcf_chr7$GT=sapply(str_split(pattern=":", informative_pos_vcf_chr7$sample),function(x) x[1])
informative_pos_vcf_chr7$AD=sapply(str_split(pattern=":", informative_pos_vcf_chr7$sample),function(x) x[2])
informative_pos_vcf_chr7$DP=as.numeric(sapply(str_split(pattern=":", informative_pos_vcf_chr7$sample),function(x) x[3]))
informative_pos_vcf_chr7$GQ=as.numeric(sapply(str_split(pattern=":", informative_pos_vcf_chr7$sample),function(x) x[4]))
informative_pos_vcf_chr7$PL=sapply(str_split(pattern=":", informative_pos_vcf_chr7$sample),function(x) x[5])
informative_pos_vcf_chr7=informative_pos_vcf_chr7[informative_pos_vcf_chr7$DP>=10,]
alleles=lapply(lapply(str_split(pattern=",", informative_pos_vcf_chr7$AD),function(x) x), as.numeric)
alleles=do.call('rbind',alleles)
informative_pos_vcf_chr7$allele_0=alleles[,1]
informative_pos_vcf_chr7$allele_1=alleles[,2]
informative_pos_vcf_chr7$allele_maj_num=apply(informative_pos_vcf_chr7,1,function(x){max(as.numeric(x['allele_0']),as.numeric(x['allele_1']))})
informative_pos_vcf_chr7$allele_min_num=apply(informative_pos_vcf_chr7,1,function(x){min(as.numeric(x['allele_0']),as.numeric(x['allele_1']))})
informative_pos_vcf_chr7=detect_bi(informative_pos_vcf_chr7)

sample_annot=read.csv("dataset_2/SampleSheet.csv", header=T, comment.char = '#')
rownames(sample_annot)=paste0("D1507",sample_annot$Sample_ID)
sample_annot$sample_id=paste0("D1507",sample_annot$Sample_ID)
sample_annot$group=sub('_[^_]*$', '', sample_annot$Sample_Name)
sample_annot$group=str_replace_all(sample_annot$group,"CTL","1")
sample_annot$group=str_replace_all(sample_annot$group,"Mix","2")
sample_annot$replicate=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$replicate[sample_annot$replicate=="primed"]<-"Rep1"
sample_annot$cell_type=sub('.*_', '', sample_annot$Sample_Name)
sample_annot$cell_type[sample_annot$cell_type!="primed"]<-"pxgl"

sample_list=NULL
sample_list_bi=NULL
sample_list_mono=NULL
percent_biallelism=NULL
alleles=NULL
alts_list=NULL
strand_list=NULL
bq_list=NULL
mq_list=NULL

# ! different VCF format from bcftools. e.g:Number of 1) forward ref alleles; 2) reverse ref; 3) forward alt; 4) reverse alt alleles

for (sample in rn(sample_annot)){
  
  sample_list[[sample]]=open_file(sample)
  # nrow(sample_list[[sample]])==nrow(sample_list[[sample]][sample_list[[sample]]$pos %in% informative_pos_vcf_chr7$pos,]) # TRUE ; no need, because samtools pileup performed on these pos
  sample_list[[sample]]=sample_list[[sample]][sample_list[[sample]]$dp>=10,]
  
  alts_list=strsplit(sample_list[[sample]]$alts, "")
  
  alts_list=lapply(alts_list,function(x){x[!(x%in%c('*','#'))]}) #remove deletions
  
  strand_list=lapply(alts_list, function(x){
    
    x[str_detect(x,"[[:upper:]]")]<-"+"
    x[str_detect(x,"[[:lower:]]")]<-"-"
    x[x=="."]<-"+"
    x[x==","]<-"-"
    
  })
  
  for (i in seq(length(alts_list))){
    
    alts_list[[i]][alts_list[[i]]%in%c(".",",")]<-sample_list[[sample]][i,'ref']
    
  }  
  
  alts_list=lapply(alts_list,function(x){toupper(x)})
  
  alts_list=lapply(alts_list, function(x){
    
    if(dim(table(x))>2){x[x%in%(names(table(x)[order(table(x),decreasing = T)][1:2]))]}else{x}
    
  }) # remove putative sequencing errors (more than 2 base types detected)
  
  sample_list[[sample]]$dp=unlist(lapply(alts_list,function(x){length(x)}))
  
  sample_list[[sample]]$allele_maj_num=unlist(lapply(alts_list,function(x){table(x)[order(table(x),decreasing = T)][1]}))
  sample_list[[sample]]$allele_maj=names(unlist(lapply(alts_list,function(x){table(x)[order(table(x),decreasing = T)][1]})))
  sample_list[[sample]]$allele_min_num=unlist(lapply(alts_list,function(x){tail(table(x)[table(x)==min(table(x))],n=1)}))
  sample_list[[sample]]$allele_min=names(unlist(lapply(alts_list,function(x){tail(table(x)[table(x)==min(table(x))],n=1)})))
  
  sample_list[[sample]]$ref=toupper(sample_list[[sample]]$ref)
  
  informative_pos_vcf_chr7_subset=informative_pos_vcf_chr7[informative_pos_vcf_chr7$pos%in%sample_list[[sample]]$pos,]
  # informative_pos_vcf_chr7_subset$pos==sample_list[[sample]]$pos
  # all(informative_pos_vcf_chr7_subset$pos==sample_list[[sample]]$pos) # because ordered
  
  # CAUTION: avoid this situation, from putative sequencing error:
  # "chrom" "pos" "ref" "dp" "alts" "bq" "mq" "allele_maj_num" "allele_maj" "allele_min_num" "allele_min" "allelism"      
  # chr7 15702622 G 13 ,,,....,,,,., FFFFFkFFFFFFk 51 13 G 0 A mono
  # chr7 15702622 G 12 ,.,.,t,,,,., FkFFFFFFFFFF 52 11 G 1 T mono
  # NOT same allele_min. Take the one from the info_chr7 that is different from allele_maj
  
  for (i in seq(nrow(sample_list[[sample]]))){
    
    if(sample_list[[sample]]$allele_min[i]==sample_list[[sample]]$allele_maj[i] |
       !(sample_list[[sample]]$allele_min[i]%in%c(informative_pos_vcf_chr7_subset$ref[i],
                                                  informative_pos_vcf_chr7_subset$alt[i]))){
      
      sample_list[[sample]]$allele_min[i]=c(informative_pos_vcf_chr7_subset$ref[i],informative_pos_vcf_chr7_subset$alt[i])[
        c(informative_pos_vcf_chr7_subset$ref[i],informative_pos_vcf_chr7_subset$alt[i])!=sample_list[[sample]]$allele_maj[i]][1]
      
      
      sample_list[[sample]]$allele_min_num[i]=0
      
    }else{
      
      sample_list[[sample]]$allele_min[i]=sample_list[[sample]]$allele_min[i]
      sample_list[[sample]]$allele_min_num[i]=sample_list[[sample]]$allele_min_num[i]
      
    }
    
  }
  
  mq_list=strsplit(sample_list[[sample]]$mq, "") # convert ascii bq & mq to numeric values, and compute median
  sample_list[[sample]]$mq=unlist(lapply(mq_list,function(x){median(qual_conv(x),na.rm=T)}))
  sample_list[[sample]]=sample_list[[sample]][sample_list[[sample]]$mq>=50,]
  
  sample_list[[sample]]$allelism=rep(NA,nrow(sample_list[[sample]]))
  sample_list[[sample]]$allelism[sample_list[[sample]]$pos%in%detect_bi(sample_list[[sample]])$pos]<-"bi"
  sample_list[[sample]]$allelism[!(sample_list[[sample]]$pos%in%detect_bi(sample_list[[sample]])$pos)]<-"mono"
  
  sample_list[[sample]]=sample_list[[sample]][nchar(sample_list[[sample]]$allele_min)<2 & nchar(sample_list[[sample]]$allele_maj)<2,] # remove indels positions that were not printed for those reads but yet reported for remaining reads with single base. Induces bias mono
  # nchar(sample_list[[sample]]$allele_maj)<2 : I think not needed but not a problem.
  
  
  sample_list_bi[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='bi',]
  
  sample_list_bi[[sample]][["chrom"]] <- factor(x = sample_list_bi[[sample]][["chrom"]], 
                                                levels = chrom_order)
  
  
  sample_list_mono[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='mono',]
  # length(sample_list[[sample]]$pos)==length(unique(sample_list[[sample]]$pos)) #TRUE
  
  sample_list_mono[[sample]][["chrom"]] <- factor(x = sample_list_mono[[sample]][["chrom"]], 
                                                  levels = chrom_order)
  
  df_plot=rbind(sample_list_bi[[sample]],sample_list_mono[[sample]])
  
  df_plot_mono=sample_list_mono[[sample]]
  # df_plot_mono$allelism=rep("mono",nrow(sample_list_mono[[sample]]))
  
  df_plot_bi=sample_list_bi[[sample]]
  # df_plot_bi$allelism=rep("bi",nrow(sample_list_bi[[sample]]))
  
  ggplot() + 
    # base rectangles for the chroms, with numeric value for each chrom on the x-axis
    statebins:::geom_rrect(data = chrom_sizes, aes(xmin = as.numeric(chromosome) - 0.1, 
                                                   xmax = as.numeric(chromosome) + 0.1, 
                                                   ymax = size, ymin = 0), 
                           colour="black", fill = "white") + 
    # rotate the plot 90 degrees
    coord_flip() +
    # black & white color theme 
    theme(axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) + 
    # give the appearance of a discrete axis with chrom labels
    scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
    # add bands for centromeres
    geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.1, 
                                      xmax = as.numeric(chromosome) + 0.1, 
                                      ymax = end, ymin = start)) +
    # add bands for CNA value
    
    geom_rect(data = df_plot, aes(xmin = as.numeric(chrom) - 0.097,
                                  xmax = as.numeric(chrom) + 0.097,
                                  ymax = pos+150000, ymin = pos-150000, fill=allelism)) +
    
    geom_rect(data = df_plot_mono, aes(xmin = as.numeric(chrom) - 0.097,
                                       xmax = as.numeric(chrom) + 0.097,
                                       ymax = pos+150000, ymin = pos-150000), fill='#3beae4') +
    
    geom_rect(data = df_plot_bi, aes(xmin = as.numeric(chrom) - 0.097, 
                                     xmax = as.numeric(chrom) + 0.097, 
                                     ymax = pos+150000, ymin = pos-150000), fill='#e75b4e') + 
    ylab("region (bp)") + 
    # scale_color_gradient(low="gold", high="red", limits = c(-10,10)) +
    # scale_fill_gradient(low="gold", high="red", limits = c(-10,10)) +
    ggtitle(paste(sep="_",sample_annot[sample, c("Sample_Name","cell_type")], collapse ="_"))
  
  # ggsave(paste0("dataset_2/",paste(sep="_",paste(sep="_",sample_annot[sample, c("Sample_Name","cell_type")], collapse ="_"),"mono_bi_snp_heterozygous_h9_samtools_pileup_chr7.pdf")))
  
  
  # % biallelism chr7 primed vs naive
  percent_biallelism[[sample]]=nrow(df_plot_bi)/(nrow(df_plot_bi)+nrow(df_plot_mono))*100
  print(sample)
}

sample=NULL

for(sample in names(sample_list)){
  df_to_be_written=sample_list[[sample]][,c("chrom","pos","dp","allele_maj","allele_maj_num","allele_min","allele_min_num","allelism")]
  write.tsv(df_to_be_written,file=paste0("dataset_2/",paste(sep = "_",sample,"_chr7_snp_table.tsv")),row.names = F)
}
