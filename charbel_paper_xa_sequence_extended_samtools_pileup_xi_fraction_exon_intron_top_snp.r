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
    snp <- read.table(paste0("vcf/", sample, "_samtools_mpileup_informative_from_wgs_h9_chrX.pileup"), header = FALSE,
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

sample_annot=fastRead("sample_annot_all.txt", header=T, sep=",",as.matrix = F)
sample_annot$sample_id=rn(sample_annot)
sample_annot$cell_type_condition=paste(sep="_",sample_annot$cell_type,sample_annot$condition)

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
      ggtitle(paste(sep="_",sample_annot[sample, c("Sample_Name","cell_type","condition")], collapse ="_"))
    
    # ggsave(paste(sep="_",paste(sep="_",sample_annot[sample, c("Sample_Name","cell_type","condition")], collapse ="_"),"mono_bi_snp_heterozygous_h9_samtools_pileup_chrX.pdf"))

    
    # % biallelism chrX primed vs naive
    percent_biallelism[[sample]]=nrow(df_plot_bi)/(nrow(df_plot_bi)+nrow(df_plot_mono))*100
    print(sample)
}


percent_biallelism=unlist(percent_biallelism)
sample_annot_naive_primed=sample_annot[sample_annot$cell_type_condition%in%c("primed_control","naive_control") & !(grepl("gro_seq",sample_annot$dataset)),]

bi_percent_df=data.frame(row.names=names(percent_biallelism[rn(sample_annot_naive_primed)]),
                         perbi=as.numeric(percent_biallelism[rn(sample_annot_naive_primed)]),
                         cell_type=sample_annot_naive_primed$cell_type)


ggplot(bi_percent_df, aes(x=cell_type, y=perbi, fill=cell_type)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  geom_point(size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  # facet_wrap(~genes, scales="free", nrow=6, ncol=3)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  scale_x_discrete(expand = c(0, 0.5)) +
  ylab("% biallelism") +
  # xlab("Cell Type")+
  geom_blank(aes(y = 0)) + #y_axis start at 0
  # scale_fill_manual("Cell Type :", values = colors)+
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, linetype = 'dashed',colour = "grey"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  stat_boxplot(geom='errorbar', linetype=1, width=0.8, size=0.25)+
  # scale_y_continuous(trans="pseudo_log",breaks =c(0,5,10,50,100), limits=c(0,max(dataDeSeqForPlot$value)+100))+
  stat_compare_means(vjust=1,label="p.signif",size=3, method="wilcox.test", paired=FALSE,label.y = 102)  

# ggsave("percentage_biallelism_naive_primed_25_percent_thres_dp_10_samtools_mpileup.pdf")

xa_xi_rep=read.tsv("h9_xa_xi_sequence.tsv")

# View(xa_xi_rep[xa_xi_rep$pos==53614686,])
# View(informative_pos_vcf_chrX_subset[informative_pos_vcf_chrX_subset$pos==53614686,]) # We do have this HUWE1 exonic pos now

sample_list_filtered=NULL
xi_fraction_list=NULL

for (sample in rn(sample_annot)){
  
  
  sample_list_filtered[[sample]]=sample_list[[sample]][sample_list[[sample]]$pos%in%xa_xi_rep$pos,]
  xa_xi_rep_filtered=xa_xi_rep[xa_xi_rep$pos%in%sample_list_filtered[[sample]]$pos,]
  sample_list_filtered[[sample]]=sample_list_filtered[[sample]][match(xa_xi_rep_filtered$pos,sample_list_filtered[[sample]]$pos),]
  # sample_list_filtered[[sample]]$pos==xa_xi_rep_filtered$pos #TRUE
  
  sample_list_filtered[[sample]]$xi_fraction=NULL
  
  for (i in seq(nrow(sample_list_filtered[[sample]]))){
    
    if(sample_list_filtered[[sample]]$allele_min[i]==xa_xi_rep_filtered$xi_final[i]){
      
      sample_list_filtered[[sample]]$xi_fraction[i]=sample_list_filtered[[sample]]$allele_min_num[i]/
        (sample_list_filtered[[sample]]$allele_maj_num[i]+sample_list_filtered[[sample]]$allele_min_num[i]) #allele_maj cannot be 0 due to our settings (dp>=10, so allele_maj_num >=5)
      
    }else if(sample_list_filtered[[sample]]$allele_min[i]==xa_xi_rep_filtered$xa_final[i]){
      
      sample_list_filtered[[sample]]$xi_fraction[i]=sample_list_filtered[[sample]]$allele_maj_num[i]/
        (sample_list_filtered[[sample]]$allele_maj_num[i]+sample_list_filtered[[sample]]$allele_min_num[i]) #allele_maj cannot be 0 due to our settings (dp>=10, so allele_maj_num >=5)
      
    }else{sample_list_filtered[[sample]]$xi_fraction[i]=NA}
  } # allele_min[1] != Xa[i] & != Xi[i]

  xi_fraction_list[[sample]]=mean(sample_list_filtered[[sample]]$xi_fraction, na.rm=T)
  
}

xi_fraction_vector=unlist(xi_fraction_list)

# all(names(xi_fraction_vector)==rn(sample_annot)) # TRUE

sample_annot_charbel=sample_annot[!(grepl('vallot',sample_annot$dataset)),]
xi_fraction_vector_charbel=xi_fraction_vector[rn(sample_annot_charbel)]

df_plot=NULL
df_plot$sample=rn(sample_annot_charbel)
df_plot$cell_type=sample_annot_charbel$cell_type
df_plot$cell_type_condition_dataset=paste(sep="_",sub('_d.*','',sample_annot_charbel$dataset),sample_annot_charbel$cell_type_condition)
df_plot$xi_fraction=xi_fraction_vector_charbel
df_plot$sample_name=sample_annot_charbel$Sample_Name
df_plot$cell_type_condition_dataset[grep('PXGL_WT',df_plot$sample_name)]<- 'rna_seq_naive_wt'
df_plot$cell_type_condition_dataset[grep('PXGL_KO',df_plot$sample_name)]<- 'rna_seq_naive_ko'

df_plot=as.data.frame(df_plot)

df_plot %>%
  group_by(cell_type_condition_dataset) %>%
  mutate(sd=sd(xi_fraction), mean=mean(xi_fraction)) %>%
  as.data.frame() -> df_plot


df_plot$cell_type_condition_dataset=factor(df_plot$cell_type_condition_dataset,levels = c("rna_seq_naive_cult_control",
                                                                                          "rna_seq_primed_control",
                                                                                          'rna_seq_naive_wt',
                                                                                          'rna_seq_naive_ko',
                                                                                          "rna_seq_naive_control",
                                                                                          "rna_seq_naive_treatment",
                                                                                          "gro_seq_naive_control",
                                                                                          "gro_seq_naive_treatment"
                                                                                          # ,
                                                                                          # "rna_seq_vallot_primed_non_eroded",
                                                                                          # "rna_seq_vallot_primed_eroded"
                                                                                          )) # to place in appropriate order on the plot





my_comparisons=list(c("rna_seq_naive_cult_control","rna_seq_primed_control"),c('rna_seq_naive_wt','rna_seq_naive_ko'),
                    c("rna_seq_naive_control","rna_seq_naive_treatment"),c("gro_seq_naive_control","gro_seq_naive_treatment"))

ggplot(data=df_plot,aes(x=cell_type_condition_dataset,y=xi_fraction,fill=cell_type_condition_dataset))+
  geom_bar(stat = 'summary')+
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position="dodge") +
  stat_compare_means(comparisons = my_comparisons, vjust=0.2, label="p.signif",size=3, method="t.test", paired=F)  +
  stat_compare_means(comparisons = my_comparisons, vjust=1, label="p.format",size=3, method="t.test", paired=F)  + 
  ggtitle('Xi fraction') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave("charbel_xi_fraction_samtools_pileup.pdf")


# charbel_snp_all_rep

gtf_chrX_filtered=gtf_chrX[-grep('_dup',gtf_chrX$attribute),]

charbel_genes=sub(';.*','',sub('.*gene_id ','',gtf_chrX_filtered$feature))

sample_list_gene_info=NULL
sample_list_gene_fraction=NULL

for (sample in rn(sample_annot_charbel)){
  
  sample_list_gene_info[[sample]]=apply(gtf_chrX_filtered,1,function(x){
    
    return(sample_list_filtered[[sample]][sample_list_filtered[[sample]]$pos > as.numeric(x['start']) & sample_list_filtered[[sample]]$pos < as.numeric(x['end']),])
    
  })
  names(sample_list_gene_info[[sample]])=charbel_genes
  
  sample_list_gene_info[[sample]]=do.call('rbind',sample_list_gene_info[[sample]])
  sample_list_gene_info[[sample]]$gene=sub('\\..*','',rn(sample_list_gene_info[[sample]]))
  
  # xa_xi_subset=xa_xi[xa_xi$pos%in%sample_list_gene_info[[sample]]$pos,]
  
  sample_list_gene_info[[sample]] %>%
    summarise(gene=gene,pos=pos,xi_fraction=xi_fraction,dp=dp) -> sample_list_gene_fraction[[sample]]
  
  
  sample_list_gene_fraction[[sample]]$sample=rep(sample,nrow(sample_list_gene_fraction[[sample]]))
  sample_list_gene_fraction[[sample]]$sample_name=rep(sample_annot_charbel[sample,'Sample_Name'],nrow(sample_list_gene_fraction[[sample]]))
  sample_list_gene_fraction[[sample]]$cell_type_condition_dataset=rep(paste(sep="_",sub('_d.*','',sample_annot_charbel[sample,'dataset']),sample_annot_charbel[sample,'cell_type_condition']),nrow(sample_list_gene_fraction[[sample]]))
  sample_list_gene_fraction[[sample]]$cell_type=rep(sample_annot_charbel[sample,'cell_type'],nrow(sample_list_gene_fraction[[sample]]))
  sample_list_gene_fraction[[sample]]$condition=rep(sample_annot_charbel[sample,'condition'],nrow(sample_list_gene_fraction[[sample]]))
}


sample_list_gene_fraction=do.call('rbind',sample_list_gene_fraction)

sample_list_gene_fraction$cell_type_condition_dataset[grep('PXGL_WT',sample_list_gene_fraction$sample_name)]<- 'rna_seq_naive_wt'
sample_list_gene_fraction$cell_type_condition_dataset[grep('PXGL_KO',sample_list_gene_fraction$sample_name)]<- 'rna_seq_naive_ko'
sample_list_gene_fraction$condition[grep('PXGL_WT',sample_list_gene_fraction$sample_name)]<- 'wt'
sample_list_gene_fraction$condition[grep('PXGL_KO',sample_list_gene_fraction$sample_name)]<- 'ko'

sample_list_gene_fraction$cell_type_condition_dataset=factor(sample_list_gene_fraction$cell_type_condition_dataset,levels = c("rna_seq_naive_cult_control",
                                                                                                                            "rna_seq_primed_control",
                                                                                                                            "rna_seq_naive_wt",
                                                                                                                            "rna_seq_naive_ko",
                                                                                                                            "rna_seq_naive_control",
                                                                                                                            "rna_seq_naive_treatment",
                                                                                                                            "gro_seq_naive_control",
                                                                                                                            "gro_seq_naive_treatment"))

sample_list_gene_fraction$group=rep(1,nrow(sample_list_gene_fraction))
sample_list_gene_fraction$group[sample_list_gene_fraction$cell_type_condition_dataset%in%c("rna_seq_naive_wt","rna_seq_naive_ko")] <- 1
sample_list_gene_fraction$group[sample_list_gene_fraction$cell_type_condition_dataset%in%c("rna_seq_naive_control","rna_seq_naive_treatment")] <- 2
sample_list_gene_fraction$group[sample_list_gene_fraction$cell_type_condition_dataset%in%c("gro_seq_naive_control","gro_seq_naive_treatment")] <- 3
sample_list_gene_fraction$group[sample_list_gene_fraction$cell_type_condition_dataset%in%c("rna_seq_naive_cult_control","rna_seq_primed_control")] <- 4


sample_list_gene_fraction$gene_snp=paste(sep = '_',sample_list_gene_fraction$gene,sample_list_gene_fraction$pos)

# select genes with SNPs detected in all replicates

sample_list_gene_fraction_filtered=sample_list_gene_fraction[!(grepl('GRO_seq_B4_C',sample_list_gene_fraction$sample_name)),] #remove fucked-up GRO_seq_B4_C

sample_list_gene_fraction_filtered %>%
  group_by(group,gene_snp) %>%
  filter(group!=3,n()>=6) %>%
  arrange(gene) -> sample_list_snp_all_rep_rna_seq

sample_list_gene_fraction_filtered %>%
  group_by(group,gene_snp) %>%
  filter(group==3,n()>=5) %>%
  arrange(gene) -> sample_list_snp_all_rep_gro_seq

sample_list_snp_all_rep=rbind(sample_list_snp_all_rep_rna_seq,
                              sample_list_snp_all_rep_gro_seq)

sample_list_snp_all_rep %>%
  group_by(cell_type_condition_dataset,gene_snp) %>%
  summarise(xi_mean_fraction=mean(xi_fraction),group=unique(group), snp=unique(pos),dp_mean=round(mean(dp))) %>%
  as.data.frame()-> sample_list_snp_all_rep_mean

sample_list_snp_all_rep_mean %>%
  group_split(group)  -> snp_all_rep_group

xi_fraction_mean=lapply(snp_all_rep_group,function(x){
  
  x %>%
    group_by(gene_snp) %>%
    filter(n()>1) %>% # filter genes observed in both control and test conditions. Facultative because already filtered above
    summarise(gene_snp=gene_snp,xi_mean_fraction=xi_mean_fraction,snp=snp,cell_type_condition_dataset=cell_type_condition_dataset,dp_mean=dp_mean) %>%
    ungroup() %>%
    group_by(cell_type_condition_dataset, .add = TRUE) %>%
    group_split() -> xi_fraction_mean_split
  
  x=do.call('cbind',xi_fraction_mean_split)
  colnames(x)=c('gene_snp',paste(sep = "_",cn(x)[2:5],unique(x[,4])),paste(sep = "_",cn(x)[6:10],unique(x[,9])))
  x=x[,c(1,2,5,7,10)]
  return(x)
})

names_snp_all_rep_group=lapply(lapply(snp_all_rep_group,function(x){unique(x$cell_type_condition_dataset)}),function(x){paste(sep="_",x[1],x[2])})
names(snp_all_rep_group)=unlist(names_snp_all_rep_group)
names(xi_fraction_mean)=names(snp_all_rep_group)


sample_list_snp_all_rep %>%
  group_by(group) %>%
  group_split() -> sample_list_snp_all_rep_group

xi_pval=lapply(sample_list_snp_all_rep_group,function(x){
  
  x %>%
    group_by(gene_snp) %>%
    group_split() -> pval_list
  
  names(pval_list)=unlist(lapply(pval_list,function(x){unique(x$gene_snp)}))
  
  pval_per_gene_list=unlist(lapply(pval_list,function(x){
    x %>%
      compare_means(formula=xi_fraction~cell_type_condition_dataset,method='t.test') %>%
      summarise(p=p)
  }))
  
  sd_per_gene_list_sd1=unlist(lapply(pval_list,function(x){
    x %>%
      group_by(cell_type_condition_dataset)%>%
      summarise(sd=sd(xi_fraction)) -> y
    
    sd1=as.numeric(y[1,2])
    names(sd1)=y[1,1]
    
    return(sd1)
    
  }))
  
  sd_per_gene_list_sd2=unlist(lapply(pval_list,function(x){
    x %>%
      group_by(cell_type_condition_dataset)%>%
      summarise(sd=sd(xi_fraction)) -> y
    
    sd2=as.numeric(y[2,2])
    names(sd2)=y[2,1]
    
    return(sd2)
    
  }))
  
  #remove genes that were all 0 and thus pval could not be computed
  names(sd_per_gene_list_sd1)=sub('\\..*','',names(sd_per_gene_list_sd1))
  names(sd_per_gene_list_sd2)=sub('\\..*','',names(sd_per_gene_list_sd2))
  
  pval_list=pval_list[sub("\\.p.*","",names(pval_per_gene_list))]
  sd_per_gene_list_sd1=sd_per_gene_list_sd1[sub("\\.p.*","",names(pval_per_gene_list))]
  sd_per_gene_list_sd2=sd_per_gene_list_sd2[sub("\\.p.*","",names(pval_per_gene_list))]
  
  data.frame(row.names =sub("\\.p.*","",names(pval_per_gene_list)),gene_snp=sub("\\.p.*","",names(pval_per_gene_list)),pval=pval_per_gene_list,sd_1=sd_per_gene_list_sd1,sd_2=sd_per_gene_list_sd2)
  
})

names(sample_list_snp_all_rep_group)=names(xi_fraction_mean)
names(xi_pval)=names(sample_list_snp_all_rep_group)

xi_fraction_mean_pval=NULL

for(name in names(xi_fraction_mean)){
  
  df_2_write=merge(xi_fraction_mean[[name]],xi_pval[[name]], by='gene_snp', all=T)
  df_2_write$gene=sub('_.*','',df_2_write$gene_snp)
  df_2_write$snp=as.numeric(sub('.*_','',df_2_write$gene_snp))
  # write.tsv(df_2_write[,c(9,10,2:8)],file = paste(sep="_",name,"xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd.tsv"))
  df_2_write=df_2_write[,-c(9,10)]
  colnames(df_2_write)=c(cn(xi_fraction_mean[[name]]),paste0('pval_',name),paste0('sd1_',name),paste0('sd2_',name))
  xi_fraction_mean_pval[[name]]=df_2_write[,c(cn(xi_fraction_mean[[name]]),paste0('pval_',name),paste0('sd1_',name),paste0('sd2_',name))]
  
}

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=sub('_.*','',gene_snp),snp=sub('.*_','',gene_snp)) %>% dplyr::select(-gene_snp) %>% as.data.frame() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29:30,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_sd.tsv")

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=sub('_.*','',gene_snp),snp=sub('.*_','',gene_snp)) %>% dplyr::select(-gene_snp) %>% as.data.frame() %>% na.omit() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29:30,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_common_sd.tsv")

# CHARBEL PER-GENE "TOP" SNP

# ADD EXON/INTRON INFO
gtf_x=gtf[gtf$chrom=="chrX",]
gtf_x=gtf_x[-grep('_dup',gtf_x$attribute),]
gtf_x$gene=sub(';.*','',sub('.*gene_id ','',gtf_x$attribute))

xi_fraction_mean_pval=NULL
df_1_gene_list=NULL

for(name in names(xi_fraction_mean)){
  
  df_1=merge(xi_fraction_mean[[name]],xi_pval[[name]], by='gene_snp', all=T)
  df_1$gene=sub('_.*','',df_1$gene_snp)
  df_1$snp=as.numeric(sub('.*_','',df_1$gene_snp))
  
  df_2=df_1[,c(1,3,9,10)]
  colnames(df_2)=c('gene_snp','dp','gene','snp')
  df_3=df_1[,c(1,5,9,10)]
  colnames(df_3)=c('gene_snp','dp','gene','snp')

  df_2 %>%
    group_by(gene) %>%
    filter(max(dp)/dp<2) -> df_2 # select 'TOP' SNP PER-GENE
 
  df_3 %>%
    group_by(gene) %>%
    filter(max(dp)/dp<2) -> df_3 # select 'TOP' SNP PER-GENE
  
  top_snp=intersect(df_2$gene_snp,df_3$gene_snp)
  
  df_1=df_1[df_1$gene_snp%in%top_snp,]
  
  df_1$feature=apply(df_1,1,function(x){
    
    gtf_x_subset=gtf_x[gtf_x$gene==x['gene'] & gtf_x$feature=='exon',]
  
    test_exon_intron=data.table::between(as.numeric(x['snp']),as.numeric(gtf_x_subset$start),as.numeric(gtf_x_subset$end))
    
    if(any(test_exon_intron)){return('exon')}else{return('intron')}

  })
  

  df_2_write=df_1
  
  # write.tsv(df_2_write[,c(9,10,11,2:8)],file = paste(sep="_",name,"xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp.tsv"))
  df_2_write$gene_snp=paste(sep="_",df_2_write$gene_snp,df_2_write$feature)
  df_2_write=df_2_write[,-c(9,10,11)]
  colnames(df_2_write)=c(cn(xi_fraction_mean[[name]]),paste0('pval_',name),paste0('sd1_',name),paste0('sd2_',name))
  xi_fraction_mean_pval[[name]]=df_2_write
  
}

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=sapply(gene_snp,function(x){unlist(str_split(x,pattern='_'))[1]}),
                                                                        snp=as.numeric(sapply(gene_snp,function(x){unlist(str_split(x,pattern='_'))[2]})),
                                                                        feature=sapply(gene_snp,function(x){unlist(str_split(x,pattern='_'))[3]})) %>% dplyr::select(-gene_snp) %>% as.data.frame() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29:31,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_sd_top_snp.tsv")

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=sapply(gene_snp,function(x){unlist(str_split(x,pattern='_'))[1]}),
                                                                        snp=as.numeric(sapply(gene_snp,function(x){unlist(str_split(x,pattern='_'))[2]})),
                                                                        feature=sapply(gene_snp,function(x){unlist(str_split(x,pattern='_'))[3]})) %>% dplyr::select(-gene_snp) %>% as.data.frame() %>% na.omit() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29:31,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_common_sd_top_snp.tsv")


# MEAN PER GENE

xi_fraction_mean_pval=NULL
df_1_gene_list=NULL

for(name in names(xi_fraction_mean)){
  
  df_1=merge(xi_fraction_mean[[name]],xi_pval[[name]], by='gene_snp', all=T)
  df_1$gene=sub('_.*','',df_1$gene_snp)
  df_1$snp=as.numeric(sub('.*_','',df_1$gene_snp))
  
  df_2=df_1[,c(1,2,3,6,7,9,10)]
  colnames(df_2)=c('gene_snp','xi_fraction','dp','pval','sd','gene','snp')
  df_3=df_1[,c(1,4,5,6,8,9,10)]
  colnames(df_3)=c('gene_snp','xi_fraction','dp','pval','sd','gene','snp')
  
  df_2 %>%
    group_by(gene) %>%
    filter(max(dp)/dp<2) %>%
    summarise(xi_fraction=mean(xi_fraction,na.rm=T),dp=round(mean(dp,na.rm=T)),
              pval=mean(pval,na.rm=T),sd=mean(sd,na.rm=T),gene=unique(gene))-> df_2
  
  df_3 %>%
    group_by(gene) %>%
    filter(max(dp)/dp<2) %>%
    summarise(xi_fraction=mean(xi_fraction,na.rm=T),dp=round(mean(dp,na.rm=T)),
              pval=mean(pval,na.rm=T),sd=mean(sd,na.rm=T),gene=unique(gene))-> df_3 # select TOP SNP PER-GENE
  
  df_4=merge(df_2,df_3,by=c('gene'),all=T)
  df_4$pval=as.numeric(apply(df_4,1,function(x){mean(c(as.numeric(x['pval.x']),
                                                    as.numeric(x['pval.y'])),
                                                    na.rm=T)}))
    
  
  df_4$pval.x=NULL
  df_4$pval.y=NULL
  
  df_2_write=df_4[,c(1,2,3,5,6,8,4,7)]
  colnames(df_2_write)=c(cn(xi_fraction_mean[[name]]),'pval','sd1','sd2')
  write.tsv(df_2_write,file = paste(sep="_",name,"xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv"))
  colnames(df_2_write)=c(cn(xi_fraction_mean[[name]]),paste0('pval_',name),paste0('sd1_',name),paste0('sd2_',name))
  xi_fraction_mean_pval[[name]]=df_2_write
  
}

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=gene_snp) %>% dplyr::select(-gene_snp) %>% as.data.frame() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_sd_top_snp_mean_gene.tsv")

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=gene_snp) %>% dplyr::select(-gene_snp) %>% as.data.frame() %>% na.omit() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_common_sd_top_snp_mean_gene.tsv")









sample_list_gene_fraction_filtered %>%
  group_by(cell_type_condition_dataset,gene) %>%
  filter(max(dp)/dp<2) %>% # select 'TOP' SNP PER-GENE
  ungroup() %>%
  group_by(group,gene_snp) %>%
  filter(group!=3,n()>=6) %>%
  arrange(gene) -> sample_list_snp_all_rep_rna_seq

sample_list_gene_fraction_filtered %>%
  group_by(cell_type_condition_dataset,gene) %>%
  filter(max(dp)/dp<2) %>% # select 'TOP' SNP PER-GENE
  ungroup() %>%
  group_by(group,gene_snp) %>%
  filter(group==3,n()>=5) %>%
  arrange(gene) -> sample_list_snp_all_rep_gro_seq

sample_list_snp_all_rep=rbind(sample_list_snp_all_rep_rna_seq,
                              sample_list_snp_all_rep_gro_seq)

sample_list_snp_all_rep %>%
  group_by(cell_type_condition_dataset,gene_snp) %>%
  summarise(xi_mean_fraction=mean(xi_fraction),group=unique(group), snp=unique(pos),dp_mean=round(mean(dp))) %>%
  as.data.frame()-> sample_list_snp_all_rep_mean

sample_list_snp_all_rep_mean %>%
  group_split(group)  -> snp_all_rep_group

xi_fraction_mean=lapply(snp_all_rep_group,function(x){
  
  x %>%
    group_by(gene_snp) %>%
    filter(n()>1) %>% # filter genes observed in both control and test conditions. Facultative because already filtered above
    summarise(gene_snp=gene_snp,xi_mean_fraction=xi_mean_fraction,snp=snp,cell_type_condition_dataset=cell_type_condition_dataset,dp_mean=dp_mean) %>%
    ungroup() %>%
    group_by(cell_type_condition_dataset, .add = TRUE) %>%
    group_split() -> xi_fraction_mean_split
  
  x=do.call('cbind',xi_fraction_mean_split)
  colnames(x)=c('gene_snp',paste(sep = "_",cn(x)[2:5],unique(x[,4])),paste(sep = "_",cn(x)[6:10],unique(x[,9])))
  x=x[,c(1,2,5,7,10)]
  return(x)
})

names_snp_all_rep_group=lapply(lapply(snp_all_rep_group,function(x){unique(x$cell_type_condition_dataset)}),function(x){paste(sep="_",x[1],x[2])})
names(snp_all_rep_group)=unlist(names_snp_all_rep_group)
names(xi_fraction_mean)=names(snp_all_rep_group)


sample_list_snp_all_rep %>%
  group_by(group) %>%
  group_split() -> sample_list_snp_all_rep_group

xi_pval=lapply(sample_list_snp_all_rep_group,function(x){
  
  x %>%
    group_by(gene_snp) %>%
    group_split() -> pval_list
  
  names(pval_list)=unlist(lapply(pval_list,function(x){unique(x$gene_snp)}))
  
  pval_per_gene_list=unlist(lapply(pval_list,function(x){
    x %>%
      compare_means(formula=xi_fraction~cell_type_condition_dataset,method='t.test') %>%
      summarise(p=p)
  }))
  
  sd_per_gene_list_sd1=unlist(lapply(pval_list,function(x){
    x %>%
      group_by(cell_type_condition_dataset)%>%
      summarise(sd=sd(xi_fraction)) -> y
    
    sd1=as.numeric(y[1,2])
    names(sd1)=y[1,1]
    
    return(sd1)
    
  }))
  
  sd_per_gene_list_sd2=unlist(lapply(pval_list,function(x){
    x %>%
      group_by(cell_type_condition_dataset)%>%
      summarise(sd=sd(xi_fraction)) -> y
    
    sd2=as.numeric(y[2,2])
    names(sd2)=y[2,1]
    
    return(sd2)
    
  }))
  
  #remove genes that were all 0 and thus pval could not be computed
  names(sd_per_gene_list_sd1)=sub('\\..*','',names(sd_per_gene_list_sd1))
  names(sd_per_gene_list_sd2)=sub('\\..*','',names(sd_per_gene_list_sd2))
  
  pval_list=pval_list[sub("\\.p.*","",names(pval_per_gene_list))]
  sd_per_gene_list_sd1=sd_per_gene_list_sd1[sub("\\.p.*","",names(pval_per_gene_list))]
  sd_per_gene_list_sd2=sd_per_gene_list_sd2[sub("\\.p.*","",names(pval_per_gene_list))]
  
  data.frame(row.names =sub("\\.p.*","",names(pval_per_gene_list)),gene_snp=sub("\\.p.*","",names(pval_per_gene_list)),pval=pval_per_gene_list,sd_1=sd_per_gene_list_sd1,sd_2=sd_per_gene_list_sd2)
  
})

names(sample_list_snp_all_rep_group)=names(xi_fraction_mean)
names(xi_pval)=names(sample_list_snp_all_rep_group)

xi_fraction_mean_pval=NULL

for(name in names(xi_fraction_mean)){
  
  df_2_write=merge(xi_fraction_mean[[name]],xi_pval[[name]], by='gene_snp', all=T)
  df_2_write$gene=sub('_.*','',df_2_write$gene_snp)
  df_2_write$snp=as.numeric(sub('.*_','',df_2_write$gene_snp))
  write.tsv(df_2_write[,c(9,10,2:8)],file = paste(sep="_",name,"xi_fraction_mean_snp_25_percent_thres_dp_10_pval_sd_top_snp.tsv"))
  df_2_write=df_2_write[,-c(9,10)]
  colnames(df_2_write)=c(cn(xi_fraction_mean[[name]]),paste0('pval_',name),paste0('sd1_',name),paste0('sd2_',name))
  xi_fraction_mean_pval[[name]]=df_2_write[,c(cn(xi_fraction_mean[[name]]),paste0('pval_',name),paste0('sd1_',name),paste0('sd2_',name))]
  
}



xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=sub('_.*','',gene_snp),snp=sub('.*_','',gene_snp)) %>% dplyr::select(-gene_snp) %>% as.data.frame() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29:30,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_sd_top_snp.tsv")

xi_fraction_mean_pval %>% reduce(left_join, by = "gene_snp") %>% mutate(gene=sub('_.*','',gene_snp),snp=sub('.*_','',gene_snp)) %>% dplyr::select(-gene_snp) %>% as.data.frame() %>% na.omit() -> xi_fraction_mean_pval_df
xi_fraction_mean_pval_df=xi_fraction_mean_pval_df[,c(29:30,1:28)]
write.tsv(xi_fraction_mean_pval_df,file = "xi_fraction_mean_snp_all_25_percent_thres_dp_10_pval_common_sd_top_snp.tsv")
