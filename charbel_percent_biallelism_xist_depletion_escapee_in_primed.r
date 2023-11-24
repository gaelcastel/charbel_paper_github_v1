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
# library(biomaRt)
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


constant_genes=read.table("constant.gtf", header=F, sep="\t")
colnames(constant_genes)=c('chrom','source','feature','start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute')

up_genes=read.table("up.gtf", header=F, sep="\t")
colnames(up_genes)=c('chrom','source','feature','start',
                           'end',
                           'score',
                           'strand',
                           'frame',
                           'attribute')

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
sample_annot_xist_depletion=sample_annot[!(sample_annot$dataset%in%c("rna_seq_d165","rna_seq_dvallot")),]

bi_percent_df=data.frame(row.names=names(percent_biallelism[rn(sample_annot_xist_depletion)]),
                         perbi=as.numeric(percent_biallelism[rn(sample_annot_xist_depletion)]),
                         cell_type_condition=sample_annot_xist_depletion$cell_type_condition,
                         dataset=sample_annot_xist_depletion$dataset)


bi_percent_df$dataset[grep("gro_seq",bi_percent_df$dataset)]<-"GRO-Seq DOX -/+"
bi_percent_df$dataset[bi_percent_df$dataset=="rna_seq_d1269"]<-"DOX -/+"
bi_percent_df$dataset[bi_percent_df$dataset=="rna_seq_d1269_bis"]<-"WT / KO"

bi_percent_df$dataset=factor(bi_percent_df$dataset,levels=c("GRO-Seq DOX -/+",
                                                            "DOX -/+",
                                                            "WT / KO"))

ggplot(bi_percent_df, aes(x=cell_type_condition, y=perbi, fill=cell_type_condition)) +
  geom_boxplot(lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  geom_point(size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  facet_wrap(~dataset, scales="free", nrow=1, ncol=3)+
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
  stat_compare_means(vjust=1,label="p.signif",size=3, method="t.test", paired=FALSE,label.y = 102)  

# ggsave("percentage_biallelism_samtools_mpileup_xist_depletion.pdf")

# dev.off()



# charbel_primed_escapees_inactivated top snps per gene

# ADD EXON/INTRON INFO

gtf_chrX_filtered=gtf_chrX[-grep('_dup',gtf_chrX$attribute),] # transcript level
charbel_genes=sub(';.*','',sub('.*gene_id ','',gtf_chrX_filtered$feature))

gtf_x=gtf[gtf$chrom=="chrX",] # exon level
gtf_x=gtf_x[-grep('_dup',gtf_x$attribute),]
gtf_x$gene=sub(';.*','',sub('.*gene_id ','',gtf_x$attribute))


sample_annot_primed=sample_annot[sample_annot$cell_type_condition=="primed_control",]

sample_list_genes=NULL

for(sample in rn(sample_annot_primed)){
  
  
  sample_list_genes[[sample]]=apply(gtf_chrX_filtered,1,function(x){
    
    return(sample_list[[sample]][sample_list[[sample]]$pos > as.numeric(x['start']) & sample_list[[sample]]$pos < as.numeric(x['end']),])
    
  })
  names(sample_list_genes[[sample]])=charbel_genes
  
  sample_list_genes[[sample]]=do.call('rbind',sample_list_genes[[sample]])
  sample_list_genes[[sample]]$gene=sub('\\..*','',rn(sample_list_genes[[sample]]))
  
  sample_list_genes[[sample]] %>%
    group_by(gene) %>%
    filter(max(dp)/dp<2) -> sample_list_genes[[sample]] # select 'TOP' SNP PER-GENE
  
  sample_list_genes[[sample]]$feature=apply(sample_list_genes[[sample]],1,function(x){
    
    gtf_x_subset=gtf_x[gtf_x$gene==x['gene'] & gtf_x$feature=='exon',]
    
    test_exon_intron=data.table::between(as.numeric(x['pos']),as.numeric(gtf_x_subset$start),as.numeric(gtf_x_subset$end))
    
    if(any(test_exon_intron)){return('exon')}else{return('intron')}
    
  })

}

# frequency of allelism pos to determine majo allelism per gene

primed_final_allelism=lapply(sample_list_genes,function(x){
  
  x %>%
    group_by(gene) %>%
    filter(allelism=='mono') %>%
    summarise(mono=n()) -> x_mono
    
  x %>%
    group_by(gene) %>%
    filter(allelism=='bi') %>%
    summarise(bi=n()) -> x_bi
  
  x_mono_bi=merge(x_mono,x_bi,by="gene",all=T)
  x_mono_bi[is.na(x_mono_bi)]=0
  
  x_allelism=apply(x_mono_bi,1,function(x){
    
    ifelse(as.numeric(x["bi"])>=as.numeric(x["mono"]),"bi","mono")
    
  })
  
  x_allelism=data.frame(gene=x_mono_bi$gene,allelism=x_allelism)
  
  return(x_allelism)
  
})


primed_final_allelism=do.call("rbind",primed_final_allelism)

primed_final_allelism %>%
  group_by(gene) %>%
  summarise(mono=sum(allelism=='mono'),bi=sum(allelism=='bi')) -> primed_final_allelism

primed_mono_bi_genes=apply(primed_final_allelism,1,function(x){
  
  ifelse(as.numeric(x["bi"])>=as.numeric(x["mono"]),"bi","mono")
  
})

names(primed_mono_bi_genes)=primed_final_allelism$gene

# escapees=names(primed_mono_bi_genes[primed_mono_bi_genes=='bi'])
# inactivated=names(primed_mono_bi_genes[primed_mono_bi_genes=='mono'])

up_genes=unique(sub('.*gene_name ','',unlist(lapply(str_split(up_genes$attribute,pattern = ";"),function(x){x[grep("gene_name",x)]}))))
constant_genes=unique(sub('.*gene_name ','',unlist(lapply(str_split(constant_genes$attribute,pattern = ";"),function(x){x[grep("gene_name",x)]}))))

up_genes_in_list=up_genes[up_genes%in%c(names(primed_mono_bi_genes))]
constant_genes_in_list=constant_genes[constant_genes%in%c(names(primed_mono_bi_genes))]

up_genes_status=sapply(up_genes_in_list,function(x){primed_mono_bi_genes[x]})
up_genes_status=data.frame(cat=rep("up",length(up_genes_status)),status=up_genes_status)
rownames(up_genes_status)=sub("\\..*","",rn(up_genes_status))
up_genes_status$gene=rn(up_genes_status)

constant_genes_status=sapply(constant_genes_in_list,function(x){primed_mono_bi_genes[x]})
constant_genes_status=data.frame(cat=rep("constant",length(constant_genes_status)),status=constant_genes_status)
rownames(constant_genes_status)=sub("\\..*","",rn(constant_genes_status))
constant_genes_status$gene=rn(constant_genes_status)

up_constant_genes_status=rbind(up_genes_status,constant_genes_status)

up_constant_genes_status$status[up_constant_genes_status$status=="bi"]="escapee"
up_constant_genes_status$status[up_constant_genes_status$status=="mono"]="inactivated"

pval <- stats::chisq.test(up_constant_genes_status$cat,up_constant_genes_status$status,correct = T,simulate.p.value = TRUE)$p.value
# CAUTION !!! PVAL CHANGES (BUT WITHIN GIVEN RANGE), AS IT IS COMPUTED FROM SIMULATION TO "INCREASE" N OBS

up_constant_genes_status%>%
  ggplot(aes(x = cat, y = frequency(status), fill = status)) +
  geom_bar(stat = "identity",position = "fill")+
  ggtitle('Frequency of escapees among gene categories\n(XIST depletion in naive)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, by = 0.1))+
  annotate("text", x=1.5, y=1.1, label=round(pval,5))

# ggsave("frequency_of_escapees_gene_categories.pdf")

table(up_constant_genes_status$status[up_constant_genes_status$cat=="up"])
table(up_constant_genes_status$status[up_constant_genes_status$cat=="constant"])

# or test hypergeometric test? With phyper?