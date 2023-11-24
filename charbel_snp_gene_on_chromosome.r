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

up_and_constant=rbind(up_gtf,
                      constant_gtf)

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

up_and_constant[["chrom"]] <- factor(x = up_and_constant[["chrom"]], 
                                              levels = chrom_order)


# pdf("chrX_gene_categories_gene_depletion_chrom_plot.pdf",width=25,height=10)

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


  geom_rect(data = up_and_constant, aes(xmin = as.numeric(chrom) - 0.097,
                                xmax = as.numeric(chrom) + 0.097,
                                ymax = end+150000, ymin = start-150000, fill=category))+
  geom_text(data = up_and_constant, aes(x = as.numeric(chrom) - 0.12,
                                        y = start, label=gene),angle=90, vjust=0.5, hjust=1)+

  xlab("position") +
  # scale_color_gradient(low="gold", high="red", limits = c(-10,10)) +
  # scale_fill_gradient(low="gold", high="red", limits = c(-10,10)) +
  ggtitle("chrX gene categories\nupon XIST depletion")


# dev.off()

table(up_and_constant$category)
