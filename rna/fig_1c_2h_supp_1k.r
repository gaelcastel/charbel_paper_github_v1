setwd('C:/Users/surfi/Desktop/charbel_paper_2024/2024_04/charbel_paper_github_v1/rna/')

source("https://gitlab.univ-nantes.fr/E114424Z/veneR/-/raw/master/loadFun.R")

options(install.packages.compile.from.source = "always")

# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))

needed_func=c("BiocManager",
"dplyr",
"tidyr",
"ggplot2",
"SingleCellExperiment",
"Matrix",
"stringr",
"devtools",
"pcaMethods",
"data.table",
"parallel",
"DESeq2",
"ggpubr",
"tidyverse",
"EnhancedVolcano",
"biomaRt",
"sva",
"gtools",
"statebins",
"wesanderson")

sapply(needed_func,function(x){
    if (x!="BiocManager" & !requireNamespace(x, quietly = T)){
      BiocManager::install(x, dependencies = TRUE, INSTALL_opts = '--no-lock')
    }

    else if (x=="BiocManager" & !requireNamespace(x, quietly = T)){
      install.packages(x,repos="http://cran.r-project.org", dependencies = TRUE, INSTALL_opts = '--no-lock')
    }
  library(x,character.only = TRUE)
}
)

# chrX

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


gtf=read.table("../../Homo_sapiens.GRCh38.90.gtf", header=F, sep="\t")
colnames(gtf)=c('chrom','source','feature','start',
                'end',
                'score',
                'strand',
                'frame',
                'attribute')

gtf=gtf[!grepl('_dup',gtf$attribute),]

gtf %>%
  filter(chrom=="X" & feature=="gene") -> gtf_chrX_gene


sample_annot=fastRead("sample_annot_all.txt", header=T, sep=",",as.matrix = F)
sample_annot=sample_annot[sample_annot$dataset!="rna_seq_d1507",]
sample_annot$sample_id=rn(sample_annot)
sample_annot$cell_type_condition=paste(sep="_",sample_annot$cell_type,sample_annot$condition)

sample_annot$sample_category=paste(sep="_",sample_annot$cell_type_condition,
                                   sample_annot$dox_ko,
                                   sample_annot$dataset)

sample_annot$sample_category=str_replace(sample_annot$sample_category,"_control","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_rna","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_treatment","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_seq","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d1269","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_bis","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d165","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d1271","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_dvallot","_vallot")

sample_annot$sample_category=factor(sample_annot$sample_category,
                                    levels = c("primed_WT",
                                               "primed_non_eroded_WT_vallot",
                                               "primed_eroded_WT_vallot",
                                               "primed_KO_XACT",
                                               "naive_cult_WT",
                                               "naive_cult_KO_XACT",
                                               "naive_WT",
                                               "naive_KO",
                                               "naive_C",
                                               "naive_T",
                                               "naive_C_gro",
                                               "naive_T_gro"))      

sample_annot$group=as.vector(sample_annot$sample_category)
sample_annot$group[sample_annot$group%in%c("primed_WT","primed_KO_XACT")]<-"group1"
sample_annot$group[sample_annot$group%in%c("naive_cult_WT","naive_cult_KO_XACT")]<-"group2"
sample_annot$group[sample_annot$group%in%c("primed_non_eroded_WT_vallot","primed_eroded_WT_vallot")]<-"group3"
sample_annot$group[sample_annot$group%in%c("naive_WT","naive_KO")]<-"group4"
sample_annot$group[sample_annot$group%in%c("naive_C","naive_T")]<-"group5"
sample_annot$group[sample_annot$group%in%c("naive_C_gro","naive_T_gro")]<-"group6"

sample_annot$cell_type_condition=paste(sep="_",sample_annot$cell_type,sample_annot$condition)

ascii_qual=read.tsv("ascii_qual.tsv",comment.char = "")
ascii_qual$symbol=rn(ascii_qual)

informative_pos_vcf_chrX=read.table('H9_WGS_hg38_filtered_chrX.vcf')
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

# # write.tsv(data.frame(informative_pos_vcf_chrX$chrom,
#                      informative_pos_vcf_chrX$pos,
#                      informative_pos_vcf_chrX$pos),file="chrX_heterozygous_snp_h9.bed",row.names = F, col.names = F)

# chromosome_plots_biallelism

# remember to add chr7 as a control

chrom_sizes <- structure(list(chromosome = c("chrX"), size = c(156040895L)), .Names = c("chromosome",
                                                                                        "size"), class = "data.frame", row.names = c(NA, -1L))
# chrom_sizes

# hg38 centromere locations
centromeres <- structure(list(chromosome = c("chrX"), start = c(58100000L), end = c(63800000L)),
                         .Names = c("chromosome", "start", "end"), class = "data.frame", row.names = c(NA, -1L
                         ))
# centromeres

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

sample_list=NULL
sample_list_bi=NULL
sample_list_mono=NULL
percent_biallelism=NULL
alleles=NULL
alts_list=NULL
strand_list=NULL
bq_list=NULL
mq_list=NULL
chrom_plot=NULL

for (sample in rn(sample_annot)){
  
  sample_list[[sample]]=open_file(sample)
  
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
  
  sample_list_bi[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='bi',]
  
  sample_list_bi[[sample]][["chrom"]] <- factor(x = sample_list_bi[[sample]][["chrom"]],
                                                levels = chrom_order)
  
  sample_list_mono[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='mono',]
  
  sample_list_mono[[sample]][["chrom"]] <- factor(x = sample_list_mono[[sample]][["chrom"]],
                                                  levels = chrom_order)
  
  df_plot_mono=sample_list_mono[[sample]]
  
  df_plot_bi=sample_list_bi[[sample]]
  
  # compute % biallelism in samples
  percent_biallelism[[sample]]=nrow(sample_list_bi[[sample]])/(nrow(sample_list_bi[[sample]])+nrow(sample_list_mono[[sample]]))*100
}

snp_gene_stat=data.frame(row.names=names(sample_list),
           snp=unlist(lapply(sample_list,function(x){nrow(x)})),
           gene=unlist(lapply(sample_list,function(x){
             x=x
             x=unlist(apply(x,1,function(x){
               sub(";.*","",sub(".*gene_id ","",gtf_chrX_gene[data.table::between(as.numeric(x["pos"]),as.numeric(gtf_chrX_gene$start),as.numeric(gtf_chrX_gene$end)),"attribute"]))}))
             length(unique(x))
           })),
           sample_category=sample_annot$sample_category)

# write.tsv(snp_gene_stat,file="xact/chrX_snp_gene_stat_biallelism.tsv")

df_plot_list=NULL

for(sample in names(sample_list)){
  
  df_plot_list[[sample]]=rbind(sample_list_bi[[sample]],sample_list_mono[[sample]])
  df_plot_list[[sample]]$sample_category=rep(sample_annot[sample,"sample_category"],nrow(df_plot_list[[sample]]))
  df_plot_list[[sample]]$group=rep(sample_annot[sample,"group"],nrow(df_plot_list[[sample]]))
  df_plot_list[[sample]]$cell_type_condition=rep(sample_annot[sample,"cell_type_condition"],nrow(df_plot_list[[sample]]))
}

df_plot=do.call("rbind",df_plot_list)

df_plot$sample_category=factor(df_plot$sample_category,levels = c("primed_non_eroded_WT_vallot",
                                                                  "primed_eroded_WT_vallot",
                                                                  "primed_WT",
                                                                  "primed_KO_XACT",
                                                                  "naive_cult_WT",
                                                                  "naive_cult_KO_XACT",
                                                                  "naive_WT",
                                                                  "naive_KO",
                                                                  "naive_C",
                                                                  "naive_T",
                                                                  "naive_C_gro",
                                                                  "naive_T_gro"))

chrom_sizes_df=data.frame(chromosome=rep(chrom_sizes$chromosome,nlevels(df_plot$sample_category)),
                          size=rep(chrom_sizes$size,nlevels(df_plot$sample_category)),
                          sample_category=levels(df_plot$sample_category))

df_plot %>%
  group_by(sample_category) %>%
  ggplot()+
  geom_rect(data=chrom_sizes_df,aes(xmin = as.numeric(chromosome) - 0.1,
                                    xmax = as.numeric(chromosome) + 0.1,
                                    ymax = size, ymin = 0),
            colour="black", fill = "white") +
  coord_flip() +
  theme(axis.text.x = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
  # add centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.1,
                                    xmax = as.numeric(chromosome) + 0.1,
                                    ymax = end, ymin = start)) +
  geom_rect(aes(xmin = as.numeric(chrom) - 0.097,
                xmax = as.numeric(chrom) + 0.097,
                ymax = pos+150000, ymin = pos-150000, fill=allelism)) +
  facet_wrap(~factor(sample_category),ncol = 2) +
  ylab("region (bp)") +
  ggtitle("Allelic status SNP location on chromosomes")
# # base rectangles for the chroms

# ggsave("claire_2023_xact_snp_location_allelic_status_chrX.pdf",width = 7,height = 7)

percent_biallelism=unlist(percent_biallelism)
# sample_annot_naive_primed=sample_annot[sample_annot$cell_type_condition%in%c("primed_control","naive_control") & !(grepl("gro_seq",sample_annot$dataset)),]

bi_percent_df=data.frame(row.names=names(percent_biallelism),
                         perbi=as.numeric(percent_biallelism),
                         cell_type_condition=sample_annot$cell_type_condition,
                         dataset=sample_annot$dataset,clone=sample_annot$clone_or_replicate
                         
)

bi_percent_df %>%
  filter(dataset !="rna_seq_dvallot") %>%
  filter(cell_type_condition %in% c("naive_control","naive_cult_control","primed_control") & dataset!="gro_seq_d1271") %>%
  ggplot() +
  geom_boxplot(aes(x=cell_type_condition, y=perbi, fill=cell_type_condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  geom_point(aes(x=cell_type_condition, y=perbi, fill=cell_type_condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  scale_x_discrete(expand = c(0, 0.5)) +
  ylab("% biallelism\nall heterozygous SNPs detected") +
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
  # stat_boxplot(geom='errorbar', linetype=1, width=0.8, size=0.25)+
  scale_y_continuous(limits=c(0,100))+
  ggtitle("chrX %biallelism in stem cells") 
# -> p
# 
# res.stat <- bi_percent_df %>%
#   filter(dataset !="rna_seq_dvallot") %>%
#   filter(cell_type_condition %in% c("naive_control","naive_cult_control","primed_control") & dataset!="gro_seq_d1271") %>%
#   t_test(data = ., perbi ~ cell_type_condition,paired = F) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj") %>%
#   add_significance("p") %>%
#   add_xy_position(fun="max",scales = "free",formula = perbi ~ cell_type_condition,dodge = 0.9)
# 
# # pdf("xist_and_markers_htseq_gene_expression.pdf",width = 9,height = 6)
# 
# print(p +
#         stat_pvalue_manual(
#           res.stat, label = "p.signif", y.position = res.stat$y.position+5
#         )+
#         stat_pvalue_manual(
#           res.stat, label = "p", y.position = res.stat$y.position+8
#         )
# 
# )
# 
# #ggsave("xact/naive_cult_wt_vs_xact_ko_percent_biallelism_chrX.pdf",width = 8, height = 6)

# chr7

open_file <- function(sample){
  snp <- read.table(paste0("vcf/", sample, "_samtools_mpileup_informative_from_wgs_h9_chr7.pileup"), header = FALSE,
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


sample_annot=fastRead("sample_annot_all.txt", header=T, sep=",",as.matrix = F)
sample_annot=sample_annot[sample_annot$dataset!="rna_seq_d1507",]
sample_annot$sample_id=rn(sample_annot)
sample_annot$cell_type_condition=paste(sep="_",sample_annot$cell_type,sample_annot$condition)

sample_annot$sample_category=paste(sep="_",sample_annot$cell_type_condition,
                                   sample_annot$dox_ko,
                                   sample_annot$dataset)

sample_annot$sample_category=str_replace(sample_annot$sample_category,"_control","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_rna","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_treatment","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_seq","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d1269","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_bis","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d165","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_d1271","")
sample_annot$sample_category=str_replace(sample_annot$sample_category,"_dvallot","_vallot")

sample_annot$sample_category=factor(sample_annot$sample_category,
                                    levels = c("primed_WT",
                                               "primed_non_eroded_WT_vallot",
                                               "primed_eroded_WT_vallot",
                                               "primed_KO_XACT",
                                               "naive_cult_WT",
                                               "naive_cult_KO_XACT",
                                               "naive_WT",
                                               "naive_KO",
                                               "naive_C",
                                               "naive_T",
                                               "naive_C_gro",
                                               "naive_T_gro"))      

sample_annot$group=as.vector(sample_annot$sample_category)
sample_annot$group[sample_annot$group%in%c("primed_WT","primed_KO_XACT")]<-"group1"
sample_annot$group[sample_annot$group%in%c("naive_cult_WT","naive_cult_KO_XACT")]<-"group2"
sample_annot$group[sample_annot$group%in%c("primed_non_eroded_WT_vallot","primed_eroded_WT_vallot")]<-"group3"
sample_annot$group[sample_annot$group%in%c("naive_WT","naive_KO")]<-"group4"
sample_annot$group[sample_annot$group%in%c("naive_C","naive_T")]<-"group5"
sample_annot$group[sample_annot$group%in%c("naive_C_gro","naive_T_gro")]<-"group6"

sample_annot$cell_type_condition=paste(sep="_",sample_annot$cell_type,sample_annot$condition)

ascii_qual=read.tsv("ascii_qual.tsv",comment.char = "")
ascii_qual$symbol=rn(ascii_qual)

informative_pos_vcf_chr7=read.table('H9_WGS_hg38_filtered_chr7.vcf')
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

# # write.tsv(data.frame(informative_pos_vcf_chr7$chrom,
#                      informative_pos_vcf_chr7$pos,
#                      informative_pos_vcf_chr7$pos),file="chr7_heterozygous_snp_h9.bed",row.names = F, col.names = F)

# chromosome_plots_biallelism

# remember to add chr7 as a control

chrom_sizes <- structure(list(chromosome = c("chr7"), size = c(156040895L)), .Names = c("chromosome",
                                                                                        "size"), class = "data.frame", row.names = c(NA, -1L))
# chrom_sizes

# hg38 centromere locations
centromeres <- structure(list(chromosome = c("chr7"), start = c(58100000L), end = c(63800000L)),
                         .Names = c("chromosome", "start", "end"), class = "data.frame", row.names = c(NA, -1L
                         ))
# centromeres

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

sample_list=NULL
sample_list_bi=NULL
sample_list_mono=NULL
percent_biallelism=NULL
alleles=NULL
alts_list=NULL
strand_list=NULL
bq_list=NULL
mq_list=NULL
chrom_plot=NULL

for (sample in rn(sample_annot)){
  
  sample_list[[sample]]=open_file(sample)
  
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
  
  sample_list_bi[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='bi',]
  
  sample_list_bi[[sample]][["chrom"]] <- factor(x = sample_list_bi[[sample]][["chrom"]],
                                                levels = chrom_order)
  
  sample_list_mono[[sample]]=sample_list[[sample]][sample_list[[sample]]$allelism=='mono',]
  
  sample_list_mono[[sample]][["chrom"]] <- factor(x = sample_list_mono[[sample]][["chrom"]],
                                                  levels = chrom_order)
  
  df_plot_mono=sample_list_mono[[sample]]
  
  df_plot_bi=sample_list_bi[[sample]]
  
  # compute % biallelism in samples
  percent_biallelism[[sample]]=nrow(sample_list_bi[[sample]])/(nrow(sample_list_bi[[sample]])+nrow(sample_list_mono[[sample]]))*100
}

df_plot_list=NULL

for(sample in names(sample_list)){
  
  df_plot_list[[sample]]=rbind(sample_list_bi[[sample]],sample_list_mono[[sample]])
  df_plot_list[[sample]]$sample_category=rep(sample_annot[sample,"sample_category"],nrow(df_plot_list[[sample]]))
  df_plot_list[[sample]]$group=rep(sample_annot[sample,"group"],nrow(df_plot_list[[sample]]))
  df_plot_list[[sample]]$cell_type_condition=rep(sample_annot[sample,"cell_type_condition"],nrow(df_plot_list[[sample]]))
}

df_plot=do.call("rbind",df_plot_list)

df_plot$sample_category=factor(df_plot$sample_category,levels = c("primed_non_eroded_WT_vallot",
                                                                  "primed_eroded_WT_vallot",
                                                                  "primed_WT",
                                                                  "primed_KO_XACT",
                                                                  "naive_cult_WT",
                                                                  "naive_cult_KO_XACT",
                                                                  "naive_WT",
                                                                  "naive_KO",
                                                                  "naive_C",
                                                                  "naive_T",
                                                                  "naive_C_gro",
                                                                  "naive_T_gro"))

chrom_sizes_df=data.frame(chromosome=rep(chrom_sizes$chromosome,nlevels(df_plot$sample_category)),
                          size=rep(chrom_sizes$size,nlevels(df_plot$sample_category)),
                          sample_category=levels(df_plot$sample_category))

df_plot %>%
  group_by(sample_category) %>%
  ggplot()+
  geom_rect(data=chrom_sizes_df,aes(xmin = as.numeric(chromosome) - 0.1,
                                    xmax = as.numeric(chromosome) + 0.1,
                                    ymax = size, ymin = 0),
            colour="black", fill = "white") +
  coord_flip() +
  theme(axis.text.x = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
  # add centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.1,
                                    xmax = as.numeric(chromosome) + 0.1,
                                    ymax = end, ymin = start)) +
  geom_rect(aes(xmin = as.numeric(chrom) - 0.097,
                xmax = as.numeric(chrom) + 0.097,
                ymax = pos+150000, ymin = pos-150000, fill=allelism)) +
  facet_wrap(~factor(sample_category),ncol = 2) +
  ylab("region (bp)") +
  ggtitle("Allelic status SNP location on chromosomes")
# # base rectangles for the chroms

# ggsave("claire_2023_xact_snp_location_allelic_status_chr7.pdf",width = 7,height = 7)

percent_biallelism=unlist(percent_biallelism)
# sample_annot_naive_primed=sample_annot[sample_annot$cell_type_condition%in%c("primed_control","naive_control") & !(grepl("gro_seq",sample_annot$dataset)),]

bi_percent_df=data.frame(row.names=names(percent_biallelism),
                         perbi=as.numeric(percent_biallelism),
                         cell_type_condition=sample_annot$cell_type_condition,dataset=sample_annot$dataset,clone=sample_annot$clone_or_replicate
                         
)

bi_percent_df %>%
  filter(dataset !="rna_seq_dvallot") %>%
  filter(cell_type_condition %in% c("naive_control","naive_cult_control","primed_control") & dataset!="gro_seq_d1271") %>%
  ggplot() +
  geom_boxplot(aes(x=cell_type_condition, y=perbi, fill=cell_type_condition),lwd=0.15,position=position_dodge(width = 0.8), outlier.stroke = 0, outlier.size = 0)+
  geom_point(aes(x=cell_type_condition, y=perbi, fill=cell_type_condition),size=0.25,position=position_jitterdodge(jitter.width=0,dodge.width=0.8),alpha=0.8)+
  theme_bw()+
  theme(strip.text.x = element_text( face = "bold.italic", size = 10))+
  scale_x_discrete(expand = c(0, 0.5)) +
  ylab("% biallelism\nall heterozygous SNPs detected") +
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
  # stat_boxplot(geom='errorbar', linetype=1, width=0.8, size=0.25)+
  scale_y_continuous(limits=c(0,120))+
  ggtitle("chr7 %biallelism in stem cells") 

# -> p
# 
# res.stat <- bi_percent_df %>%
#   filter(dataset !="rna_seq_dvallot") %>%
#   filter(cell_type_condition %in% c("naive_control","naive_cult_control","primed_control") & dataset!="gro_seq_d1271") %>%
#   t_test(data = ., perbi ~ cell_type_condition,paired = F) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj") %>%
#   add_significance("p") %>%
#   add_xy_position(fun="max",scales = "free",formula = perbi ~ cell_type_condition,dodge = 0.9)
# 
# # pdf("xist_and_markers_htseq_gene_expression.pdf",width = 9,height = 6)
# 
# print(p +
#         stat_pvalue_manual(
#           res.stat, label = "p.signif", y.position = res.stat$y.position+1
#         )+
#         stat_pvalue_manual(
#           res.stat, label = "p", y.position = res.stat$y.position+10
#         )
# 
# )

#ggsave("xact/naive_cult_wt_vs_xact_ko_percent_biallelism_chr7.pdf",width = 8, height = 6)
