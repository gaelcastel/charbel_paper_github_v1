#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

.libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))

library(ggplot2)
library(stringr)
library(dplyr)

for (i in seq(length(args))){

  setwd(paste0("/shared/projects/xci/",args[i],"/stem_cells/charbel_2022/cutrun_2023/"))

  if(!dir.exists("Figures")) dir.create("Figures")

  meta <- read.table(paste0("metadata_",args[i],".csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

  list_files <- list.files("stats", pattern = ".coli")

  read_number <- lapply(list_files, function(x){
    fil <- read.table(paste0("stats/", x), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    rea <- fil$V3[1]
    print(rea)
    name <- unlist(str_split(x, "\\."))[1]
    return(c(name, rea))
  })


  read_coli <- data.frame(matrix(unlist(read_number), ncol = 2, byrow = 1))

  colnames(read_coli) <- c("sample", "read")

  read_coli$Sample_ID <- sapply(read_coli$sample, function(x){
    paste0("C",unlist(str_split(x, "C"))[2])
  })

  read_coli <- merge(read_coli, meta, by = "Sample_ID")

  read_coli$read <- as.numeric(read_coli$read)
  read_coli$mark <- sapply(read_coli$Sample_Name, function(x){
    unlist(str_split(x, "_"))[2]
  })

  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ggplot(read_coli, aes(x = Sample_Name, y = read, fill = mark)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("reads number") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size =15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=cbPalette)
  ggsave("Figures/Number_reads_coli.png", width = 15, height = 15)
  ggsave("Figures/Number_reads_coli.pdf", width = 15, height = 15)

  read_coli$scale_factor <- 10000/read_coli$read                  

  ggplot(read_coli, aes(x = Sample_Name, y = scale_factor, fill = mark)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("scale factor") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size =15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=cbPalette)
  ggsave("Figures/Scale_factor_coli.png", width = 15, height = 15)
  ggsave("Figures/Scale_factor_coli.pdf", width = 15, height = 15)

  marks <- unique(read_coli$mark)



  #Calculate mean for each replicate
  means <- read_coli %>%
    group_by(mark) %>%
    summarise(mean = mean(read))

  read_coli <- merge(read_coli, means, by = "mark")

  read_coli$scale_factor_by_mark <- read_coli$mean / read_coli$read

  ggplot(read_coli, aes(x = Sample_Name, y = scale_factor_by_mark, fill = mark)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("scale factor") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size =15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=cbPalette)
  ggsave("Figures/Scale_factor_by_mark_coli.png", width = 15, height = 15)
  ggsave("Figures/Scale_factor_by_mark_coli.pdf", width = 15, height = 15)

  ##################################################
  #   For the yeasts
  #################################################

  list_files <- list.files("stats", pattern = ".yeast")

  read_number <- lapply(list_files, function(x){
    fil <- read.table(paste0("stats/", x), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    rea <- sum(fil$V3)
    # rea <- fil$V3[fil$V1 == "chrXIII"]
    print(rea)
    name <- unlist(str_split(x, "\\."))[1]
    return(c(name, rea))
  })


  read_yeast <- data.frame(matrix(unlist(read_number), ncol = 2, byrow = 1))

  colnames(read_yeast) <- c("sample", "read")

  read_yeast$Sample_ID <- sapply(read_yeast$sample, function(x){
    paste0("C",unlist(str_split(x, "C"))[2])
  })

  read_yeast <- merge(read_yeast, meta, by = "Sample_ID")

  #read_yeast$Sample_Name <- sapply(read_yeast$Sample_Name, function(x){
  #  unlist(str_split(x, "CUTandRUN_"))[2]
  #})

  read_yeast$read <- as.numeric(read_yeast$read)
  read_yeast$mark <- sapply(read_yeast$Sample_Name, function(x){
    unlist(str_split(x, "_"))[2]
  })

  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ggplot(read_yeast, aes(x = Sample_Name, y = read, fill = mark)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("reads number") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size =15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=cbPalette)
  ggsave("Figures/Number_reads_yeast_all_chromosomes.png", width = 15, height = 15)
  ggsave("Figures/Number_reads_yeast_all_chromosomes.pdf", width = 15, height = 15)

  read_yeast$scale_factor <- 10000/read_yeast$read                  

  ggplot(read_yeast, aes(x = Sample_Name, y = scale_factor, fill = mark)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("scale factor") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size =15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=cbPalette)
  ggsave("Figures/Scale_factor_yeast_all_chromosomes.png", width = 15, height = 15)
  ggsave("Figures/Scale_factor_yeast_all_chromosomes.pdf", width = 15, height = 15)

  marks <- unique(read_yeast$mark)


  #Calculate mean for each replicate
  means <- read_yeast %>%
    group_by(mark) %>%
    summarise(mean = mean(read))

  read_yeast <- merge(read_yeast, means, by = "mark")

  read_yeast$scale_factor_by_mark <- read_yeast$mean / read_yeast$read

  ggplot(read_yeast, aes(x = Sample_Name, y = scale_factor_by_mark, fill = mark)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("scale factor") +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size =15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=cbPalette)
  ggsave("Figures/Scale_factor_by_mark_yeast_all_chromosomes.png", width = 15, height = 15)
  ggsave("Figures/Scale_factor_by_mark_yeast_all_chromosomes.pdf", width = 15, height = 15)


  corr_yeast_coli <- data.frame(Sample_Name = read_yeast$Sample_Name,
                                yeast = read_yeast$read,
                                coli = read_coli$read,
                                mark = read_yeast$mark)

  require(dplyr)
  require(broom)

  corr <- corr_yeast_coli  %$%
    cor.test(coli, yeast) %>%
    tidy %>%
    mutate_if(is.numeric, round, 4)

  corr

  text <- paste0('r = ', corr$estimate, ', p-value = ', corr$p.value)

  text

  ggplot(corr_yeast_coli, aes(x = coli, y = yeast)) +
  geom_point(aes(color = mark)) +
  #xlim(0, max(corr_yeast_coli$coli) + 100) +
  #ylim(0, max(corr_yeast_coli$yeast) + 100) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate('text',  x = 2000, y = 35, label=text) +
  theme_bw() +
  ylab("Yeast (reads number)") +
  xlab("E. coli (reads number)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size =15))
  ggsave("Figures/correlation_yeast_coli_all_chromosomes.png", width = 15, height = 15)
  ggsave("Figures/correlation_yeast_coli_all_chromosomes.pdf", width = 15, height = 15)

  corrscale_yeast_coli <- data.frame(Sample_Name = read_yeast$Sample_Name,
                                    Sample_ID = paste0("D1249", read_yeast$Sample_ID),
                                yeast = read_yeast$scale_factor,
                                coli = read_coli$scale_factor,
                                mark = read_yeast$mark)

  require(dplyr)
  require(broom)

  corr <- corrscale_yeast_coli  %$%
    cor.test(coli, yeast) %>%
    tidy %>%
    mutate_if(is.numeric, round, 4)

  corr

  text <- paste0('r = ', corr$estimate, ', p-value = ', corr$p.value)

  text

  ggplot(corrscale_yeast_coli, aes(x = coli, y = yeast)) +
  geom_point(aes(color = mark)) +
  #xlim(0, max(corr_yeast_coli$coli) + 100) +
  #ylim(0, max(corr_yeast_coli$yeast) + 100) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate('text',  x = 20, y = 2, label=text) +
  theme_bw() +
  ylab("Yeast (scale factor)") +
  xlab("E. coli (scale factor)") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size =20),
        axis.title.x = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size =15))
  ggsave("Figures/correlation_scale_factor_yeast_coli_all_chromosomes.png", width = 15, height = 15)
  ggsave("Figures/correlation_scale_factor_yeast_coli_all_chromosomes.pdf", width = 15, height = 15)


  write.table(corrscale_yeast_coli, file = paste0("stats/",args[i],"_scale_factors_yeast_ecoli_all_chromosomes.txt"), col.names = TRUE, sep = "\t",
              row.names = FALSE, quote = FALSE)

  #For each histone mark
  each_histone_scale_yeast_coli <- data.frame(Sample_Name = read_yeast$Sample_Name,
                                    Sample_ID = paste0("D1249", read_yeast$Sample_ID),
                                    yeast = read_yeast$scale_factor_by_mark,
                                    coli = read_coli$scale_factor_by_mark,
                                    mark = read_yeast$mark)

  write.table(each_histone_scale_yeast_coli, file = paste0("stats/",args[i],"_scale_factors_by_mark_yeast_ecoli_all_chromosomes.txt"), col.names = TRUE, sep = "\t",
              row.names = FALSE, quote = FALSE)
}
