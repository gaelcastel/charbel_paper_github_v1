#!/usr/bin/env Rscript

.libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.1.1/lib/R/library/"))
# .libPaths(c("/shared/projects/xci/r_packages/","/shared/software/conda/envs/r-4.0.3/lib/R/library/"))

args = commandArgs(trailingOnly=TRUE)

human <- args[1]
mouse <- args[2]

output.names <- unlist(strsplit(human, ".hg38"))[1]

if (file.exists(paste0(output.names,"/Filtered_bams"))){
  unlink(paste0(output.names,"/Filtered_bams/"), recursive = TRUE)
}else{
  dir.create(output.names)
}
 
library(XenofilteR, lib="/shared/projects/xci/r_packages/")
library(BiocParallel)

sample.list <- data.frame(human, mouse)

bp.param <- MulticoreParam() #use MulticoreParam so that all workers know about global functions. DO NOT use SnowParam, otherwise biocparallel ERROR: cannot find load functions


XenofilteR(sample.list, destination.folder = output.names, bp.param = bp.param, output.names)
