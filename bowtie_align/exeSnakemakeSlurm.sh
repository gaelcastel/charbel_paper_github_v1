#!/bin/bash

module load bowtie2/2.4.4 #load this version on ipop-up, otherwise fails because lack of libtbb.so.2
module load samtools/1.13
module load snakemake/5.7.4
module load deeptools/3.5.0
module load picard/2.23.5
module load r/4.1.1 #load this version on ipop-up, otherwise fails with r/4.0.3
module load htseq/0.13.5
module load fastqc/0.11.9
module load multiqc/1.9
module load trim-galore/0.6.5

snakemake --unlock
snakemake -rp -j 200 --resources load=100 --cluster "sbatch -p {cluster.partition} --cpus-per-task {cluster.cpu} --mem {cluster.ram}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.json

