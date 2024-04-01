#!/bin/bash

#SBATCH -o RAP_H9.%N.%j.out
#SBATCH -e RAP_H9.%N.%j.err
#SBATCH --partition long
#SBATCH --cpus-per-task 16
#SBATCH --ntasks-per-node 1
#SBATCH --mem 32GB
#SBATCH --time 2-12:00:00


module load bowtie2/2.4.1
module load samtools/1.13
module load snakemake/5.7.4
module load deeptools/3.5.0
module load picard/2.23.5
module load r/4.0.3

grep -A5000 -m1 -e 'Data' SampleSheet.csv |tail -n+2 > metadata.csv

snakemake --unlock
snakemake -r -j 8 --rerun-incomplete --cluster "sbatch --mem=32G --nodes 1 --ntasks-per-node 1 --cpus-per-task 16"

