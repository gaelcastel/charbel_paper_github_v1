#!/bin/bash

module load bowtie2/2.4.4
module load samtools/1.13
module load snakemake/5.7.4
module load deeptools/3.5.0
module load picard/2.23.5
module load r/4.1.1

# grep -A5000 -m1 -e 'Data' SampleSheet_homo_sapiens.csv |tail -n+2 > metadata_homo_sapiens.csv
# grep -A5000 -m1 -e 'Data' SampleSheet_callitrix_jacchus.csv |tail -n+2 > metadata_callitrix_jacchus.csv
# CAUTION !! ADD column Lane, (1), manually

snakemake --unlock
snakemake -rp -j 100 --resources load=100 --cluster "sbatch -p {cluster.partition} --cpus-per-task {cluster.cpu} --mem {cluster.ram}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.json
