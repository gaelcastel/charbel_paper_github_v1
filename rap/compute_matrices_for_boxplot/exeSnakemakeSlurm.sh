#!/bin/bash

module load snakemake/5.7.4
module load deeptools/3.5.0
module load wiggletools/1.2.11

snakemake --unlock
snakemake -rp -j 200 --resources load=100 --cluster "sbatch -p {cluster.partition} --cpus-per-task {cluster.cpu} --mem {cluster.ram}" --cluster-config cluster_config.json --latency-wait 1800 --max-jobs-per-second 1 --configfile config.json

