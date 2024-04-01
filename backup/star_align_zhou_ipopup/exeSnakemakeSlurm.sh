#!/bin/bash

module load snakemake/7.7.0
module load star/2.7.9a
module load multiqc/1.9
module load fastqc/0.11.9

snakemake -rp -j 32 --cluster "sbatch -p {cluster.partition} --cpus-per-task {cluster.cpu} --mem {cluster.ram}" --cluster-config cluster_config.json --latency-wait 3600 --max-jobs-per-second 1 --configfile config.json

exit
