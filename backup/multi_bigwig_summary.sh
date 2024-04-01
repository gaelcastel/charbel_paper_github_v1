#!/bin/bash

#SBATCH -o charbel_rap_multi_bw_summary.%N.%j.out
#SBATCH -e charbel_rap_multi_bw_summary.%N.%j.err
#SBATCH --partition fast
#SBATCH --cpus-per-task 8
#SBATCH --ntasks-per-node 1
#SBATCH --mem 64GB

module load deeptools/3.5.0

srun multiBigwigSummary bins -b /shared/projects/dubii2021/gcastel/charbel_2022/RAP_H9/bw_log2ratio/* -o charbel_rap_h9_mulitbigwigsummary.npz --outRawCounts charbel_rap_h9_scores_per_bin.tab

