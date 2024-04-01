#!/bin/bash

#SBATCH -o charbel_rap_multi_bw_summary.%N.%j.out
#SBATCH -e charbel_rap_multi_bw_summary.%N.%j.err
#SBATCH --partition ipop-up
#SBATCH --cpus-per-task 8
#SBATCH --ntasks-per-node 1
#SBATCH --mem 64GB

module load deeptools/3.5.0

srun multiBigwigSummary bins -b /shared/projects/xci/homo_sapiens/stem_cells/charbel_2022/rap_h9/bigwig/* -o rap_h9_no_igg_norm_mulitbigwigsummary.npz --outRawCounts rap_h9_no_igg_norm_scores_per_bin.tab

