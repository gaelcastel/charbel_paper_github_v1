#!/bin/bash

#SBATCH -o rap_h9_multi_bw_summary_downsampled.%N.%j.out
#SBATCH -e rap_h9_multi_bw_summary_downsampled.%N.%j.err
#SBATCH --partition ipop-up
#SBATCH --cpus-per-task 8
#SBATCH --ntasks-per-node 1
#SBATCH --mem 48GB

module load deeptools/3.5.0

# srun multiBigwigSummary BED-file --BED rmsk_x.bed --numberOfProcessors max --verbose -b /mnt/c/Users/gael/charbel_2022/rap_h9/bw/*.hg38.BPM.all.bw -o rap_h9_bw_summary_transposable_elements.npz --outRawCounts rap_h9_bw_summary_transposable_elements.tab
# srun multiBigwigSummary BED-file --BED rmsk_x.bed --numberOfProcessors 8 --verbose -b /shared/projects/xci/homo_sapiens/stem_cells/charbel_2022/rap_h9/bigwig/*_downsampled.bw -o rap_h9_bw_summary_transposable_elements_downsampled.npz --outRawCounts rap_h9_bw_summary_transposable_elements_downsampled.tab
srun multiBigwigSummary bins -b /shared/projects/xci/homo_sapiens/stem_cells/charbel_2022/rap_h9/bigwig/*_downsampled.bw -o rap_h9_no_igg_norm_mulitbigwigsummary_downsampled.npz --outRawCounts rap_h9_no_igg_norm_scores_per_bin_downsampled.tab

exit
