#!/bin/bash

#SBATCH -o correlated_bins.%N.%j.out
#SBATCH -e correlated_bins.%N.%j.err
#SBATCH --partition fast
#SBATCH --cpus-per-task 16
#SBATCH --ntasks-per-node 1
#SBATCH --mem 192GB

module load r/4.0.3

srun Rscript projection_r_script_sbatch.R

exit
