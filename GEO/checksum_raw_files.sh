#!/bin/bash

#SBATCH -o checksum.%N.%j.out
#SBATCH -e checksum.%N.%j.err
#SBATCH --partition ipop-up
#SBATCH --cpus-per-task 4
#SBATCH --ntasks-per-node 1
#SBATCH --mem 16GB

