#!/bin/bash

#SBATCH -o big.%N.%j.out
#SBATCH -e big.%N.%j.err
#SBATCH --partition ipop-up
#SBATCH --cpus-per-task 8
#SBATCH --ntasks-per-node 1
#SBATCH --mem 48GB

samp="$1"
scale="$2"
input="$3"
normalizeUsing="$4"
binSize="$5"
smoothLength="$6"
output="$7"
logs="$8"
specie="$9"


module load deeptools/3.5.0

if [ ${specie} = "coli" ]
then
	sf=($(awk -F "\t" -v sam=${samp} '$2 == sam {print $4}' ${scale}))
else
	sf=($(awk -F "\t" -v sam=${samp} '$2 == sam {print $3}' ${scale}))
fi

echo ${sf}

# bamCoverage -b $input --normalizeUsing $normalizeUsing --scaleFactor $sf --binSize $binSize --smoothLength $smoothLength -o $output 2> $logs

bamCoverage -b $input --normalizeUsing $normalizeUsing --scaleFactor $sf --binSize $binSize --smoothLength $smoothLength -o $output 2> $logs