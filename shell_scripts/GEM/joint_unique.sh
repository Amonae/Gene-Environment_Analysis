#!/bin/bash
#SBATCH --job-name=joint_unique
#SBATCH --partition=sunlab
#SBATCH --nodes=1

hormone="Testo"
interaction=(BMI CAD T2D)
sex=(Female Male)

for x in "${interaction[@]}";
do for y in "${sex[@]}";
do dir="GEM_$hormone/$x/$y";
cd $dir;

gunzip *FUMA.txt.gz

ml R/4.0.3
Rscript --vanilla scripts/R/joint_unique.R

gzip *FUMA.txt 

done;
done
