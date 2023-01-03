#!/bin/bash
#SBATCH --job-name=joint_unique
#SBATCH --partition=sunlab
#SBATCH --nodes=1

hormone="SHBG"
interaction=(BMI T2D CAD)
sex=(Female Male)

for x in "${interaction[@]}";
do for y in "${sex[@]}";
do dir="/projects/sunlab/Students_Work/Amonae_work/GEM_$hormone/$x/$y";
cd $dir;

gunzip *FUMA.txt.gz

ml R/4.0.3
Rscript --vanilla joint_unique.R

gzip *FUMA.txt 
done;
done
