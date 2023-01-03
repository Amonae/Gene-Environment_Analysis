#!/bin/bash
#SBATCH --job-name=pre_FUMA
#SBATCH --partition=sunlab
#SBATCH --nodes=1

hormone="SHBG"
interaction=(BMI T2D CAD)
sex=(Female Male)

for x in "${interaction[@]}";
do for y in "${sex[@]}";
do dir="/projects/sunlab/Students_Work/Amonae_work/GEM_$hormone/$x/$y";
cd $dir;

ml R/4.0.3
Rscript --vanilla pre_FUMA.R

gzip \*.txt 
done;
done
