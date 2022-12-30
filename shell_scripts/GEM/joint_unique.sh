#!/bin/bash
#SBATCH --job-name=joint_unique
#SBATCH --partition=sunlab
#SBATCH --nodes=1

hormone="SHBG"
interaction="BMI"
sex="Male"

dir="/projects/sunlab/Students_Work/Amonae_work/GEM_$hormone/$interaction/$sex"
cd $dir

gunzip *FUMA.txt.gz

ml R/4.0.3
joint_unique.R

gzip *FUMA.txt 