#!/bin/bash
#SBATCH --job-name=OR2prs
#SBATCH --partition=sunlab
#SBATCH --nodes=1

hormone="Testo"
interaction=(CAD T2D)
sex=(Female Male)

for x in "${interaction[@]}";
do for y in "${sex[@]}";
do dir="/projects/sunlab/Students_Work/Amonae_work/GEM_$hormone/$x/$y";
cd $dir;


ml R/4.0.3
Rscript --vanilla /projects/sunlab/Students_Work/Amonae_work/scripts/R/oddsratio_2prs.R
 

done;
done
