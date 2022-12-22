#!/bin/bash
#SBATCH --job-name=plink_female
#SBATCH --partition=sunlab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6

for ((i=1; i<=22; i++)); do

imputedir="/projects/sunlab/UKB/Data_raw"
outdir="/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/out_plink"

input=$imputedir"/ukb_imp_chr"$i"_v3"

/projects/sunlab/bin/plink2 \
--memory 30000 \
--threads 6 \
--pfile $input \
--maf 0.01 \
--no-psam-pheno \
--pheno pheno_plink.txt \
--pheno-name BMI \
--covar covar.txt \
--linear cols=+a1freq,+machr2,+beta hide-covar \
--covar-variance-standardize \
--out $outdir"/plink_BMI_female_C"$i""

done