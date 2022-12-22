#!/bin/bash
#SBATCH --job-name=GEM
#SBATCH --partition=sunlab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

ml intel

imputedir="/projects/sunlab/UKB/Data_raw"
dir="/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female"
pheno=$dir"/BMI_testo_f.pheno"
outdir="/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/out"

for ((i=1; i<=22; i++)); do

input=$imputedir"/ukb_imp_chr"$i"_v3"

/projects/sunlab/bin/GEM-1.4.2/src/GEM \
--pfile $input \
--pheno-file $pheno \
--sampleid-name FID \
--pheno-name BMI \
--covar-names Age_enrolled PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
--exposure-names Testosterone \
--output-style meta \
--robust 1 \
--missing-value NA \
--threads 16 \
--delim \0 \
--out $outdir"/testo_bmi_f_chr"$i".txt"

done
