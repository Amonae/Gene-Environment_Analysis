# Gene-Environment_Analysis
Using the GEM tool to conduct GxE analysis

**Gene-Envorinment Million (GEM) was created by Westerman et al. 2021**  
[GEM: scalable and flexible gene–environment interaction analysis in millions of samples](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8545347/)

This is a repository to store my shell and R scripts:

### R-  
1. covar_pheno_create.R- creates the covar_pheno.txt file that includes covariates and the exposure (phenotype). Outcome polygenic risk scores (PRS) and phenotypes should be added individually as each requires a different protocol. Once the outcome is added to covar_pheno.txt, it should be saved as *.pheno <br/><br/> 
2. GLM_create.R- creates a glm or lm model of **outcome~ age + PC1-10 + PRS*exposure** from the *.pheno file depending on whether the outcome is binary (glm) or continuous (lm). Once run on the cluster, the models are printed in the log file. <br/><br/> 
3. LM_plots.R- Used to created bar plots and linear regression plots for continuous outcomes to better view the interaction between PRS and exposure. This should be used for continuous outcomes. <br/><br/>
4. oddsratio_2/4prs.R- These are used to calculate odds ratios in either 2 or 4 PRS groups and 2 or 4 exposure groups for binary outcomes. It also produces forest and bar plots for visualization. <br/><br/>
5. GEMtoFUMA.R- This concatenates the chromosome output files generated after GEM analysis, removes SNPs with MAF <= 0.01, and saves files with all effects as well as joint, interaction, and marginal effects individually. The individual effects files are compressed and submitted to FUMA. <br/><br/>
6. interaction_loci.R - This uses the interaction effects GenomicRiskLoci.txt output file from FUMA and compares the p values of SNPs included to their counterparts in marginal and joint effects files.  <br/><br/>
7. joint_unique.R - This uses the joint and marginal effects GenomicRiskLoci.txt output files from FUMA and identifies risk loci that are identified using joint effects but do not overlap marginal effects genomic risk loci. Loci here are defined as 250kb regions on either side of SNPs included in the GenomicRiskLoci files. 
<br/><br/>

### Shell- 

1. plink_gwasl.sh- This is used to run gwas via plink. In this analysis, I compared its output to GEM marginal effects. It was confirmed that marginal effects are the same as regular GWAS for outcome~ covariates+ genotype. <br/><br/>
2. BMI_testo_GEM_chr1_22.sh- This is used to run GEM analysis on chromosomes 1-22. <br/><br/>
3. unzip_FUMA.sh- This just unzips the FUMA results files. <br/><br/>
4. joint_unique.sh- This is used in conjunction with the joint_unique.R file to unzip the FUMA.txt.gz files, run R, then rezip the files. <br/><br/>
5. generic.sh- This is a versatile file for running R scripts.
