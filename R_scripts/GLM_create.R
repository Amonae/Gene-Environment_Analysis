#### This creates lm/ glm models for binary and continuous outcome variables
##### Model is printed in the logfile, but not saved to a file

rm(list=ls())
suppressPackageStartupMessages(library(data.table,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(crayon,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(dplyr,lib.loc="/projects/sunlab/R.lib"))

# Loading data

pheno_file = fread(Sys.glob("*.pheno"))

# removing FID column
if("FID" %in% colnames(pheno_file))
{
  print("Removing FID column")
  pheno_file = pheno_file[,-c("FID")]
}


# Changing names for standardization

hormone_name = names(pheno_file)[13]
outcome_name = names(pheno_file)[14]
names(pheno_file)[13:15] = c("hormone", "outcome", "PRS_outcome")

### Keeping track of hormone and outcome for logfile
print(getwd())
print(paste( "Hormone:", hormone_name)) 
print(paste( "Outcome:", outcome_name))



#### Creating glm model

binary_check = apply(pheno_file,2,function(x) { all(x %in% 0:1) }) # If true, outcome is binary

if(binary_check[14] == FALSE){
	model = lm(formula = outcome ~ Age_enrolled + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
	PC9 + PC10 + PRS_outcome * hormone, data = pheno_file)} else {model = glm(formula = outcome ~ Age_enrolled + 
	PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PRS_outcome * hormone, data = pheno_file, 
	family = binomial())}

cat("Regression model for", paste0("PRS_", outcome_name, "x", hormone_name), "\n") 
print(summary(model))

