##### This prepares GEM results for FUMA analysis. It combines the chromosomes output, eliminates entries with 
###   MAF<=0.01, and saves files with all data, as well as joint, interaction, and marginal data seperately  

# remove everything in environment
rm(list=ls())

cat("\n \n \n", getwd(), "\n")

# Setting Paths

r.libs = "/projects/sunlab/R.lib"

# Load libraries

suppressPackageStartupMessages(library(data.table,lib.loc= r.libs))
suppressPackageStartupMessages(library(crayon,lib.loc= r.libs))
suppressPackageStartupMessages(library(dplyr,lib.loc= r.libs))
suppressPackageStartupMessages(library(cli,lib.loc= r.libs))

# Function to extract data

multmerge = function(path){
  filenames=list.files(path=path, full.names=TRUE)
  rbindlist(lapply(filenames, fread))
}

## Applying function to extract data

cat("Combining chromosome data \n")
path = "out"
data = multmerge(path)

cat("Done! \nDimensions: \n", dim(data), "Removing SNPS with MAF <= 0.01 \n")
before = nrow(data)

data$MAF = ifelse(data$AF<0.5,data$AF,1-data$AF) # creating a MAF column
data = data[data$MAF>0.01] # Filtering by MAF >0.01

after = nrow(data)
cat("Done! \n", before-after, "SNPS removed \n")

### Viewing data

n_samples = unlist(data[1,6])
cat("There are", n_samples, "samples. YOU WILL NEED TO NOTE THIS FOR FUMA \n")

cat("There are", nrow(data[data$robust_P_Value_Marginal<5e-08]), 
	"SNPS with robust Marginal P values <5e-08 \n")

cat("There are", nrow(data[data$robust_P_Value_Interaction<5e-08]), 
	"SNPS with robust Interaction P values <5e-08 \n")

cat("There are", nrow(data[data$robust_P_Value_Joint<5e-08]), 
	"SNPS with robust Joint P values <5e-08 \n")

#### Saving data

cat("Saving data to", getwd(), "as GEM_results_MAF01.txt \n")

write.table(data, file = "GEM_results_MAF01.txt", row.names = F, quote = F)


#### Creating marginal, joint, and interaction effects dataframes for FUMA analysis

cat("Creating marginal, joint, and interaction effects dataframes for FUMA analysis \n")

marginal = data[,c(1:5, 8, 9, 19)] 
interaction = data[,c(1:5,12,14,20)]
joint = data[,c(1:5,11,13,21)]

colnames(marginal)[c(4:8)] = c("A0","A1","BETA","SE","P")
colnames(interaction)[c(4:8)] = c("A0","A1","BETA","SE","P")
colnames(joint)[c(4:8)] = c("A0","A1","BETA","SE","P")


write.table(marginal, file = "marginal_FUMA.txt", row.names = F, quote = F)
cat("\n Marginal file saved to", paste0(getwd(), "/marginal_FUMA.txt \n"))

write.table(interaction, file = "interaction_FUMA.txt", row.names = F, quote = F)
cat("\n Interaction file saved to", paste0(getwd(), "/interaction_FUMA.txt \n"))

write.table(joint, file = "joint_FUMA.txt", row.names = F, quote = F)
cat("\n Joint file saved to", paste0(getwd(), "/joint_FUMA.txt \n DONE! \n"))
q()