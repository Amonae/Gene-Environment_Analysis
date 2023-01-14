#### Script to generate covariate/ pheno file using only white Europeans
#### Age, Sex, and PC1-10 are default covariates. Add more if needed in the appropriate section 
# Note that from UK Biobank, 1 = Male and 0 = Female
#### Only tested using hormone phenotypes so far

rm(list=ls())
suppressPackageStartupMessages(library(data.table,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(crayon,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(dplyr,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(cli,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(withr,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(labeling,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(farver,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(digest,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(ggplot2,lib.loc="/projects/sunlab/R.lib"))
suppressPackageStartupMessages(library(viridisLite,lib.loc="/projects/sunlab/R.lib"))



#### Labeling locations. Edit as needed

cov_loc = "/projects/sunlab/UKB/Data_clean/ukb_core_data_669473.rda" # This file has Age, ethnicity, sex, and PC1-10 data. It will be loaded as 'bd'
pheno_loc = "/projects/sunlab/UKB/Data_clean//ukb29860.csv"  # This has SHBG, Testo, and Estradiol data 

#### Setting covariates and phenotype
covars = c("IID","Age_enrolled","Gen_Ethnic","Ethnic","Sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
pheno_names = c("SHBG", "Estradiol", "Testosterone")
field_ids = c(30830, 30800, 30850)
instances = c(0)
arrays =  c(0)
pheno_ids = paste(paste(field_ids, instances, sep = "-"), arrays, sep = ".")
names(pheno_ids) = pheno_names
cat( "Covariates:\n", covars, "\n")
cat( "Phenotypes:\n", names(pheno_ids), "\n")
cat( "Phenotype Field IDs:\n", pheno_ids, "\n")


#### ********* Do you want the final file to only include complete cases? 
#############  Yes means all covariates and phenotypes must be present for all patients
complete_cases = "yes"  


#### ********* Do you want to create a bar plot that shows phenotype by sex?
####             (note this has only been tested with continuous phenotypes) 
make_plot = "yes"

# Loading covariate file

cat("Loading covariate file\n This takes a while\n")
load(file = cov_loc) # Takes forever to load
cat("File loaded\n")

### Renaming covariates. ADD TO THIS SECTION FOR MORE COVARIATES***************

names(bd)[grepl("eid", names(bd))] = "IID"
names(bd)[names(bd)=="f.21022.0.0"]<-"Age_enrolled"
names(bd)[names(bd)=="f.22006.0.0"]<-"Gen_Ethnic"
names(bd)[names(bd)=="f.21000.0.0"]<-"Ethnic"
names(bd)[names(bd)=="f.31.0.0"]<-"Sex"
names(bd)[names(bd)=="f.22009.0.1"]<-"PC1"
names(bd)[names(bd)=="f.22009.0.2"]<-"PC2"
names(bd)[names(bd)=="f.22009.0.3"]<-"PC3"
names(bd)[names(bd)=="f.22009.0.4"]<-"PC4"
names(bd)[names(bd)=="f.22009.0.5"]<-"PC5"
names(bd)[names(bd)=="f.22009.0.6"]<-"PC6"
names(bd)[names(bd)=="f.22009.0.7"]<-"PC7"
names(bd)[names(bd)=="f.22009.0.8"]<-"PC8"
names(bd)[names(bd)=="f.22009.0.9"]<-"PC9"
names(bd)[names(bd)=="f.22009.0.10"]<-"PC10"


# make df with only covariate columns

covars = c("IID","Age_enrolled","Gen_Ethnic","Ethnic","Sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
cat( "Covariates:\n", covars, "\n")
covar_df = bd[,covars]
rm(bd)



### Summary of covariates
print(summary(covar_df))
cat("\n")
print(table(covar_df$Sex))


###### ******Loading phenotype file

cat("\nLoading phenotype file\n")
pheno_file = fread(pheno_loc, stringsAsFactors=F)
cat("\nPhenotype file loaded\n")

## Renaming columns 
names(pheno_file)[grepl("eid", names(pheno_file))] = "IID"  

for(i in 1:length(field_ids)){
  names(pheno_file)[names(pheno_file)==pheno_ids[i]] = names(pheno_ids)[i]
}

# Making a df with phenotype info
pheno_df = pheno_file[,c("IID", ..pheno_names)]

### Viewing summary stats

print(summary(pheno_df)) 
cat("\n")

# merging pheno and covars dfs by "IID". IDs that are not in both will be dropped

cat("Merging covariate and phenotype dataframes by IID \n Starting dimensions of covariate file: \n",
	dim(covar_df), "\nStarting dimensions of phenotype file: \n", dim(pheno_df), "\n")

covar_pheno = merge(covar_df, pheno_df, by = "IID")

cat("Dimensions of combined dataframe: \n", dim(covar_pheno), "\n")




#Only interested in white British patients (ethnicity 1)

cat("Removing patients not of white British ancestry \n")
before = nrow(covar_pheno)

covar_pheno = subset(covar_pheno, Ethnic %in% c("1","1001","1002","1003"))
after = nrow(covar_pheno)

cat("Removed", before-after, "patients\n")


# Making sure all covariates and phenotypes are present

if(tolower(complete_cases) == "yes"){
	cat("Complete cases set to 'yes' \nRemoving rows with missing values")
	before = nrow(covar_pheno)
	covar_pheno = covar_pheno[complete.cases(covar_pheno),]
	after = nrow(covar_pheno)
	cat("\nRemoved", before-after, "patients\n")}


## Remove Relatedness (entire section from Ellen)

cat("Removing patients with Kinship >=0.0884\n")

kin<-read.table("/projects/sunlab/UKB/Data_raw/ukb34031_rel_s488363.dat",header=T)
# dim(kin) #107156      5
kin<-kin[kin$Kinship>=0.0884,] #2nd degree
# dim(kin) #40230     5
# length(unique(kin$ID1)) #36179
# length(unique(kin$ID2)) #36159
ID<-c(kin$ID1,kin$ID2) #80460
dup<-names(table(ID)[table(ID)>1]) #9776
kin<-kin[(!(kin$ID1%in%dup))&(!(kin$ID2%in%dup)),] #27653

set.seed(1234)
kin$temp<-sample(1:2,dim(kin)[1],replace=T)
kin$exclude<-ifelse(kin$temp==1,kin$ID1,kin$ID2)
# length(unique(kin$exclude)) #27653
exclude<-c(kin$exclude,dup)
# length(exclude) # 37429 = 25653+9776
# length(unique(exclude)) # 37429

before = nrow(covar_pheno)
covar_pheno = subset(covar_pheno, !covar_pheno$IID %in% exclude)
after = nrow(covar_pheno)
cat("\nRemoved", before-after, "patients\n")

cat("Removing ethnicity columns\n")
covar_pheno = covar_pheno[,-c(3,4)] # getting rid of ethnicity columns
cat("Final dataset dimensions:\n", dim(covar_pheno), "\nColumn names:\n", 
	paste(names(covar_pheno), "\n"), "Saving file to", getwd(), "\n")

write.table(covar_pheno, file = "covar_pheno.txt",row.names=F,quote=F)


###### Need to find a way to read the names as variables

if(tolower(make_plot) == "yes"){
	covar_pheno = data.frame(covar_pheno)
	print("Creating phenotype x Sex plots")
	for( name in pheno_names){
		x = covar_pheno[,name]
		Plot = ggplot(subset(covar_pheno, !is.na(x)), aes(x= x, fill = as.factor(Sex)))+ 
		geom_histogram(color = "#e9ecef", alpha= 0.5, position ="identity", bins = 50) + labs(fill = "Sex")
		ggsave(paste0(name,"xSex.png"))
		cat("\n", name, "plot saved")}}

