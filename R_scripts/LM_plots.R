#### This creates bar and linear regression plots for continuous outcome variables

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


# Creating dataframes by PRS group (PRS_groups are based on quartiles)
hormone = summary(pheno_file$hormone)  
PRS = summary(pheno_file$PRS_outcome)      

pheno_file$PRS_group = factor(ifelse(pheno_file$PRS_outcome<= PRS[3],"low", "high"), levels = c("low", "high"), ordered = T) # cutoffs = median


pheno_file$hormone_group = factor(ifelse(pheno_file$hormone <= hormone[2],"very low", ifelse(pheno_file$hormone >hormone[2] & pheno_file$hormone <= hormone[3],"low", ifelse(pheno_file$hormone > hormone[3] & pheno_file$hormone <=hormone[5],"high", "very high"))), levels = c("very low", "low", "high", "very high"), ordered = T) # cutoffs = 1st Qu, median, 3rd Qu 


#### Organizing data for plotting

df = group_by(pheno_file,PRS_group, hormone_group)
df = summarise(df, Median = median(outcome), count = n())

cat("Median", outcome_name, "by PRS and", hormone_name, "group \n") 
print(df)


cat("Saving bar and linear regression plots /n")
ggplot(data = subset(df,PRS_group!="NA"),                                      
       aes(x = PRS_group,
           y = Median,
           fill = hormone_group)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.75), width = 0.70)

ggsave(paste0(hormone_name, "x", outcome_name, "barplot.png"))


ggplot(data = subset(pheno_file,PRS_group!="NA"), aes(x = hormone,  y = outcome, color = PRS_group)) +
  geom_smooth(method = "lm")

ggsave(paste0(hormone_name, "x", outcome_name, "LMplot.png"))
