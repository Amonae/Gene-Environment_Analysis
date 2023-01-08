## Calculate Odds Ratios for different groups

# load all R libraries and remove objects from environment

rm(list=ls())
library(data.table,lib.loc="/projects/sunlab/R.lib")
library(crayon,lib.loc="/projects/sunlab/R.lib")
library(dplyr,lib.loc="/projects/sunlab/R.lib")
library(cli,lib.loc="/projects/sunlab/R.lib")
library(epitools, lib.loc = "/projects/sunlab/Students_Work/Amonae_work/R.lib") # Note different lib loc


# Loading data

pheno_file = fread(Sys.glob("*.pheno"))

# Changing names for standardization

hormone_name = names(pheno_file)[13]
outcome_name = names(pheno_file)[14]
names(pheno_file)[13:15] = c("hormone", "outcome", "PRS_outcome")

# Creating dataframes by PRS group (PRS_groups are based on quartiles)
hormone = summary(pheno_file$hormone)  
PRS = summary(pheno_file$PRS_outcome)      

pheno_file$PRS_group = factor(ifelse(pheno_file$PRS_outcome<= PRS[2],"very low", ifelse(pheno_file$PRS_outcome> PRS[2] & pheno_file$PRS_outcome<= PRS[3],"low", ifelse(pheno_file$PRS_outcome> PRS[3] & pheno_file$PRS_outcome<= PRS[5],"high", "very high"))), levels = c("very low", "low", "high", "very high"), ordered = T) # cutoffs = 1st Qu, median, 3rd Qu


pheno_file$hormone_group = factor(ifelse(pheno_file$hormone <= hormone[2],"very low", ifelse(pheno_file$hormone >hormone[2] & pheno_file$hormone <= hormone[3],"low", ifelse(pheno_file$hormone > hormone[3] & pheno_file$hormone <=hormone[5],"high", "very high"))), levels = c("very low", "low", "high", "very high"), ordered = T) # cutoffs = 1st Qu, median, 3rd Qu 

low2 = pheno_file[pheno_file$PRS_group == "very low" | pheno_file$PRS_group == "low",] # for comparing 2 PRS levels and 2 hormone levels
low2$hormone_group = factor(ifelse(low2$hormone_group == "very low" | low2$hormone_group == "low", "low", "high"), levels = c("low", "high"), ordered = T) 

high2 = pheno_file[pheno_file$PRS_group == "high" | pheno_file$PRS_group == "very high",] # for comparing 2 PRS levels
high2$hormone_group = factor(ifelse(high2$hormone_group == "very low" | high2$hormone_group == "low", "low", "high"), levels = c("low", "high"), ordered = T) 


very.low4 = pheno_file[pheno_file$PRS_group == "very low",]
low4 = pheno_file[pheno_file$PRS_group == "low",]
high4 = pheno_file[pheno_file$PRS_group == "high",]
very.high4 = pheno_file[pheno_file$PRS_group == "very high",]

# Creating a list of the dfs
 
dfs =  list( low2, high2, very.low4, low4, high4, very.high4)
names(dfs) = c("low2", "high2", "very.low4", "low4", "high4", "very.high4")

# creating results object
results_df = data.frame(matrix(nrow = length(dfs), ncol = 5))
columns = c("Odds Ratio", "95% C.I. Lower", "95% C.I. Upper", "Fisher pval", "Chi2 pval") 
colnames(results_df) = columns
rownames(results_df) = names(dfs)

print(paste( "Hormone:", hormone_name)) 
print(paste( "Outcome:", outcome_name))

for(i in 1:length(dfs)){
	cont_table =  xtabs(~hormone_group + outcome, data = dfs[[i]])
	result = oddsratio(cont_table, y = NULL,
          method = c("wald"),
          conf.level = 0.95,
          rev = c("neither"),
          correction = FALSE,
          verbose = FALSE)
	print(paste("Results for", names(dfs)[i]))
	print(result)
	results_df[i,1] = result$measure[2,1]
	results_df[i,2] = result$measure[2,2]
   	results_df[i,3] = result$measure[2,3]
    	results_df[i,4] = result$p.value[2,2]
    	results_df[i,5] = result$p.value[2,3]	
	}
	
write.csv(results_df, file = "Odds_ratios.csv")	
