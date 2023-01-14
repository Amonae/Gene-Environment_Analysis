## Calculate Odds Ratios for different groups

# print working directory
getwd()

# load all R libraries and remove objects from environment

rm(list=ls())
library(data.table)
library(crayon)
library(dplyr)
library(cli)
library(ggplot2)
library(epitools) # Note different lib loc


# Loading data

pheno_file = fread(Sys.glob("*.pheno"))

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

# Creating a list of the PRS_group dfs

dfs =  list( low2, high2, very.low4, low4, high4, very.high4)
names(dfs) = c("low2", "high2", "very.low4", "low4", "high4", "very.high4")

### Creating empty df which will be filled with results

results_df = data.frame(matrix(nrow = 20, ncol = 7))
columns = c("PRS_group", "Hormone_group","Odds_Ratio", "C.I.95_Lower", "C.I.95_Upper", "Fisher_pval", "Chi2_pval") 
colnames(results_df) = columns



### Function to get OR
OR_results = function(x){
  cont_table =  xtabs(~hormone_group + outcome, data = x)
  result = oddsratio(cont_table, y = NULL,
                     method = c("wald"),
                     conf.level = 0.95,
                     rev = c("neither"),
                     correction = FALSE,
                     verbose = FALSE)
  result}

### Functions to extract data from OR_results function


measures = function(list){
   (data.frame(list$measure))
}

p.values = function(list){
  data.frame((list$p.value))
}


#### Applying functions and getting results
RESULTS = lapply(dfs, OR_results) # This gives a list of result tables
MEASURES = lapply(RESULTS, measures)
P.VALUES = lapply(RESULTS, p.values)

### Filling in results_df
results_df$PRS_group= rep(names(RESULTS), times = c(2,2,4,4,4,4))
results_df$Hormone_group = factor(unlist(lapply(MEASURES, row.names)), levels = c("very low", "low", "high", "very high"))
results_df$Odds_Ratio= unlist(lapply(MEASURES,"[",1))
results_df$C.I.95_Lower = unlist(lapply(MEASURES,"[",2))
results_df$C.I.95_Upper = unlist(lapply(MEASURES,"[",3))
results_df$Fisher_pval= unlist(lapply(P.VALUES,"[",2))
results_df$Chi2_pval = unlist(lapply(P.VALUES,"[",3))


#Writing to File

write.csv(results_df, file = "Odds_ratios.csv")	
print("Odds Ratio CSV saved")



############## Creating plot


OR_plot = results_df %>%
  ggplot(aes(y=Hormone_group, x=Odds_Ratio, label=Hormone_group)) +
  geom_point(size=1, shape=19) +
  geom_errorbarh(aes(xmin=C.I.95_Lower, xmax=C.I.95_Upper), height=.3) +
  geom_vline(xintercept=1, linetype='longdash') +
  facet_wrap(~factor(PRS_group, levels = c("low2", "high2", "very.low4", "low4", "high4", "very.high4")), ncol=2)

ggsave('OR_plot.png', OR_plot, device='png', width=5, units="in")

print("Plot Saved")


