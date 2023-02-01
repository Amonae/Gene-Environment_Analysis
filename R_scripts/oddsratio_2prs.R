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
suppressPackageStartupMessages(library(epitools, lib.loc = "/projects/sunlab/Students_Work/Amonae_work/R.lib")) # Note different lib loc


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

pheno_file$PRS_group = factor(ifelse(pheno_file$PRS_outcome<= PRS[3],"low", "high"), levels = c("low", "high"), ordered = T) # cutoffs = median


pheno_file$hormone_group = factor(ifelse(pheno_file$hormone <= hormone[2],"very low", ifelse(pheno_file$hormone >hormone[2] & pheno_file$hormone <= hormone[3],"low", ifelse(pheno_file$hormone > hormone[3] & pheno_file$hormone <=hormone[5],"high", "very high"))), levels = c("very low", "low", "high", "very high"), ordered = T) # cutoffs = 1st Qu, median, 3rd Qu 

low_PRS = pheno_file[pheno_file$PRS_group == "low",] # for comparing 2 PRS levels and 4 hormone levels

high_PRS = pheno_file[pheno_file$PRS_group == "high",] # for comparing 2 PRS levels and 4 hormone levels

dfs =  list(low_PRS, high_PRS)
names(dfs) = c("low_PRS", "high_PRS")

### Creating empty df which will be filled with results

results_df = data.frame(matrix(nrow = 8, ncol = 7))
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

print(RESULTS)

### Filling in results_df
results_df$PRS_group= rep(names(RESULTS), times = c(4,4))
results_df$PRS_group = factor(results_df$PRS_group, levels = c("low_PRS", "high_PRS"), ordered = T)
results_df$Hormone_group = factor(unlist(lapply(MEASURES, row.names)), levels = c("very low", "low", "high", "very high"), ordered = T)
results_df$Odds_Ratio= unlist(lapply(MEASURES,"[",1))
results_df$C.I.95_Lower = unlist(lapply(MEASURES,"[",2))
results_df$C.I.95_Upper = unlist(lapply(MEASURES,"[",3))
results_df$Fisher_pval= unlist(lapply(P.VALUES,"[",2))
results_df$Chi2_pval = unlist(lapply(P.VALUES,"[",3))

write.csv(results_df, file = "Odds_ratios_2prs.csv")	
print("Odds Ratio CSV saved")

############## Creating forest plot


OR_FP = ggplot(data=results_df,
           aes(x = Hormone_group,y = Odds_Ratio, ymin = C.I.95_Lower, ymax = C.I.95_Upper ))+
  geom_pointrange(aes(col=Hormone_group))+
  geom_hline(aes(fill=Hormone_group),yintercept =1, linetype=2)+
  xlab('Hormone group')+ ylab("Odds Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=C.I.95_Lower, ymax=C.I.95_Upper,col=Hormone_group),width=0.5,cex=1)+ 
  facet_wrap(~PRS_group,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()

ggsave('OR_Forestplot_2PRS.png', OR_FP, device='png', width=5, units="in")

print("Forest Plot Saved")


############## Bar plot
	

OR_BP = ggplot(data = results_df,                                      
       aes(x = PRS_group,
           y = Odds_Ratio,
           fill = Hormone_group)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.75), width = 0.70)

ggsave('OR_Barplot_2PRS.png', OR_BP, device='png', width=5, units="in")

print("Bar Plot Saved")
