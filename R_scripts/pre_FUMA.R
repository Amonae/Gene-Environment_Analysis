### preparing GEM output for FUMA analysis

multmerge = function(path){
  filenames=list.files(path=path, full.names=TRUE)
  rbindlist(lapply(filenames, fread))
}

path = paste(getwd(), "out", sep = "/")
df = multmerge(path)
 

df$MAF = ifelse(df$AF<0.5,df$AF,1-df$AF) # creating a MAF column
df = df[df$MAF>0.01] # Filtering by MAF >0.01
dim(df) # 9820376 25

dim(df[df$robust_P_Value_Marginal<5e-08])  # 15713


dim(df[df$robust_P_Value_Interaction<5e-08]) # 18


dim(df[df$robust_P_Value_Joint<5e-08]) #13299

write.table(df, file = "results_MAF01.txt", row.names = F, quote = F)


#### Plots made in FUMA

#For FUMA

marginal = df[,c(1:5, 8, 9, 19)] 

interaction = df[,c(1:5,12,14,20)]


joint = df[,c(1:5,11,13,21)]


colnames(marginal)[c(4:8)] = c("A0","A1","BETA","SE","P")
colnames(interaction)[c(4:8)] = c("A0","A1","BETA","SE","P")
colnames(joint)[c(4:8)] = c("A0","A1","BETA","SE","P")


write.table(marginal, file = "marginal_FUMA.txt", row.names = F, quote = F)
write.table(interaction, file = "interaction_FUMA.txt", row.names = F, quote = F)
write.table(joint, file = "joint_FUMA.txt", row.names = F, quote = F)
