# Running GEM on Body Mass Index (BMI) as the phenotype and testosterone (female) as the exposure
## This code should be run in R within the HPC environment 


### Terminal
#$ srun -p sunlab --pty bash # interactive node 


### open R
#$ ml R/4.0.3  # module load
#$ R           # open module 

# load packages --Again this is all done in the HPC terminal
rm(list=ls())
library(data.table)
library(crayon)
library(dplyr)
library(cli)


# using cov file and pheno file that were made in "Hormones.R"

pheno = fread("/data/Hormones/pheno_file_f.txt")
covar = fread("/data/Hormones/cov_file_f.txt")
names(pheno)
names(covar)

# Need to combine these
testo_GEM = merge(covar, pheno, by = "IID") 
dim(testo_GEM)
# [1] 249694     17

testo_GEM = testo_GEM[!is.na(testo_GEM$Testosterone),] # removing entries with no Testo data
dim(testo_GEM)
# [1] 198429     17


## Remove Relatedness
kin<-read.table("/projects/sunlab/UKB/Data_raw/ukb34031_rel_s488363.dat",header=T)
dim(kin) #107156      5
kin<-kin[kin$Kinship>=0.0884,] #2nd degree
dim(kin) #40230     5
length(unique(kin$ID1)) #36179
length(unique(kin$ID2)) #36159
ID<-c(kin$ID1,kin$ID2) #80460
dup<-names(table(ID)[table(ID)>1]) #9776
kin<-kin[(!(kin$ID1%in%dup))&(!(kin$ID2%in%dup)),] #27653

set.seed(1234)
kin$temp<-sample(1:2,dim(kin)[1],replace=T)
kin$exclude<-ifelse(kin$temp==1,kin$ID1,kin$ID2)
length(unique(kin$exclude)) #27653
exclude<-c(kin$exclude,dup)
length(exclude) # 37429 = 25653+9776
length(unique(exclude)) # 37429

testo_GEM<- subset(testo_GEM, !testo_GEM$IID %in% exclude)
dim(testo_GEM) # 181536     17


# Adding in UKB BMI data

BMI = fread("/data/UKB_BMI/Enrol_BMI.csv")
head(BMI)

# Adding primary care BMI data
pricare <- fread("/data/Longitudinal_BMI/primarycare_BMI_weight_height.txt",header=T,stringsAsFactors=F)
dim(pricare) #4157208       8
head(pricare)  
table(pricare$read_2)
#          229..   22A..   22K..      
#2972794  260751  489647  434016
table(pricare$read_3)
#            229..   22A..   22K..
# 1184414  735636 1156041 1081117
# "22K.." is the code for BMI. "22A.." is the code for Weight. "229.." is the code for Weight

pricare$event_dt<-as.Date(pricare$event_dt,format="%Y-%m-%d")
pricare<-pricare[!is.na(pricare$event_dt),]
dim(pricare) # 4153919   8
length((unique(pricare$eid))) # 213856


#exploring the data
head(sort(table(pricare$value1),decreasing=T)) 
head(sort(table(pricare$value2),decreasing=T)) # 4153813 are NA
head(sort(table(pricare$value3),decreasing=T))
# NA           Kg      cm   kg/m2      kg    27.7
#3902111   35045   27426   11908    6678    1737

# Selecting only BMI entries
pricare <- pricare[(pricare$read_2 == "22K.."&!is.na(pricare$read_2))|(pricare$read_3 == "22K.."&!is.na(pricare$read_3)),]
dim(pricare) # 1514188   8 This is different from the number reported before. Not sure why.

pricare$value1 <- as.numeric(pricare$value1)
pricare$value2 <- as.numeric(pricare$value2) # this turns non numeric values into NAs
pricare$value3 <- as.numeric(pricare$value3)# this turns non numeric values into NAs

pricare$dup <- ifelse(!is.na(pricare$value1) & !is.na(pricare$value3), 1,0) #check if any patients have two inputs on the same date
table(pricare$dup)# All 0. So no duplicates

#Creating BMI column from value1 and value3. Value2 is mostly NAs
pricare$bmi <- ifelse(is.na(pricare$value1), pricare$value3, pricare$value1)
pricare<-na.omit(pricare[,c("eid","event_dt","bmi")])
summary(pricare)
# Lowest BMI is 0.0, earliest event date is 1/1/1901. Going to get rid of likely incorrect entries
pricare<-pricare[pricare$bmi>=15&pricare$bmi<=79.5&pricare$event_dt>="1985-01-01"&pricare$event_dt<="2022-05-17",]
summary(pricare)





# Comparing IID in the 2 dfs
colnames(BMI)[c(1,4)] = c("IID", "BMI")
colnames(pricare)[c(1,3)] = c("IID", "BMI")
length(setdiff(pricare$IID, BMI$IID)) #0. So the pricare data isnt actually needed
length(setdiff( BMI$IID, pricare$IID)) # 293722


#Adding BMI data to covariate/ testo df
BMI_testo = merge(testo_GEM, BMI[,c("IID", "BMI")], by = "IID")
dim(BMI_testo) #  181536   

summary(BMI_testo$BMI)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  12.12   23.49   26.14   27.12   29.74   67.38     569

# Need to remove patients with 15<=BMI>=79.5  
BMI_testo = BMI_testo[BMI_testo$BMI>=15&BMI_testo$BMI<=79.5,] # 180949 patients left

#Getting rid of columns I don't need. Final file should have IID, FID, Phenotype, Covariates, and Exposure

BMI_testo = BMI_testo[,c(1:3, 5:14,16,18)]

names(BMI_testo)
#[1] "IID"          "FID.x"        "Age_enrolled" "PC1"          "PC2"
#[6] "PC3"          "PC4"          "PC5"          "PC6"          "PC7"
#[11] "PC8"          "PC9"          "PC10"         "Testosterone" "BMI"

colnames(BMI_testo)[2] = c("FID")

##### Adding BMI Polygenic Risk score (PRS) data 

load(file = "/data/ukb_core_data_52576.rda")
ls() # use this to find out what you just loaded. Name = bd

grep("26216",names(bd),value=T) # finding Standard BMI PRS data.  ""f.26216.0.0""
names(bd)[names(bd)=="f.26216.0.0"] = "PRS_BMI"

colnames(bd)[1] # looking for ID column [1] "f.eid"
colnames(bd)[1]= "IID"

BMI_testo = merge(BMI_testo, bd[,c("IID", "PRS_BMI")], by = "IID")
dim(BMI_testo) # 180949     16

summary(BMI_testo$PRS_BMI)
# Min. 1st Qu.  Median    Mean     3rd Qu.    Max.    NA's
#-4.5402 -0.8628 -0.2065 -0.2064  0.4511     3.8874     671


summary(BMI_testo$BMI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  15.02   23.49   26.14   27.12   29.74   67.38

summary(BMI_testo$Testosterone)
#   Min. 1st Qu.  Median    Mean   3rd Qu.    Max.
#  0.350   0.721   1.012   1.118   1.372    31.424


summary(BMI_testo$Age_enrolled)
# Min. 1st Qu.  Median    Mean     3rd Qu.    Max.
#  39.00   50.00   57.00   56.16   63.00    70.00

write.table(BMI_testo, file = "/data/GEM_testo_BMI/Female/BMI_testo_f.pheno", row.names = F, quote = F)

#Summary stats pre-analysis

##### Building BMI~PRSxTesto model

BMI_testo = fread("/data/GEM_testo_BMI/Female/BMI_testo_f.pheno")
head(BMI_testo)

model = lm(formula = BMI ~ Age_enrolled + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PRS_BMI * Testosterone, data = BMI_testo)

summary(model)

#####################################
#Viewing BMI vs PRS and BMI vs Testo

BMI_testo$PRS_group = ifelse(BMI_testo$PRS_BMI < -0.2065,"low","high" ) # Less than median = "low"
BMI_testo$Testo_group = ifelse(BMI_testo$Testosterone < 1.012,"low","high" ) # Less than median = "low"

df = group_by(BMI_testo,PRS_group, Testo_group)
df = summarise(df, BMI_median = median(BMI), count = n())

ggplot(data = subset(df,PRS_group!="NA"),                                      
       aes(x = PRS_group,
           y = BMI_median,
           fill = Testo_group)) +
  geom_bar(stat = "identity",
           position = position_dodge(0.75), width = 0.70)




ggplot(data = subset(BMI_testo,PRS_group!="NA"), aes(x = Testosterone,  y = BMI, color = PRS_group)) +
  geom_smooth(method = "lm")+
  labs(title="BMI by Testosterone in High and Low PRS groups")


# 4 quartiles 
BMI_testo$PRS_group = factor(ifelse(BMI_testo$PRS_BMI <= -0.8628,"very low", ifelse(BMI_testo$PRS_BMI > -0.8628 & BMI_testo$PRS_BMI <= -0.2065,"low", ifelse(BMI_testo$PRS_BMI > -0.2065 & BMI_testo$PRS_BMI <= 0.4511,"high", "very high"))), levels = c("very low", "low", "high", "very high")) # cutoffs = 1st Qu, median, 3rd Qu
class(BMI_testo$PRS_group)


BMI_testo$Testo_group = factor(ifelse(BMI_testo$Testosterone <= .721,"very low", ifelse(BMI_testo$Testosterone > 0.721 & BMI_testo$Testosterone <= 1.012,"low", ifelse(BMI_testo$Testosterone > 1.012 & BMI_testo$Testosterone <=1.372,"high", "very high"))), levels = c("very low", "low", "high", "very high")) # cutoffs = 1st Qu, median, 3rd Qu

df = group_by(BMI_testo,PRS_group, Testo_group)
df= summarise(df, BMI_median = median(BMI), count = n())


# Arbitrary cutoffs for Testo

BMI_testo$Testo_group = factor(ifelse(BMI_testo$Testosterone <=0.4,"very low", ifelse(BMI_testo$Testosterone > 0.4 & BMI_testo$Testosterone <= 1.0,"low", ifelse(BMI_testo$Testosterone > 1.0 & BMI_testo$Testosterone <=5,"high", "very high"))), levels = c("very low", "low", "high", "very high")) # cutoffs = 1st Qu, median, 3rd Qu

table(BMI_testo$Testo_group)
# very low       low      high very high 
#     4014     84653     91955       327 

#**** Results****

# Function to combine files
multmerge = function(path){
  filenames=list.files(path=path, full.names=TRUE)
  rbindlist(lapply(filenames, fread))
}

path <- "/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/out"
testo_BMI <- multmerge(path)

dim(testo_BMI) # 16622203      24
names(testo_BMI)  

testo_BMI$MAF = ifelse(testo_BMI$AF<0.5,testo_BMI$AF,1-testo_BMI$AF) # creating a MAF column
testo_BMI = testo_BMI[testo_BMI$MAF>0.01] # Filtering by MAF >0.01
dim(testo_BMI) # 9820365 25

marg_sig = testo_BMI[testo_BMI$robust_P_Value_Marginal<5e-08]  # P_Value_Marginal refers to the p_value of the main effect (aka G). This is essentially a GWAS p-value for each SNP
dim(marg_sig) # 12394   25

inter_sig = testo_BMI[testo_BMI$robust_P_Value_Interaction<5e-08] # P_Value_Interaction refers to the p_val of the interaction between G and the env (testosterone in this case). How much is the interaction between the SNP and testosterone associated with BMI outcome 
dim(inter_sig) # 28 25

joint_sig = testo_BMI[testo_BMI$robust_P_Value_Joint<5e-08] # P_Value_Joint refers to the join association of GxE + G + E.
dim(joint_sig) # 9674   25

write.table(testo_BMI, file = "testo_BMI_MAF01.txt", row.names = F, quote = F)


#### Plots made in FUMA

#For FUMA
data = fread("testo_BMI_MAF01.txt")

marginal = data[,c(1:5, 8, 9, 19)] 
names(marginal)
# [1] "SNPID"                   "CHR"
#[3] "POS"                     "Non_Effect_Allele"
#[5] "Effect_Allele"           "Beta_Marginal"
#[7] "robust_SE_Beta_Marginal" "robust_P_Value_Marginal"

interaction = data[,c(1:5,12,14,20)]
names(interaction)
#[1] "SNPID"                         "CHR"
#[3] "POS"                           "Non_Effect_Allele"
#[5] "Effect_Allele"                 "Beta_G-Testosterone"
#[7] "robust_SE_Beta_G-Testosterone" "robust_P_Value_Interaction"

joint = data[,c(1:5,11,13,21)]
names(joint)
#[1] "SNPID"                "CHR"                  "POS"
#[4] "Non_Effect_Allele"    "Effect_Allele"        "Beta_G"
#[7] "robust_SE_Beta_G"     "robust_P_Value_Joint"

colnames(marginal)[c(4:8)] = c("A0","A1","BETA","SE","P")
colnames(interaction)[c(4:8)] = c("A0","A1","BETA","SE","P")
colnames(joint)[c(4:8)] = c("A0","A1","BETA","SE","P")


write.table(marginal, file = "marginal_FUMA.txt", row.names = F, quote = F)
write.table(interaction, file = "interaction_FUMA.txt", row.names = F, quote = F)
write.table(joint, file = "joint_FUMA.txt", row.names = F, quote = F)

#Note for FUMA: there are 180596 samples

###############     PLINK
# Want to run GWAS for feamles vs BMI
# Pheno file is BMI_testo_f.pheno
# Going to extract covariates from that file

cov = fread("BMI_testo_f.pheno")
cov = cov[,c(1:13)] # removing pheno columns
write.table(cov, file = "covar.txt", row.names = F, quote = F)

#### After multmerge and testo_BMI has been created for plink results

dim(testo_BMI) # 9934690      15
testo_BMI = testo_BMI[testo_BMI$A1_FREQ>0.01]
dim(testo_BMI)
# [1] 9781051      15\
# N = 180596


##### **** Comparing Plink to GEM joint

plink = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/FUMA_BMI_F_plink/GenomicRiskLoci.txt")

joint = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/FUMA_testoxBMI_joint_f/GenomicRiskLoci.txt")
joint_sig = joint[joint$P<5E-8,]

plink = plink[,2:5] # has uniqID, rsID, chromosome,and position for the locus
joint = joint[,2:5]

# defining start and end at 250 kb on either side of POS
plink$start = plink$pos-250000
plink$end = plink$pos + 250000


joint$start = joint$pos-250000
joint$end = joint$pos+250000


# creating a function to make dfs of joint SNPS that are within plink Loci 

#find.snps = function(positions, ranges){
  positions[ranges,
    .(chr,start,end,POS=x.pos, plink_Locus = uniqID, rsID = x.rsID),
   on=.(chr,pos>=start,pos<=end),nomatch=0L]
}

# biocManager to see what loci overlap between plink and joint 

library(GenomicRanges,lib.loc="/projects/sunlab/R.lib") # I had to load a lot of other libraries individually to make this work. annoying

plink_ranges <- GRanges(
    seqnames = Rle(plink$chr),
    ranges = IRanges(start = plink$start, end = plink$end))
names(plink_ranges) = plink$rsID


joint_ranges <- GRanges(
    seqnames = Rle(joint$chr),
    ranges = IRanges(start = joint$start, end = joint$end))
names(joint_ranges) = joint$rsID

overlaps = subsetByOverlaps(joint_ranges, plink_ranges) # 106 overlaps
joint_unique = setdiff(joint$rsID, names(overlaps))
#  [1] "rs10915816"      "rs935376"        "rs10198543"      "rs7590658"
# [5] "3:93973756_CA_C" "rs374021540"     "rs9402104"       "rs10957605"
# [9] "rs76793042"      "rs61925678"      "rs1576655"       "rs3810291"
# [13] "rs16986608"      "rs73183246"

joint[joint$rsID %in% joint_unique,]

#### Comparing unique joint Pvals to marg and intersection. 

marg = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/marginal_FUMA.txt")
inter = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/interaction_FUMA.txt")
joint = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/joint_FUMA.txt")

snps = c("rs10915816", "rs935376", "rs10198543", "rs7590658", "3:93973756_CA_C", "rs374021540", "rs9402104", "rs10957605","rs76793042", "rs61925678", "rs1576655", "rs3810291", "rs16986608", "rs73183246")

joint_BP = joint[joint$SNPID %in% snps,c("SNPID","BETA","P")]
names(joint_BP)[2:3] = c("BETA_joint", "P_joint")

marg_BP = marg[marg$SNPID %in% snps,c("SNPID","BETA","P")]
names(marg_BP)[2:3] = c("BETA_marg", "P_marg")

inter_BP = inter[inter$SNPID %in% snps,c("SNPID","BETA","P")]
names(inter_BP)[2:3] = c("BETA_inter", "P_inter")

df = cbind(marg_BP, inter_BP[,c(2,3)], joint_BP[,c(2,3)])
write.csv(df, file = "joint_loci_comp.csv")

##### Comparing plink to GEM marginal
marg = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/marginal_FUMA.txt")
plink = fread("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/testo_BMI_plink.txt")

colnames(plink)[c(1,3,4)] = c("CHR", "SNPID", "A0")  #rename columns


comparison = function(df1,df2){
     df1 = df1[df1$P<5e-8] #only genome wide significant snps
     df2 = df2[df2$P<5e-8]
     merge_1 = merge(df1, df2, by = c("CHR", "POS", "A0", "A1"))
     merge_2 = merge(df1, df2, by.x = c("CHR", "POS", "A0", "A1"), by.y = c("CHR", "POS", "A1", "A0"))
     print(paste("Dimensions of first merge:", dim(merge_1)))
     print(paste("Dimensions of second merge:", dim(merge_2)))
     rbind(merge_1, merge_2)
     }


plink_marg = comparison(plink,marg)

#plots
png("/projects/sunlab/Students_Work/Amonae_work/GEM_testo_BMI/Female/plink_marg_corr.png",width = 1600, height = 1200, bg = "white",pointsize=24)
plot(plink_marg$BETA.x, plink_marg$BETA.y, pch = 19, col = "lightblue")

# Print Pearson correlation on plot
text(x= -1, y = 1, paste("Correlation:", round(cor(plink_marg$BETA.x, plink_marg$BETA.y), 4))) # 0.9999

dev.off()



