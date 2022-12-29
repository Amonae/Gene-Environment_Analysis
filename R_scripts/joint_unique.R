# remove everything in environment
rm(list=ls())

# This will create a csv file comparing Beta and P values between marginal and interaction results to joint results 
# Setting Paths

r.libs = "/projects/sunlab/R.lib"

# Load libraries

library(data.table,lib.loc= r.libs)
library(crayon,lib.loc= r.libs)
library(dplyr,lib.loc= r.libs)
library(cli,lib.loc= r.libs)
library(BiocGenerics,lib.loc= r.libs)
library(S4Vectors,lib.loc= r.libs)
library(IRanges,lib.loc= r.libs)
library(GenomeInfoDb,lib.loc= r.libs)
library(GenomicRanges,lib.loc= r.libs) 

# Loading data

marg_loci = fread("FUMA_results/FUMA_marg/GenomicRiskLoci.txt")  
joint_loci = fread("FUMA_results/FUMA_joint/GenomicRiskLoci.txt") 

marg = fread("marginal_FUMA.txt")
inter = fread("interaction_FUMA.txt")
joint = fread("joint_FUMA.txt")

# Subsetting Loci DFs
marg_loci = marg_loci[,2:5] # has uniqID, rsID, chromosome,and position for the locus
joint_loci = joint_loci[,2:5]

# defining start and end at 250 kb on either side of POS
marg_loci$start = marg_loci$pos-250000
marg_loci$end = marg_loci$pos + 250000


joint_loci$start = joint_loci$pos-250000
joint_loci$end = joint_loci$pos+250000

# biocManager to see what loci overlap between marg and joint
marg_ranges = GRanges(
    seqnames = Rle(marg_loci$chr),
    ranges = IRanges(start = marg_loci$start, end = marg_loci$end))
names(marg_ranges) = marg_loci$rsID


joint_ranges = GRanges(
    seqnames = Rle(joint_loci$chr),
    ranges = IRanges(start = joint_loci$start, end = joint_loci$end))
names(joint_ranges) = joint_loci$rsID

overlaps = subsetByOverlaps(joint_ranges, marg_ranges) 
joint_unique = setdiff(joint_loci$rsID, names(overlaps))

# Comparing unique joint Pvals to marg and intersection. 

joint_BP = joint[joint$SNPID %in% joint_unique,c("SNPID","BETA","P")]
names(joint_BP)[2:3] = c("BETA_joint", "P_joint")

marg_BP = marg[marg$SNPID %in% joint_unique,c("SNPID","BETA","P")]
names(marg_BP)[2:3] = c("BETA_marg", "P_marg")

inter_BP = inter[inter$SNPID %in% joint_unique,c("SNPID","BETA","P")]
names(inter_BP)[2:3] = c("BETA_inter", "P_inter")

df = cbind(marg_BP, inter_BP[,c(2,3)], joint_BP[,c(2,3)])
write.csv(df, file = "joint_unique_comp.csv")
