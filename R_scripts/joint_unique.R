# remove everything in environment
rm(list=ls())

print(getwd())
# This will create a csv file comparing Beta and P values between marginal and interaction results to joint results 
# Setting Paths

r.libs = "/projects/sunlab/R.lib"

# Load libraries

suppressPackageStartupMessages(library(data.table,lib.loc= r.libs))
suppressPackageStartupMessages(library(crayon,lib.loc= r.libs))
suppressPackageStartupMessages(library(dplyr,lib.loc= r.libs))
suppressPackageStartupMessages(library(cli,lib.loc= r.libs))
suppressPackageStartupMessages(library(BiocGenerics,lib.loc= r.libs))
suppressPackageStartupMessages(library(S4Vectors,lib.loc= r.libs))
suppressPackageStartupMessages(library(IRanges,lib.loc= r.libs))
suppressPackageStartupMessages(library(GenomeInfoDb,lib.loc= r.libs))
suppressPackageStartupMessages(library(GenomicRanges,lib.loc= r.libs)) 

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
print(paste(length(joint_unique), "unique joint SNPS"))

# Comparing unique joint Pvals to marg and intersection. 

joint_BP = joint[joint$SNPID %in% joint_unique,c("SNPID","CHR", "POS","BETA","P")]
names(joint_BP)[4:5] = c("BETA_joint", "P_joint")
print("Joint Beta and Pvals")
print(head(joint_BP))


marg_BP = marg[marg$SNPID %in% joint_unique,c("SNPID","BETA","P")]
names(marg_BP)[2:3] = c("BETA_marg", "P_marg")
print("Marginal Beta and Pvals")
print(head(marg_BP))

inter_BP = inter[inter$SNPID %in% joint_unique,c("SNPID","BETA","P")]
names(inter_BP)[2:3] = c("BETA_inter", "P_inter")
print("Interaction Beta and Pvals")
print(head(inter_BP))


df = cbind(joint_BP, inter_BP[,c(2,3)], marg_BP[,c(2,3)])
write.csv(df, file = "joint_unique_comp.csv")
print("joint_unique_comp saved")
# Comparing all joint_loci Pvals to marg and intersection. 

joint_BP.2 = joint[joint$SNPID %in% joint_loci$rsID,c("SNPID","CHR", "POS","BETA","P")]
names(joint_BP.2)[4:5] = c("BETA_joint", "P_joint")

marg_BP.2 = marg[marg$SNPID %in% joint_loci$rsID,c("SNPID","BETA","P")]
names(marg_BP.2)[2:3] = c("BETA_marg", "P_marg")

inter_BP.2 = inter[inter$SNPID %in% joint_loci$rsID,c("SNPID","BETA","P")]
names(inter_BP.2)[2:3] = c("BETA_inter", "P_inter")

df.2 = cbind(joint_BP.2, inter_BP.2[,c(2,3)], marg_BP.2[,c(2,3)])
write.csv(df.2, file = "joint_all_comp.csv")
print("joint_unique_all saved")

q()