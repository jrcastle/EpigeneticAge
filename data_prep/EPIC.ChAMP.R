#ChAMP on 
source("https://bioconductor.org/biocLite.R")
library(devtools)
#install_github("JoshuaTian/MyChAMP")
library(ChAMP)


idat.dir    = "/Volumes/ChunyanLab/forJames/idat2"
analysis.dir = "/Volumes/ChunyanLab/forJames/CHAMP_output"
setwd( '/Users/jrca253/Documents/EpigeneticAge/test_code/ImputeCheck_850k')

#load data==========================================================================================
myLoad <- champ.load(directory = idat.dir, methValue="B", arraytype="EPIC")
myLoad$pd$Race = gsub("\\..*", "", myLoad$pd$Match_Group)
myLoad$pd$Age = as.numeric(gsub("[EA]A\\.([0-9]+)\\..*", "\\1", myLoad$pd$Match_Group))
myLoad$pd$Subtype = gsub(".*\\.([a-zA-Z]+)", "\\1", myLoad$pd$Match_Group)
myLoad$pd$Subtype<-ifelse(myLoad$pd$Subtype=='NA',NA,myLoad$pd$Subtype)

CpG.GUI(arraytype="EPIC")
write.csv(myLoad$beta, file="M-value_11_sample.csv")

#QC=================================================================================================
#quality check
champ.QC(beta=myLoad$beta, pheno=myLoad$pd$Sample_Group) # Alternatively QC.GUI(arraytype="EPIC")
QC.GUI(beta=myLoad$beta, pheno=myLoad$pd$Sample_Group, arraytype="EPIC")

#normalization======================================================================================
myNorm <- champ.norm(beta=myLoad$beta,method="BMIQ", arraytype="EPIC", core = 6)
myNorm <- champ.norm(beta=myLoad$beta,method="BMIQ",rgSet=myLoad$rgSet,mset=myLoad$mset, arraytype="EPIC")
myNorm <- (myNorm+0.00001)/(1.00001-myNorm)
QC.GUI(beta=myNorm, arraytype="EPIC")


na_omit_B_seq = na.omit(B_seq)
rownames(na_omit_B_seq) = na_omit_B_seq[['CpG']]
na_omit_B_seq = na_omit_B_seq[c(-1,-2)]
na_omit_B_seq = na_omit_B_seq/100
myNorm_seq = champ.norm(beta = na_omit_B_seq, method="BMIQ", arraytype="EPIC")

beta_mseq <- myLoad$beta[which(rownames(myLoad$beta) %in% rownames(na_omit_B_seq)),]
myNorm_mseq <- champ.norm(beta=beta_mseq,method="BMIQ", arraytype="EPIC")
#save the normalized data
# Anno = read.table("../data/EPIC_Array_Data-IU/2_VariableAnnotation.txt", skip=7, header=T, sep="\t")
# BDAT = merge(myNorm, Anno, by.x=0, by.y=1)
# write.csv(BDAT, file="ChAMP.Normalized.beta.csv", row.names=F)

#SVD================================================================================================
champ.SVD(beta=myNorm, pd=myLoad$pd)
#none of the race or pair-setting is significant

# myCombat = champ.runCombat(beta=myNorm, pd=myLoad$pd, batchname=c("Match_Group"))

#DM=================================================================================================
#DM=================================================================================================
#DM=================================================================================================


#DM1================================================================================================
#1. White vs Black in K-----------------------------------------------------------------------------
#matched
PHE = myLoad$pd$Race
PHE[myLoad$pd$Sample_Group != "K"] = NA
table(PHE) #3 vs 3
PHE[myLoad$pd$Sample_Name %in% c("K106158", "K102645")] = NA
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)])
pair.EAvAA = c(1, 1, NA,2,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)

myDMP1 <-   champ.PairedDMP(beta=myNorm, pair =pair.EAvAA, pheno=PHE, arraytype="EPIC", compare.group=c("EA","AA"), adjPVal=1)
DMP.GUI(DMP=myDMP1, beta=myNorm, pheno=PHE)
write.csv(myDMP1, file="/Users/nanlin/Desktop/DM_sites_2EAvs2AA_all.csv")

#control for covariates
PHE = cbind(myLoad$pd$Race,myLoad$pd$Age,myLoad$pd$Slide)
PHE[which(myLoad$pd$Sample_Group != "K"),] = NA
table(PHE) #3 vs 3
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], na.omit(PHE))

myDMP1 <- champ.DMP(beta=myNorm, pheno=PHE, arraytype="EPIC", compare.group=c("EA","AA"), adjPVal=1)
DMP.GUI(DMP=myDMP1, beta=myNorm, pheno=PHE)
write.csv(myDMP1, file="/Users/nanlin/Desktop/DM_sites_3EAvs3AA_all_controled.csv")


#DM2================================================================================================
#2. T vs N : tumor vs adjacent----------------------------------------------------------------------
#same probe number as truseq with Truseq
PHE = myLoad$pd$Sample_Group
PHE[PHE=="K"] = NA
table(PHE) #4vs4
PHE      = c(NA,'T',NA,'N','N','T','N','T','T',NA,'N')
pair.TvN = c(NA,'4',NA,'2','1','2','3','1','3',NA,'4')
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvN[!is.na(PHE)])

myDMP2 <- champ.PairedDMP(beta=myNorm_mseq, pheno=PHE, pair = pair.TvN, arraytype="EPIC", compare.group=c("T","N"), adjPVal=1)
DMP.GUI(DMP=myDMP2, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP2, file="/Users/chunyanhelab/Downloads/matchedprobe_Array_DM_sites_4Tvs4N.csv")


myNorm_seq = myNorm_seq[,c(3,11,1,5,4,9,6,8,10,2,7)]
PHE = myLoad$pd$Sample_Group
PHE[PHE=="K"] = NA
table(PHE) #4vs4
PHE      = c(NA,'T',NA,'N','N','T','N','T','T',NA,'N')
pair.TvN = c(NA,'4',NA,'2','1','2','3','1','3',NA,'4')
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvN[!is.na(PHE)])

myDMP2 <- champ.PairedDMP(beta=myNorm_seq, pheno=PHE, pair = pair.TvN, arraytype="EPIC", compare.group=c("T","N"), adjPVal=1)

write.csv(myDMP2, file="/Users/chunyanhelab/Downloads/matchedprobe_seq_DM_sites_4Tvs4N.csv")

#matched with Truseq
PHE = myLoad$pd$Sample_Group
PHE[PHE=="K"] = NA
table(PHE) #4vs4
PHE      = c(NA,'T',NA,'N','N','T','N','T','T',NA,'N')
pair.TvN = c(NA,'4',NA,'2','1','2','3','1','3',NA,'4')
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvN[!is.na(PHE)])

myDMP2 <- champ.PairedDMP(beta=myNorm, pheno=PHE, pair = pair.TvN, arraytype="EPIC", compare.group=c("T","N"), adjPVal=1)
DMP.GUI(DMP=myDMP2, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP2, file="/Volumes/My Passport/matched_Truseq_DM_sites_4Tvs4N.csv")

#matched
PHE = myLoad$pd$Sample_Group
PHE[PHE=="K"] = NA
table(PHE) #5vs5
PHE      = c(NA,NA,'T',NA,NA,NA,'N','N','T','N','T',NA,'T',NA,NA,'N')
pair.TvN = c(NA,NA,'4',NA,NA,NA,'2','1','2','3','1',NA,'3',NA,NA,'4')
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvN[!is.na(PHE)])

myDMP2 <- champ.PairedDMP(beta=myNorm, pheno=PHE, pair = pair.TvN, arraytype="EPIC", compare.group=c("T","N"), adjPVal=1)
DMP.GUI(DMP=myDMP2, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP2, file="/Volumes/My Passport/DM_sites_5Tvs5N_all.csv")
#unmatched
PHE = cbind(myLoad$pd$Sample_Group, myLoad$pd$Race,myLoad$pd$Age,myLoad$pd$Slide)
PHE[which(myLoad$pd$Sample_Group == "K"),] = NA
table(PHE) #5 vs 5
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], na.omit(PHE))

myDMP2 <- champ.DMP(beta=myNorm, pheno=PHE, arraytype="EPIC", compare.group=c("T","N"), adjPVal=1)

DMP.GUI(DMP=myDMP1, beta=myNorm, pheno=PHE)
write.csv(myDMP2, file="/Volumes/My Passport/DM_sites_5Tvs5N_all_unmatched.csv")
#DM3================================================================================================
#3.T vs K : tumor vs normal-------------------------------------------------------------------------
#matched probe with Truseq
PHE = myLoad$pd$Sample_Group
PHE[PHE=="N"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','T060')] = NA
table(PHE) #4vs4
PHE      = c('K',NA,'K',NA,NA,'T',NA,'T','T','K',NA)
pair.TvK = c('3',NA,'2',NA,NA,'2',NA,'1','3','1',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvK[!is.na(PHE)])
cbind(myLoad$pd$Sample_Name, PHE)

myDMP3 <- champ.PairedDMP(beta=myNorm_mseq, pheno=PHE, pair = pair.TvK, arraytype="EPIC", compare.group=c("T","K") , adjPVal=1)

# DMP.GUI(DMP=myDMP3, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP3, file="/Users/chunyanhelab/Downloads/matchedprobe_array_DM_sites_3Tvs3K.csv")

PHE = myLoad$pd$Sample_Group
PHE[PHE=="N"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','T060')] = NA
table(PHE) #4vs4
PHE      = c('K',NA,'K',NA,NA,'T',NA,'T','T','K',NA)
pair.TvK = c('3',NA,'2',NA,NA,'2',NA,'1','3','1',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvK[!is.na(PHE)])
cbind(myLoad$pd$Sample_Name, PHE)

myDMP3 <- champ.PairedDMP(beta=myNorm_seq, pheno=PHE, pair = pair.TvK, arraytype="EPIC", compare.group=c("T","K") , adjPVal=1)

# DMP.GUI(DMP=myDMP3, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP3, file="/Users/chunyanhelab/Downloads/matchedprobe_seq_DM_sites_3Tvs3K.csv")

#matched with Truseq
PHE = myLoad$pd$Sample_Group
PHE[PHE=="N"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','T060')] = NA
table(PHE) #4vs4
PHE      = c('K',NA,'K',NA,NA,'T',NA,'T','T','K',NA)
pair.TvK = c('3',NA,'2',NA,NA,'2',NA,'1','3','1',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvK[!is.na(PHE)])
cbind(myLoad$pd$Sample_Name, PHE)

myDMP3 <- champ.PairedDMP(beta=myNorm, pheno=PHE, pair = pair.TvK, arraytype="EPIC", compare.group=c("T","K") , adjPVal=1)

# DMP.GUI(DMP=myDMP3, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP3, file="/Volumes/My Passport/matched_Truseq_DM_sites_3Tvs3K.csv")
#matched
PHE = myLoad$pd$Sample_Group
PHE[PHE=="N"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','T060')] = NA
table(PHE) #4vs4
pair.TvK = c(NA,'3',NA,'4','2',NA,NA,NA,'2',NA,'1',NA,'3','1','4',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.TvK[!is.na(PHE)])
cbind(myLoad$pd$Sample_Name, PHE)

myDMP3 <- champ.PairedDMP(beta=myNorm, pheno=PHE, pair = pair.TvK, arraytype="EPIC", compare.group=c("T","K") , adjPVal=1)

# DMP.GUI(DMP=myDMP3, beta=myNorm, pheno=myLoad$pd$Sample_Group)
write.csv(myDMP3, file="/Volumes/My Passport/DM_sites_4Tvs4K.new.csv")

#unmatched
PHE = cbind(myLoad$pd$Sample_Group, myLoad$pd$Race,myLoad$pd$Age,myLoad$pd$Slide)
PHE[which(myLoad$pd$Sample_Group == "N"),] = NA
table(PHE) #5 vs 6
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], na.omit(PHE))

myDMP3 <- champ.DMP(beta=myNorm, pheno=PHE, arraytype="EPIC", compare.group=c("T","K"), adjPVal=1)
DMP.GUI(DMP=myDMP3, beta=myNorm, pheno=PHE)
write.csv(myDMP3, file="/Volumes/My Passport/DM_sites_5Tvs6K_unmatched.csv")

#DM4================================================================================================
#N vs K : adjacent vs normal------------------------------------------------------------------------
#matched probes with Truseq
PHE = myLoad$pd$Sample_Group
PHE[PHE=="T"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','N060')] = NA
table(PHE) #4vs4
PHE      = c('K',NA,'K','N','N',NA,'N',NA,NA,'K',NA)
pair.NvK = c('3',NA,'2','2','1',NA,'3',NA,NA,'1',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.NvK[!is.na(PHE)])

myDMP4 <- champ.PairedDMP(beta=myNorm_mseq, pheno=PHE, pair = pair.NvK, arraytype="EPIC", compare.group=c("K","N"), adjPVal=1)
# DMP.GUI(DMP=myDMP, beta=myNorm, pheno=myLoad$pd$Sample_Group)

write.csv(myDMP4, file="/Users/chunyanhelab/Downloads/matchedprobe_array_DM_sites_3Kvs3N.csv")

PHE = myLoad$pd$Sample_Group
PHE[PHE=="T"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','N060')] = NA
table(PHE) #4vs4
PHE      = c('K',NA,'K','N','N',NA,'N',NA,NA,'K',NA)
pair.NvK = c('3',NA,'2','2','1',NA,'3',NA,NA,'1',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.NvK[!is.na(PHE)])

myDMP4 <- champ.PairedDMP(beta=myNorm_seq, pheno=PHE, pair = pair.NvK, arraytype="EPIC", compare.group=c("K","N"), adjPVal=1)
# DMP.GUI(DMP=myDMP, beta=myNorm, pheno=myLoad$pd$Sample_Group)

write.csv(myDMP4, file="/Users/chunyanhelab/Downloads/matchedprobe_seq_DM_sites_3Kvs3N.csv")

#matched with Truseq
PHE = myLoad$pd$Sample_Group
PHE[PHE=="T"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','N060')] = NA
table(PHE) #4vs4
PHE      = c('K',NA,'K','N','N',NA,'N',NA,NA,'K',NA)
pair.NvK = c('3',NA,'2','2','1',NA,'3',NA,NA,'1',NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.NvK[!is.na(PHE)])

myDMP4 <- champ.PairedDMP(beta=myNorm, pheno=PHE, pair = pair.NvK, arraytype="EPIC", compare.group=c("K","N"), adjPVal=1)
# DMP.GUI(DMP=myDMP, beta=myNorm, pheno=myLoad$pd$Sample_Group)

write.csv(myDMP4, file="/Volumes/My Passport/matched_Truseq_DM_sites_3Kvs3N.csv")
#matched
PHE = myLoad$pd$Sample_Group
PHE[PHE=="T"] = NA
PHE[myLoad$pd$Sample_Name %in% c("K106158",'K102484','N060')] = NA
table(PHE) #4vs4
pair.NvK = c(NA,'3',NA,'4','2',NA,'2','1',NA,'3',NA,'4',NA,'1',NA,NA)
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], PHE[!is.na(PHE)],pair.NvK[!is.na(PHE)])

myDMP4 <- champ.PairedDMP(beta=myNorm, pheno=PHE, pair = pair.NvK, arraytype="EPIC", compare.group=c("K","N"), adjPVal=1)
# DMP.GUI(DMP=myDMP, beta=myNorm, pheno=myLoad$pd$Sample_Group)

write.csv(myDMP4, file="/Volumes/My Passport/DM_sites_4Kvs4N.csv")

#unmatched
PHE = cbind(myLoad$pd$Sample_Group, myLoad$pd$Race,myLoad$pd$Age,myLoad$pd$Slide)
PHE[which(myLoad$pd$Sample_Group == "T"),] = NA
table(PHE) #5 vs 6
cbind(myLoad$pd$Sample_Name[!is.na(PHE)], na.omit(PHE))

myDMP4 <- champ.DMP(beta=myNorm, pheno=PHE, arraytype="EPIC", compare.group=c("N","K"), adjPVal=1)
DMP.GUI(DMP=myDMP1, beta=myNorm, pheno=PHE)
write.csv(myDMP4, file="/Volumes/My Passport/DM_sites_6Kvs5N_unmatched.csv")

#====================================================================================================
#DMR
PHE = myLoad$pd$Sample_Group
PHE[PHE=="K"] = NA
PHE = na.omit(PHE)
TvN.Array.DMR<-champ.DMR(beta=myNorm[,c(2,4,5,6,7,8,9,11)],
          pheno=PHE,
          arraytype="EPIC",
          method = "Bumphunter",
          minProbes=3,
          adjPvalDmr=0.05,
          cores=3,
          minDmrSep=100 ,
          ## following parameters are specifically for Bumphunter method.
          maxGap=300,
          cutoff=NULL,
          pickCutoff=TRUE,
          smooth=TRUE,
          smoothFunction=loessByCluster,
          useWeights=FALSE,
          permutations=NULL,
          B=250,
          nullMethod="bootstrap")

write.csv(TvN.Array.DMR, file="TvN_Array_DMR.csv")


PHE = c('K','K','T','T','T','K')
TvK.Array.DMR<-champ.DMR(beta=myNorm[,c(1,3,6,8,9,10)],
                         pheno=PHE,
                         arraytype="EPIC",
                         method = "Bumphunter",
                         minProbes=3,
                         adjPvalDmr=0.05,
                         cores=3,
                         minDmrSep=100 ,
                         ## following parameters are specifically for Bumphunter method.
                         maxGap=300,
                         cutoff=NULL,
                         pickCutoff=TRUE,
                         smooth=TRUE,
                         smoothFunction=loessByCluster,
                         useWeights=FALSE,
                         permutations=NULL,
                         B=250,
                         nullMethod="bootstrap")

write.csv(TvK.Array.DMR, file="TvK_Array_DMR.csv")


PHE = c('K','K','N','N','N','K')
NvK.Array.DMR<-champ.DMR(beta=myNorm[,c(1,3,4,5,7,10)],
                         pheno=PHE,
                         arraytype="EPIC",
                         method = "Bumphunter",
                         minProbes=3,
                         adjPvalDmr=0.05,
                         cores=3,
                         minDmrSep=100,
                         ## following parameters are specifically for Bumphunter method.
                         maxGap=300,
                         cutoff=NULL,
                         pickCutoff=TRUE,
                         smooth=TRUE,
                         smoothFunction=loessByCluster,
                         useWeights=FALSE,
                         permutations=NULL,
                         B=250,
                         nullMethod="bootstrap")

write.csv(NvK.Array.DMR, file="NvK_Array_DMR.csv")

citation(package = "ChAMP")
citation(package = "DSS")
citation(package = "clusterProfiler")
citation(package = "missMethyl")
