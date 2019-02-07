rm(list=ls()); gc();
setwd("/home/jrca253/EpigeneticAge/data")
#setwd("/Users/jrca253/Documents/EpigeneticAge/test_code/imputationTests")
library(tidyr)
library(gtools)
library(Rcpp)
sourceCpp("sameSampleImpute.cpp")

NEARBY.LIMIT      <- 2000 # bp
BPDIST.WEIGHT     <- TRUE

delim             <- '\t'
meth.file.in      <- "meth_K_gt10R_AddMissHorvCpGs_KNT_KnnImp.txt"
meth.file.out     <- "meth_K_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL.txt"
bad.cpg.list      <- "CpGsWithGt50pctMissingness_K.txt"


###########################################################################################
# LOAD/MERGE DATA
###########################################################################################
print(paste("Loading ", meth.file.in, " ...", sep = ""))
meth.df.in = read.table(meth.file.in, header = TRUE, sep = delim)

print(paste("Loading ", bad.cpg.list, " ...", sep = ""))
bad.cpg.list = read.table(bad.cpg.list, header = TRUE, sep = delim)

print("Merging data frames ...")
meth.df.in <- rbind(meth.df.in, bad.cpg.list)

print("Reordering")
meth.df.in <- meth.df.in[ mixedorder( as.vector(meth.df.in$position) ),]

###########################################################################################
# ADD EXTRA COLUMNS AND GET LIST OF BAD CPGS
###########################################################################################
meth.df.in$tmp   <- as.character(meth.df.in$position)
meth.df.in       <- separate(data = meth.df.in, col = tmp, into = c("chr", "bppos"), sep = "\\:")
meth.df.in$bppos <- as.numeric( meth.df.in$bppos )
meth.df.in$chr   <- as.numeric(gsub("[a-zA-Z ]", "", meth.df.in$chr))

bad.cpg.list <- as.character( bad.cpg.list$position )
bad.cpg.rows <- as.integer(which(meth.df.in$position %in% bad.cpg.list))


###########################################################################################
# IMPUTE C++
###########################################################################################
print("Imputing ...")
ptm <- proc.time()
sameSampleImpute(bad.cpg.list, bad.cpg.rows, meth.df.in, NEARBY.LIMIT, BPDIST.WEIGHT)
print("IMPUTATION FINISHED!")
print("Process time:")
print(proc.time() - ptm)

meth.df.in$chr <- NULL
meth.df.in$bppos <- NULL

print("Dropping CpGs where imputation failed ...")
failed.count <- nrow( meth.df.in[!complete.cases(meth.df.in),] )
if( failed.count > 0 ){
  print( paste("Failed to impute ", failed.count, " CpGs ...", sep = "") )
}
meth.df.in <- meth.df.in[complete.cases(meth.df.in),]


###########################################################################################
# WRITE
###########################################################################################
print( paste("Writing ", meth.file.out, " ...", sep = "") )
write.table(meth.df.in, meth.file.out, sep = delim, row.names=FALSE)
