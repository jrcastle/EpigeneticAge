library("impute")
setwd("/home/jrca253/EpigeneticAge/data")

CLEAN             <- TRUE
MISSINGNESS.LIMIT <- 0.01 # Default is 0.50 
maxp              <- 1500 # Default 1500, 150000 used for N
knn               <- 10   # Default is 10
options(expressions = 20000)

delim         <- '\t'
meth.file.in  <- "meth_K_gt10R_AddMissHorvCpGs_KNT.txt"
meth.file.out <- "TMP.txt"
bad.cpg.out   <- "TMP2.txt"

###########################################################################################
# LOAD DATA AND TRANSFORM INTO A MATRIX
###########################################################################################
print(paste("Loading ", meth.file.in, " ...", sep = ""))
meth.df.in = read.table(meth.file.in, header = TRUE, sep = delim)

if( CLEAN ){
    limit = as.integer( (ncol(meth.df.in) - 1) * MISSINGNESS.LIMIT )

    print( paste("Saving list of CpGs with >50% missing entries to ", bad.cpg.out, " ...", sep = "") )
    bad.cpgs <- meth.df.in[ rowSums(is.na(meth.df.in)) > limit, , drop=FALSE]
    write.table(bad.cpgs, bad.cpg.out, row.names=FALSE, sep = delim)
    rm(bad.cpgs); gc()

    print("Cleaning so that no rows have >50% missing entries")
    meth.df.in <- meth.df.in[ rowSums(is.na(meth.df.in)) <= limit, ]
}

meth.matrix.in <- meth.df.in[ c(2:ncol(meth.df.in)) ]
meth.matrix.in <- data.matrix(meth.matrix.in)


###########################################################################################
# IMPUTE MISSING DATA AND TRANSFORM RESULTS INTO A DATAFRAME
###########################################################################################
print("Imputing ...")
meth.df.out <- data.frame(impute.knn(meth.matrix.in, k = knn, maxp = maxp)$data)
meth.df.out$position = meth.df.in$position
meth.df.out <- meth.df.out[ colnames(meth.df.in) ]


###########################################################################################
# SAVE
###########################################################################################
print(paste("Saving ", meth.file.out, " ...", sep = ""))
write.table(meth.df.out, meth.file.out, sep = delim, row.names=FALSE)
