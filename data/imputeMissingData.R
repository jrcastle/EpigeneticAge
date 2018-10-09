library("impute")
setwd("/home/jrca253/EpigeneticAge/data")

meth.file.in  = "meth_N_lt10_missing.txt"
meth.file.out = "meth_N_imputed.txt"


###########################################################################################
# LOAD DATA AND TRANSFORM INTO A MATRIX
###########################################################################################
print(paste("Loading ", meth.file.in, " ...", sep = ""))
meth.df.in = read.table(meth.file.in, header = TRUE, sep = '\t')
meth.matrix.in <- meth.df.in[ c(2:ncol(meth.df.in)) ]
meth.matrix.in <- data.matrix(meth.matrix.in)


###########################################################################################
# IMPUTE MISSING DATA AND TRANSFORM RESULTS INTO A DATAFRAME
###########################################################################################
print("Imputing ...")
meth.df.out <- data.frame(impute.knn(meth.matrix.in)$data)
meth.df.out$position = meth.df.in$position
meth.df.out <- meth.df.out[ colnames(meth.df.in) ]


###########################################################################################
# SAVE
###########################################################################################
print(paste("Saving ", meth.file.in, " ...", sep = ""))
write.table(meth.df.out, meth.file.out, sep = "\t", row.names=FALSE)

