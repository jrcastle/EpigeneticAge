df.cov  = read.table("cov.txt", sep = '\t', header = TRUE)
df.meth = read.table("meth_noNA2.txt", sep = '\t', header = TRUE)

# Get Columns with NA values
df.tmp <- df.cov[, colSums(is.na(df.cov)) > 0]
samples.to.drop <- colnames(df.tmp)


df.cov[samples.to.drop]  <- list(NULL)
df.meth[samples.to.drop] <- list(NULL)

write.table(df.cov, "cov_noNA.txt", sep="\t", row.names = FALSE)
write.table(df.meth, "meth_noNA_addCov.txt", sep="\t", row.names = FALSE)