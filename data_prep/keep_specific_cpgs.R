setwd("/home/jrca253/EpigeneticAge/data")
cpg.file      <- "WhiteBlackDiffMethCpGs.txt"
meth.file.in  <- "meth_K_AllCpGs.txt"
meth.file.out <- "meth_K_WhiteBlackDiffMethCpGs.txt"

df.cpg  <- read.table(cpg.file, header = FALSE)
colnames(df.cpg) <- c("CpGs")
df.cpg$CpGs <- as.character(df.cpg$CpGs)

print("Reading in file ...")
df.meth <- read.table(meth.file.in, header = TRUE, sep = '\t')
df.meth$position <- as.character(df.meth$position)
print("Reducing ...")
df.meth <- df.meth[ which(df.meth$position %in% df.cpg$CpGs), ]

print("Writing file ...")
write.table(
  df.meth,
  file = meth.file.out,
  sep = "\t", 
  row.names=FALSE
)