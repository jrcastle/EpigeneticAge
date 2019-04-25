setwd( '/home/jrca253/EpigeneticAge/data' )
library(gtools)

delim            <- '\t'
meth.file.in     <- "meth_N_gt10R.txt"
meth.file.out    <- "meth_N_gt10R_AddMissHorvCpGs.txt"
missing.cpg.list <- "MissingHorvathClockCpGs_N.txt"


###########################################################################################
# LOAD DATA
###########################################################################################

# Meth file
print( paste("Loading ", meth.file.in, " ...", sep = "") )
df.meth <- read.table(meth.file.in, header = TRUE, sep = delim)

# Missing CpG List
print( paste("Loading ", missing.cpg.list, " ...", sep = "") )
df.missing.cpg <- read.table(missing.cpg.list, header = TRUE)

print("Filling missing CpG data.frame with NA values ...")
for(name in colnames(df.meth)){
  if(name == "position"){ next }
  df.missing.cpg[, name] <- NA 
}


###########################################################################################
# APPEND AND SORT
###########################################################################################

# Append
print( paste("Appending missing CpG data.frame to ", meth.file.in, " ...", sep = "") )
df.meth <- rbind(df.meth, df.missing.cpg)

# Sort
print("Sorting ...")
df.meth <- df.meth[mixedorder(as.vector(df.meth$position)),]


###########################################################################################
# WRITE
###########################################################################################
print( paste("Writing ", meth.file.out, " ...", sep = "") )
write.table(df.meth, meth.file.out, sep = delim, row.names=FALSE)
