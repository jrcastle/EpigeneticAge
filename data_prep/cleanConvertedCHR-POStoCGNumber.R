library(gtools)
setwd("/home/jrca253/EpigeneticAge/data")

REORDER = FALSE
in.file  <- "tmp.csv"
out.file <- "meth_T_gt10R_AddMissHorvCpGs_KNT_CGHorvClock.csv"
delim = ','

print("Loading ...")
dat0 <- read.table(in.file, sep = delim, header = TRUE)
dat0$position <- as.character(dat0$position)

print("Cleaning ...")
for (i in 2:dim(dat0)[[2]] ){ 
    dat0[,i]=as.numeric(as.character(dat0[,i])) 
}

if( REORDER ){
    print("Reordering ...")
    dat0 <- dat0[ mixedorder( as.vector(dat0$position) ),]
}

print("Writing ...")
write.table(dat0, file = out.file, row.names = FALSE, sep = delim)