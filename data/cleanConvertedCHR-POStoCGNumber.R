out.file <- "meth_T_850k_CGNumber.txt"

dat0 <- read.csv("tmp.csv")

for (i in 2:dim(dat0)[[2]] ){ 
    dat0[,i]=as.numeric(as.character(dat0[,i])) 
}

write.csv(dat0, file = out.file, row.names = FALSE)