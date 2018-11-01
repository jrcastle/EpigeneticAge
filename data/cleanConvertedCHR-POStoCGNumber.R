dat0 <- read.csv("tmp.csv")

for (i in 2:dim(dat0)[[2]] ){ 
    dat0[,i]=as.numeric(as.character(dat0[,i])) 
}

write.csv(dat0, file = "meth_CGNumber.csv", row.names = FALSE)