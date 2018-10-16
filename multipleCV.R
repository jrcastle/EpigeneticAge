setwd("/home/jrca253/EpigeneticAge")
library(glmnet)

cov.train  = "data/cov_K_train.txt"
meth.train = "data/meth_K_cpgs_in_KNT_imputed_train.txt"

alpha <- 0.5
adult.age <- 20
number.cvs <- 100

##########################################################################
# AGE TRANSFORM FUNCTION
##########################################################################
transform.age <- function(age){
  x = 0
  if(age <= adult.age){
    x = log(age + 1) - log(adult.age + 1)
  }else{
    x = (age - adult.age) / (adult.age + 1)
  }
  return(x)
}

transform.age.inverse <- function(tage){
  x = 0
  if(tage <= 0){
    x = exp(tage + log(adult.age + 1)) - 1
  }else{
    x = tage * (adult.age + 1) + adult.age
  }
  return(x)
}


##########################################################################
# LOAD COVARIATE FILE
##########################################################################
print("Loading covariate file ...")
df.cov <- read.table(cov.train, header = TRUE, row.names = 1, sep = '\t')
age <- unlist(df.cov[1,], use.names = FALSE)
age.transformed <- sapply(age, transform.age)


##########################################################################
# LOAD METH FILE
##########################################################################

##### TRAINING SAMPLES #####
print("Loading meth file ...")
df.meth <- read.table(meth.train, header = TRUE, row.names = 1, sep = '\t', skipNul = TRUE)
df.meth.clean <- df.meth[complete.cases(df.meth), ]
meth.training.data <- t(as.matrix(df.meth.clean))


##########################################################################
# FIND OPTIMAL LAMBDA
##########################################################################
FIRST <- TRUE
for(i in 1:number.cvs){
  print( paste( "Running cross-validation ", i, " of ", number.cvs, sep = "") )
  
  N <- 0
  while( N != 99 ){
    glmnet.Training.CV = cv.glmnet(meth.training.data, age.transformed, nfolds = 10, alpha = alpha, family="gaussian")
    N <- length(glmnet.Training.CV$lambda)
  }
  
  if( FIRST ){
    df.lambda <- data.frame(glmnet.Training.CV$lambda)
    colnames(df.lambda) <- paste("CV", i, sep = '')
    
    df.mse <- data.frame(glmnet.Training.CV$cvm)
    colnames(df.mse) <- paste("CV", i, sep = '')
    FIRST = FALSE
  }else{
    newname <- paste("CV", i, sep = '')
    df.lambda[[ newname ]] <- glmnet.Training.CV$lambda
    df.mse[[ newname ]] <- glmnet.Training.CV$cvm
  }
  
}

print( "length(unique(as.list(df.lambda))) == 1 ?" )
print( length(unique(as.list(df.lambda))) == 1 )

df.lambda$Mean <- rowMeans(df.lambda)
df.mse$Mean    <- rowMeans(df.mse)

print( "df.lambda$Mean == df.lambda$CV1 ?" )
print( df.lambda$Mean == df.lambda$CV1 )

lambda <- df.lambda$Mean
mean.mse <- df.mse$Mean 

lambda.min <- lambda[ which.min(mean.mse) ]
mse.min <- min(mean.mse)

##########################################################################
# SAVE
##########################################################################
save(df.lambda, file = "df.lambda.RData")
save(df.mse, file = "df.mse.RData")
save(lambda.min, file = "lambda.min.RData")
save(mse.min, file = "mse.min.RData")
