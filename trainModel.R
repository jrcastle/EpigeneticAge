setwd("/home/jrca253/EpigeneticAge")
library(glmnet)
library(ggplot2)

cov.train  <- "data/cov_K_train.txt"
meth.train <- "data/meth_K_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_train.txt"

alpha     <- 0.5
adult.age <- 20
one.cv    <- FALSE

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
print("Loading cov file ...")
df.cov <- read.table(cov.train, header = TRUE, row.names = 1, sep = '\t')
age <- unlist(df.cov[1,], use.names = FALSE)
age.transformed <- sapply(age, transform.age)


##########################################################################
# LOAD METH TRAINING SAMPLES
##########################################################################
print("Loading meth file ...")
df.meth <- read.table(meth.train, header = TRUE, row.names = 1, sep = '\t', skipNul = FALSE)
meth.training.data <- t(as.matrix(df.meth))


##########################################################################
# TRAIN MODEL
##########################################################################
print("Training ...")
if( one.cv ){
  # This is old code that does a single cross validation, the glmnet
  # package notes that the results of this are random and will fluctuate
  # from CV to CV. The suggested way to reduce the randomness is to run 
  # multiple CVs and average the results. This is accomplished in 
  # multipleCV.R
  glmnet.Training.CV <- cv.glmnet(meth.training.data, age.transformed, nfolds = 10, alpha = alpha, family="gaussian")
  lambda.glmnet.Training <- glmnet.Training.CV$lambda.min
}else{
  load( "lambda.min.100CV.RData" )
  lambda.glmnet.Training <- lambda.min
}

print("lambda.glmnet.Training:")
print(lambda.glmnet.Training)

glmnet.Training <- glmnet(meth.training.data, age.transformed, family = "gaussian", alpha = alpha, nlambda = 100)

save(glmnet.Training, file = "glmnet.Training.RData")
save(lambda.glmnet.Training, file = "lambda.glmnet.Training.RData")


##########################################################################
# MODEL COEFFICIENTS
##########################################################################
fit.coefficients <- coef(glmnet.Training, s = lambda.glmnet.Training)
fit.features <- rownames(fit.coefficients)

model.coefficients <- summary(fit.coefficients)
indices <- model.coefficients[,1]

model.coefficients$name = fit.features[indices]
model.coefficients <- data.frame(model.coefficients$name, model.coefficients$x)
write.table(model.coefficients, "model_coefficients.csv", sep = ",", row.names=FALSE)