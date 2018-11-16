setwd("/home/jrca253/EpigeneticAge")
library(glmnet)
library(ggplot2)

#cov.train  <- "data/cov_K_train.txt"
cov.train  <- "data/cov_K.txt"
#meth.train <- "data/meth_K_cpgs_in_KNT_imputed_train.txt"
meth.train <- "data/meth_K_WhiteBlackDiffMethCpGs_imputed.txt"

alpha     <- 0.5
adult.age <- 20
one.cv    <- TRUE

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
df.cov <- read.table(cov.train, header = TRUE, row.names = 1, sep = '\t')
age <- unlist(df.cov[1,], use.names = FALSE)
age.transformed <- sapply(age, transform.age)


##########################################################################
# LOAD METH TRAINING SAMPLES
##########################################################################
df.meth <- read.table(meth.train, header = TRUE, row.names = 1, sep = '\t', skipNul = FALSE)
meth.training.data <- t(as.matrix(df.meth))


##########################################################################
# TRAIN MODEL
##########################################################################
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

glmnet.Training <- glmnet(meth.training.data, age.transformed, family = "gaussian", alpha = alpha, nlambda = 100)

save(glmnet.Training, file = "glmnet.Training.RData")
save(lambda.glmnet.Training, file = "lambda.glmnet.Training.RData")


##########################################################################
# DNA METHYLATION AGE PREDICTION
##########################################################################
# Note: the version of R on the HPC does not support png
# do not run the remainder of this code on the HPC
q()

result <- predict(glmnet.Training, meth.training.data, type="response", s=lambda.glmnet.Training)
result <- sapply(result,transform.age.inverse)


##### QUICK CHECKS #####
residual = age - result

p <- ggplot(data.frame(res = residual), aes(x=res)) +
  geom_histogram(binwidth = 1, color="black", fill="white") +
  labs(x = "Sample Age - Meth Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  

png("residual_hist_trainsample.png")
p
dev.off()

png("MethAgevsSampleAge_trainsample.png")
plot(age, result, main="Methlyation Age vs Sample Age", xlab="Sample Age ", ylab="Methylation Age ", pch=19) 
abline(lm(result~age), col="red") # regression line (y~x) 
rsq <- summary(lm(result~age))$r.squared
text(23,70, paste("r^2 = ", round(rsq, digits = 4), sep = ""))

