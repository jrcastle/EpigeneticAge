setwd("/home/jrca253/EpigeneticAge")
library(glmnet)
library(ggplot2)

cov.train  = "data/cov_train_noNA.txt"
meth.train = "data/meth_train_noNA.txt"

alpha = 0.5
adult.age = 20

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
# LOAD METH FILE
##########################################################################

##### TRAINING SAMPLES #####
df.meth <- read.table(meth.train, header = TRUE, row.names = 1, sep = '\t', skipNul = TRUE)
df.meth.clean <- df.meth[complete.cases(df.meth), ]
meth.training.data <- t(as.matrix(df.meth.clean))


##########################################################################
# TRAIN MODEL
##########################################################################
glmnet.Training.CV = cv.glmnet(meth.training.data, age.transformed, nfolds = 10, alpha = alpha, family="gaussian")
lambda.glmnet.Training = glmnet.Training.CV$lambda.min
glmnet.Training = glmnet(meth.training.data, age.transformed, family = "gaussian", alpha = alpha, nlambda = 100)


##########################################################################
# DNA METHYLATION AGE PREDICTION
##########################################################################
result <- predict(glmnet.Training, meth.training.data, type="response", s=lambda.glmnet.Training)
result <- sapply(result,transform.age.inverse)
save(glmnet.Training, file = "glmnet.Training.RData")
save(lambda.glmnet.Training, file = "lambda.glmnet.Training.RData")

##### QUICK CHECKS #####
residual = age - result

p <- ggplot(data.frame(res = residual, weight = 1), aes(x=res)) +
  geom_histogram(binwidth = 5, color="black", fill="white") +
  labs(x = "Sample Age - Meth Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  

png("residual_hist_trainsample.png")
p
dev.off()

hist(residual, main = "Prediction Residuals", xlab = "Age - MethAge")
dev.off()

png("MethAgevsSampleAge_trainsample.png")
plot(age, res, main="Methlyation Age vs Sample Age", xlab="Sample Age ", ylab="Methylation Age ", pch=19) 
abline(lm(res~age), col="red") # regression line (y~x) 
rsq <- summary(lm(res~age))$r.squared
text(23,70, paste("r^2 = ", round(rsq, digits = 4), sep = ""))

