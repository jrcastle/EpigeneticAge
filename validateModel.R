setwd("/home/jrca253/EpigeneticAge")
library(glmnet)
library(ggplot2)

cov.vali  = "data/cov_vali_noNA.txt"
meth.vali = "data/meth_vali_noNA.txt"

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
df.cov <- read.table(cov.vali, header = TRUE, row.names = 1, sep = '\t')
age <- unlist(df.cov[1,], use.names = FALSE)
age.transformed <- sapply(age, transform.age)


##########################################################################
# LOAD METH FILE
##########################################################################

##### VALIDATION SAMPLES #####
df.meth.validate <- read.table(meth.vali, header = TRUE, row.names = 1, sep = '\t', skipNul = TRUE)
df.meth.validate.clean <- df.meth.validate[complete.cases(df.meth.validate), ]
meth.validate.data <- t(as.matrix(df.meth.validate.clean))


##########################################################################
# LOAD MODEL
##########################################################################
load("glmnet.Training.RData")
load("lambda.glmnet.Training.RData")


##########################################################################
# DNA METHYLATION AGE PREDICTION
##########################################################################
result <- predict(glmnet.Training, meth.training.data, type="response", s=lambda.glmnet.Training)
result <- sapply(result,transform.age.inverse)


##########################################################################
# VALIDATE 
##########################################################################
residual = age - result

p <- ggplot(data.frame(res = residual, weight = 1), aes(x=res)) +
  geom_histogram(binwidth = 5, color="black", fill="white") +
  labs(x = "Sample Age - Meth Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  

png("residual_hist.png")
p
dev.off()

hist(residual, main = "Prediction Residuals", xlab = "Age - MethAge")
dev.off()

png("MethAgevsSampleAge.png")
plot(age, res, main="Methlyation Age vs Sample Age", xlab="Sample Age ", ylab="Methylation Age ", pch=19) 
abline(lm(res~age), col="red") # regression line (y~x) 
rsq <- summary(lm(res~age))$r.squared
text(23,70, paste("r^2 = ", round(rsq, digits = 4), sep = ""))

