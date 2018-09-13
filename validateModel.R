#setwd("/home/jrca253/EpigeneticAge")
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(glmnet)
library(ggplot2)

cov.vali  = "data/cov_vali_noNA_split4.txt"
meth.vali = "data/meth_vali_noNA_split4.txt"
model.dir = "betas_only_split4/"
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
df.meth.validate <- read.table(meth.vali, header = TRUE, row.names = 1, sep = '\t')
#df.meth.validate.clean <- df.meth.validate[complete.cases(df.meth.validate), ]
#meth.validate.data <- t(as.matrix(df.meth.validate.clean))
meth.validate.data <- t(as.matrix(df.meth.validate))


##########################################################################
# LOAD MODEL
##########################################################################
load( paste(model.dir, "glmnet.Training.RData", sep = '') )
load( paste(model.dir, "lambda.glmnet.Training.RData", sep = '') )


##########################################################################
# DNA METHYLATION AGE PREDICTION
##########################################################################
result <- predict(glmnet.Training, meth.validate.data, type="response", s=lambda.glmnet.Training)
result <- sapply(result,transform.age.inverse)


##########################################################################
# VALIDATE 
##########################################################################
residual = age - result
median(residual)
p <- ggplot(data.frame(res = residual), aes(x=res)) +
  geom_histogram(binwidth = 5, color="black", fill="white") +
  labs(x = "Sample Age - Meth Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
  

png( paste(model.dir, "residual_hist.png", sep = '') )
p
dev.off()

png( paste(model.dir, "MethAgevsSampleAge.png", sep = '') )
plot(age, result, main="Methlyation Age vs Sample Age", xlab="Sample Age ", ylab="Methylation Age ", pch=19) 
abline(lm(result~age), col="red") # regression line (y~x) 
rsq <- summary(lm(result~age))$r.squared
text(24,70, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# MODEL COEFFICIENTS
##########################################################################
fit.coefficients <- coef(glmnet.Training, s = lambda.glmnet.Training)
fit.features <- rownames(fit.coefficients)

model.coefficients <- summary(fit.coefficients)
indices <- model.coefficients[,1]

model.coefficients$name = fit.features[indices]
model.coefficients <- data.frame(model.coefficients$name, model.coefficients$x)
write.table(model.coefficients, paste(model.dir, "model_coefficients.csv", sep = ''), sep = ",", row.names=FALSE)



