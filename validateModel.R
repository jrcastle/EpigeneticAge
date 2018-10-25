#setwd("/home/jrca253/EpigeneticAge")
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(glmnet)
library(ggplot2)

rm(list=ls()); gc();

seed      <- "123"
model.dir <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')
meth.vali <- paste("data/meth_K_cpgs_in_KNT_imputed_vali_seed", seed, ".txt", sep = "")
cov.vali  <- paste("data/cov_K_vali_seed", seed, ".txt", sep = "")

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
mean(age)
sd(age)
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
load( paste(model.dir, "lambda.min.100CV.RData", sep = '') )
lambda.glmnet.Training <- lambda.min

##########################################################################
# DNA METHYLATION AGE PREDICTION
##########################################################################
result <- predict(glmnet.Training, meth.validate.data, type="response", s=lambda.glmnet.Training)
result <- sapply(result,transform.age.inverse)


##########################################################################
# VALIDATE 
##########################################################################
residual <- result - age
median.error <- median(residual)
stdev.error <- sd(residual)
mean.error <- mean(residual)

p <- ggplot(data.frame(res = residual), aes(x=res)) +
  geom_histogram(binwidth = 5, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 45)
  ) + 
  labs(x = "Sample Age - Meth Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals") + 
  annotate("text", x = -30, y = 40, label = paste("Median Residual = ", round(median.error, digits = 2), sep = "")) + 
  annotate("text", x = -30, y = 38, label = paste("St. Dev Residual = ", round(stdev.error, digits = 2), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "residual_hist.png", sep = '') )
p
dev.off()

png( paste(model.dir, "MethAgevsSampleAge.png", sep = '') )
plot(age, 
     result, 
     main="Methlyation Age vs Sample Age", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result~age), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result~age))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
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
