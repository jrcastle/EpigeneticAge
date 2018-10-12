setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)

model.dir   <- "cpgs_in_KNT_imputed/"
meth.file.T <- "data/meth_T_cpgs_in_KNT_imputed_ClockCpGs.txt"
cov.file.T  <- "data/cov_T.txt"


###########################################################################################
# AGE TRANSFORMATION FUNCTIONS
###########################################################################################
trafo = function(x, adult.age = 20){
  x = (x + 1) / (1 + adult.age); 
  y = ifelse(x <= 1, log(x), x-1);
  y
}

anti.trafo = function(x, adult.age=20){ 
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age) 
}


###########################################################################################
# LOAD METH/COV DATA AND PREP FOR PREDICTION
###########################################################################################

##### METH T ####
meth.data.T <- read.table(meth.file.T, header = TRUE, sep = '\t')
n.samples.T <- dim(meth.data.T)[[2]] - 1
x <- c(1)
x <- rep(x, n.samples.T)
x <- c("Intercept", x)
df <- data.frame(x)
df <- t(df)
colnames(df) <- colnames(meth.data.T)
meth.data.T <- rbind(meth.data.T, df)
row.names(meth.data.T) <- NULL

##### COV T #####
cov.T <- read.table(cov.file.T, header = TRUE, row.names = 1, sep = '\t')


###########################################################################################
# SPLIT BY CANCER SUBTYPE
###########################################################################################

##### DIVIDE SAMPLES BY TYPE #####
tmp <- data.frame(t(cov.T))
LuminalA.samples <-rownames(tmp[ which(tmp$Cancer.subtypeLuminal.A == 1), ])
LuminalB.samples <-rownames(tmp[ which(tmp$Cancer.subtypeLuminal.B == 1), ])
Triplenegative.samples <-rownames(tmp[ which(tmp$Cancer.subtypeTriplenegative == 1), ])
rm(tmp); gc()

##### METH DATA #####
meth.data.LumA <- meth.data.T[,c('position', LuminalA.samples)]
meth.data.LumB <- meth.data.T[,c('position', LuminalB.samples)]
meth.data.TrpN <- meth.data.T[,c('position', Triplenegative.samples)]

##### SAMPLE AGES #####
sample.ages.LumA <- as.numeric(as.vector(cov.T[,LuminalA.samples]["Age",]))
sample.ages.LumB <- as.numeric(as.vector(cov.T[,LuminalB.samples]["Age",]))
sample.ages.TrpN <- as.numeric(as.vector(cov.T[,Triplenegative.samples]["Age",]))


###########################################################################################
# MODEL COEFFICIENTS 
###########################################################################################
clock.cpg.coef <- read.csv(paste(model.dir, "model_coefficients.csv", sep = ''), stringsAsFactors = FALSE)
clock.cpg.coef <- clock.cpg.coef[ c("model.coefficients.name", "model.coefficients.x") ]
clock.cpg.coef[clock.cpg.coef$model.coefficients.name == "(Intercept)", "model.coefficients.name"] <- "Intercept"

##### LumA #####
meth.data.LumA <- meth.data.LumA[tolower(order(meth.data.LumA$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.LumA$position == clock.cpg.coef$model.coefficients.name

##### LumB #####
meth.data.LumB <- meth.data.LumB[tolower(order(meth.data.LumB$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.LumB$position == clock.cpg.coef$model.coefficients.name

##### TrpN #####
meth.data.TrpN <- meth.data.TrpN[tolower(order(meth.data.TrpN$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.TrpN$position == clock.cpg.coef$model.coefficients.name


###########################################################################################
# PREDICT
###########################################################################################

##### LumA #####
meth.data.LumA$position <- NULL
X <- data.matrix(meth.data.LumA)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.LumA <- t(X) %*% beta
result.LumA <- sapply(result.LumA, anti.trafo)

##### LumB #####
meth.data.LumB$position <- NULL
X <- data.matrix(meth.data.LumB)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.LumB <- t(X) %*% beta
result.LumB <- sapply(result.LumB, anti.trafo)

##### TrpN #####
meth.data.TrpN$position <- NULL
X <- data.matrix(meth.data.TrpN)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.TrpN <- t(X) %*% beta
result.TrpN <- sapply(result.TrpN, anti.trafo)


##########################################################################
# LumA PLOTS 
##########################################################################
residual.LumA <- result.LumA - sample.ages.LumA 
mean.error.LumA <- mean(residual.LumA)
stdev.error.LumA <- sd(residual.LumA)

p <- ggplot(data.frame(res = residual.LumA), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(-5, 92)
  ) + 
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 16)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Luminal A Tumor Tissue") + 
  annotate("text", x = 75, y = 15, label = paste("Mean Residual = ", round(mean.error.LumA, digits = 1), sep = "")) + 
  annotate("text", x = 75, y = 14, label = paste("St. Dev Residual = ", round(stdev.error.LumA, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "residual_hist_LumA.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "MethAgevsSampleAge_LumA.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.LumA, 
     result.LumA, 
     main="Methlyation Age vs Sample Age in Luminal A Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,150), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.LumA~sample.ages.LumA), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.LumA~sample.ages.LumA))$r.squared
text(90,30, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# LumB PLOTS 
##########################################################################
residual.LumB <- result.LumB - sample.ages.LumB 
mean.error.LumB <- mean(residual.LumB)
stdev.error.LumB <- sd(residual.LumB)

p <- ggplot(data.frame(res = residual.LumB), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 5)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Luminal B Tumor Tissue") + 
  annotate("text", x = 80, y = 4.8, label = paste("Mean Residual = ", round(mean.error.LumB, digits = 1), sep = "")) + 
  annotate("text", x = 80, y = 4.5, label = paste("St. Dev Residual = ", round(stdev.error.LumB, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "residual_hist_LumB.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "MethAgevsSampleAge_LumB.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.LumB, 
     result.LumB, 
     main="Methlyation Age vs Sample Age in Luminal B Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,150), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.LumB~sample.ages.LumB), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.LumB~sample.ages.LumB))$r.squared
text(90,30, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# TrpN PLOTS 
##########################################################################
residual.TrpN <- result.TrpN - sample.ages.TrpN 
mean.error.TrpN <- mean(residual.TrpN)
stdev.error.TrpN <- sd(residual.TrpN)

p <- ggplot(data.frame(res = residual.TrpN), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 9)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Triple-Negative Tumor Tissue") + 
  annotate("text", x = 30, y = 8, label = paste("Mean Residual = ", round(mean.error.TrpN, digits = 1), sep = "")) + 
  annotate("text", x = 30, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.TrpN, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "residual_hist_TrpN.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "MethAgevsSampleAge_TrpN.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.TrpN, 
     result.TrpN, 
     main="Methlyation Age vs Sample Age in Triple-Negative Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.TrpN~sample.ages.TrpN), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.TrpN~sample.ages.TrpN))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# MERGED RESIDUAL PLOT 
##########################################################################

tmp.LumA <- data.frame(residual.LumA)
tmp.LumA$type <- "Luminal A"
colnames(tmp.LumA) <- c("res", "ttype")

tmp.LumB <- data.frame(residual.LumB)
tmp.LumB$type <- "Luminal B"
colnames(tmp.LumB) <- c("res", "ttype")

tmp.TrpN <- data.frame(residual.TrpN)
tmp.TrpN$type <- "Triple Negative"
colnames(tmp.TrpN) <- c("res", "ttype")

df.ABN <- rbind(tmp.LumA, tmp.LumB, tmp.TrpN)
rm(tmp.LumA); rm(tmp.LumB); rm(tmp.TrpN); gc();


p <- ggplot(df.ABN, aes(x = res)) +
  geom_histogram(data=subset(df.ABN, ttype == 'Luminal A'),       aes(fill = "Luminal A"),       binwidth = 10, alpha = 1) +
  geom_histogram(data=subset(df.ABN, ttype == 'Luminal B'),       aes(fill = "Luminal B"),       binwidth = 10, alpha = 1) +
  geom_histogram(data=subset(df.ABN, ttype == 'Triple Negative'), aes(fill = "Triple Negative"), binwidth = 10, alpha = 1) +
  scale_fill_manual(
    name = "Cancer Subtype", 
    values = c("deepskyblue4","darkorange3","seagreen4"), 
    labels=c("Luminal A", "Luminal B", "Triple Negative")
  ) +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(-40, 125)
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 20)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals for Cancer Subtypes") + 
  annotate(
    "text", x = 90, y = 15, 
    label = paste("Mean Residual (LumA) = ", round(mean.error.LumA,  digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) + 
  annotate(
    "text", x = 90, y = 14, 
    label = paste("St. Dev Residual (LumA) = ", round(stdev.error.LumA, digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) +
  annotate(
    "text", x = 90, y = 12, 
    label = paste("Mean Residual (LumB) = ", round(mean.error.LumB,  digits = 1), sep = ""), 
    color = "darkorange3"
  ) + 
  annotate(
    "text", x = 90, y = 11, 
    label = paste("St. Dev Residual (LumB) = ", round(stdev.error.LumB, digits = 1), sep = ""), 
    color = "darkorange3"
  ) +
  annotate(
    "text", x = 90, y = 9,  
    label = paste("Mean Residual (TrpN) = ", round(mean.error.TrpN,  digits = 1), sep = ""), 
    color = "seagreen4"
  ) + 
  annotate(
    "text", x = 90, y = 8,
    label = paste("St. Dev Residual (TrpN) = ", round(stdev.error.TrpN, digits = 1), sep = ""), 
    color = "seagreen4"
  ) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(legend.position=c(0.8, 0.9)) +
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "residual_hist_CancerSub.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


wilcox.test(residual.LumA, mu=0, conf.int = TRUE)
wilcox.test(residual.LumB, mu=0, conf.int = TRUE)
wilcox.test(residual.TrpN, mu=0, conf.int = TRUE)
