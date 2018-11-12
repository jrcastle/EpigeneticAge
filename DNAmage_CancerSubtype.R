setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(RColorBrewer)

rm(list=ls()); gc();

seed        <- "123"
model.dir   <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')
meth.file.T <- paste("data/meth_T_cpgs_in_KNT_imputed_ClockCpGs_seed", seed, ".txt", sep = "")
cov.file.T  <- paste("data/cov_T_seed", seed, ".txt", sep = "")
 

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
LuminalA.samples       <-rownames(tmp[ which(tmp$Cancer.subtypeLuminal.A == 1), ])
LuminalB.samples       <-rownames(tmp[ which(tmp$Cancer.subtypeLuminal.B == 1), ])
Triplenegative.samples <-rownames(tmp[ which(tmp$Cancer.subtypeTriplenegative == 1), ])
ERPRposHer2pos.samples <- rownames(tmp[ which(tmp$ERPRposHer2pos == 1), ])
ERPRposHer2neg.samples <- rownames(tmp[ which(tmp$ERPRposHer2neg == 1), ])
ERPRnegHer2pos.samples <- rownames(tmp[ which(tmp$ERPRnegHer2pos == 1), ])
ERPRnegHer2neg.samples <- rownames(tmp[ which(tmp$ERPRnegHer2neg == 1), ])

##### DIVIDE SAMPLES BY TUMOR GRADE #####
grade3.5.samples <- rownames(tmp[ which(tmp$Cancer.grade6to7 == 0 & tmp$Cancer.grade8to9 == 0), ])
grade6.7.samples <- rownames(tmp[ which(tmp$Cancer.grade6to7 == 1), ])
grade8.9.samples <- rownames(tmp[ which(tmp$Cancer.grade8to9 == 1), ])

##### DIVIDE SAMPLES BY CANCER STAGE #####
stageII.samples <- rownames(tmp[ which(tmp$Cancer.stageIII == 0 & tmp$Cancer.stageIV == 0), ])
stageIII.IV.samples <- rownames(tmp[ which(tmp$Cancer.stageIII == 1 | tmp$Cancer.stageIV == 1), ])
#stageIV.samples <- rownames(tmp[ which(tmp$Cancer.stageIV == 1), ])
rm(tmp); gc()

##### METH DATA #####
meth.data.LumA <- meth.data.T[,c('position', LuminalA.samples)]
meth.data.LumB <- meth.data.T[,c('position', LuminalB.samples)]
meth.data.TrpN <- meth.data.T[,c('position', Triplenegative.samples)]
meth.data.Tpp     <- meth.data.T[,c('position', ERPRposHer2pos.samples)]
meth.data.Tpn     <- meth.data.T[,c('position', ERPRposHer2neg.samples)]
meth.data.Tnp     <- meth.data.T[,c('position', ERPRnegHer2pos.samples)]
meth.data.Tnn     <- meth.data.T[,c('position', ERPRnegHer2neg.samples)]
meth.data.g3.5    <- meth.data.T[,c('position', grade3.5.samples)]
meth.data.g6.7    <- meth.data.T[,c('position', grade6.7.samples)]
meth.data.g8.9    <- meth.data.T[,c('position', grade8.9.samples)]
meth.data.sII     <- meth.data.T[,c('position', stageII.samples)]
meth.data.sIII.IV <- meth.data.T[,c('position', stageIII.IV.samples)]
#meth.data.sIV  <- meth.data.T[,c('position', stageIV.samples)]

##### SAMPLE AGES #####
sample.ages.LumA <- as.numeric(as.vector(cov.T[,LuminalA.samples]["Age",]))
sample.ages.LumB <- as.numeric(as.vector(cov.T[,LuminalB.samples]["Age",]))
sample.ages.TrpN <- as.numeric(as.vector(cov.T[,Triplenegative.samples]["Age",]))
sample.ages.Tpp     <- as.numeric(as.vector(cov.T[,ERPRposHer2pos.samples]["Age",]))
sample.ages.Tpn     <- as.numeric(as.vector(cov.T[,ERPRposHer2neg.samples]["Age",]))
sample.ages.Tnp     <- as.numeric(as.vector(cov.T[,ERPRnegHer2pos.samples]["Age",]))
sample.ages.Tnn     <- as.numeric(as.vector(cov.T[,ERPRnegHer2neg.samples]["Age",]))
sample.ages.g3.5    <- as.numeric(as.vector(cov.T[,grade3.5.samples]["Age",]))
sample.ages.g6.7    <- as.numeric(as.vector(cov.T[,grade6.7.samples]["Age",]))
sample.ages.g8.9    <- as.numeric(as.vector(cov.T[,grade8.9.samples]["Age",]))
sample.ages.sII     <- as.numeric(as.vector(cov.T[,stageII.samples]["Age",]))
sample.ages.sIII.IV <- as.numeric(as.vector(cov.T[,stageIII.IV.samples]["Age",]))
#sample.ages.sIV  <- as.numeric(as.vector(cov.T[,stageIV.samples]["Age",]))


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

##### Tpp #####
meth.data.Tpp <- meth.data.Tpp[tolower(order(meth.data.Tpp$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.Tpp$position == clock.cpg.coef$model.coefficients.name

##### Tpn #####
meth.data.Tpn <- meth.data.Tpn[tolower(order(meth.data.Tpn$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.Tpn$position == clock.cpg.coef$model.coefficients.name

##### Tnp #####
meth.data.Tnp <- meth.data.Tnp[tolower(order(meth.data.Tnp$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.Tnp$position == clock.cpg.coef$model.coefficients.name

##### Tnn #####
meth.data.Tnn <- meth.data.Tnn[tolower(order(meth.data.Tnn$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.Tnn$position == clock.cpg.coef$model.coefficients.name

##### g3.5 #####
meth.data.g3.5 <- meth.data.g3.5[tolower(order(meth.data.g3.5$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.g3.5$position == clock.cpg.coef$model.coefficients.name

##### g6.7 #####
meth.data.g6.7 <- meth.data.g6.7[tolower(order(meth.data.g6.7$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.g6.7$position == clock.cpg.coef$model.coefficients.name

##### g8.9 #####
meth.data.g8.9 <- meth.data.g8.9[tolower(order(meth.data.g8.9$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.g8.9$position == clock.cpg.coef$model.coefficients.name

##### sII #####
meth.data.sII <- meth.data.sII[tolower(order(meth.data.sII$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.sII$position == clock.cpg.coef$model.coefficients.name

##### sIII.IV #####
meth.data.sIII.IV <- meth.data.sIII.IV[tolower(order(meth.data.sIII.IV$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.sIII.IV$position == clock.cpg.coef$model.coefficients.name

##### sIV #####
#meth.data.sIV <- meth.data.sIV[tolower(order(meth.data.sIV$position)),]
#clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
#meth.data.sIV$position == clock.cpg.coef$model.coefficients.name


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

##### Tpp #####
meth.data.Tpp$position <- NULL
X <- data.matrix(meth.data.Tpp)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.Tpp <- t(X) %*% beta
result.Tpp <- sapply(result.Tpp, anti.trafo)

##### Tpn #####
meth.data.Tpn$position <- NULL
X <- data.matrix(meth.data.Tpn)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.Tpn <- t(X) %*% beta
result.Tpn <- sapply(result.Tpn, anti.trafo)

##### Tnp #####
meth.data.Tnp$position <- NULL
X <- data.matrix(meth.data.Tnp)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.Tnp <- t(X) %*% beta
result.Tnp <- sapply(result.Tnp, anti.trafo)

##### Tnn #####
meth.data.Tnn$position <- NULL
X <- data.matrix(meth.data.Tnn)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.Tnn <- t(X) %*% beta
result.Tnn <- sapply(result.Tnn, anti.trafo)

##### g3.5 #####
meth.data.g3.5$position <- NULL
X <- data.matrix(meth.data.g3.5)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.g3.5 <- t(X) %*% beta
result.g3.5 <- sapply(result.g3.5, anti.trafo)

##### g6.7 #####
meth.data.g6.7$position <- NULL
X <- data.matrix(meth.data.g6.7)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.g6.7 <- t(X) %*% beta
result.g6.7 <- sapply(result.g6.7, anti.trafo)

##### g8.9 #####
meth.data.g8.9$position <- NULL
X <- data.matrix(meth.data.g8.9)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.g8.9 <- t(X) %*% beta
result.g8.9 <- sapply(result.g8.9, anti.trafo)

##### sII #####
meth.data.sII$position <- NULL
X <- data.matrix(meth.data.sII)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.sII <- t(X) %*% beta
result.sII <- sapply(result.sII, anti.trafo)

##### sIII.IV #####
meth.data.sIII.IV$position <- NULL
X <- data.matrix(meth.data.sIII.IV)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.sIII.IV <- t(X) %*% beta
result.sIII.IV <- sapply(result.sIII.IV, anti.trafo)

##### sIV #####
#meth.data.sIV$position <- NULL
#X <- data.matrix(meth.data.sIV)
#beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

#result.sIV <- t(X) %*% beta
#result.sIV <- sapply(result.sIV, anti.trafo)


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

png( paste(model.dir, "CancerStudies/residual_hist_LumA.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_LumA.png", sep = ''), width = 500, height = 500, units = "px" )
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


png( paste(model.dir, "CancerStudies/residual_hist_LumB.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_LumB.png", sep = ''), width = 500, height = 500, units = "px" )
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


png( paste(model.dir, "CancerStudies/residual_hist_TrpN.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_TrpN.png", sep = ''), width = 500, height = 500, units = "px" )
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
# Tpp PLOTS 
##########################################################################
residual.Tpp <- result.Tpp - sample.ages.Tpp 
mean.error.Tpp <- mean(residual.Tpp)
stdev.error.Tpp <- sd(residual.Tpp)

p <- ggplot(data.frame(res = residual.Tpp), aes(x=res)) +
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
  labs(title = "Prediction Residuals in ER/PR+ Her2+ Tumor Tissue") + 
  annotate("text", x = 75, y = 15, label = paste("Mean Residual = ", round(mean.error.Tpp, digits = 1), sep = "")) + 
  annotate("text", x = 75, y = 14, label = paste("St. Dev Residual = ", round(stdev.error.Tpp, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_Tpp.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_Tpp.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.Tpp, 
     result.Tpp, 
     main="Methlyation Age vs Sample Age in ER/PR+ Her2+ Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,150), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.Tpp~sample.ages.Tpp), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.Tpp~sample.ages.Tpp))$r.squared
text(90,30, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# Tpn PLOTS 
##########################################################################
residual.Tpn <- result.Tpn - sample.ages.Tpn 
mean.error.Tpn <- mean(residual.Tpn)
stdev.error.Tpn <- sd(residual.Tpn)

p <- ggplot(data.frame(res = residual.Tpn), aes(x=res)) +
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
  labs(title = "Prediction Residuals in ER/PR+ Her2- Tumor Tissue") + 
  annotate("text", x = 75, y = 15, label = paste("Mean Residual = ", round(mean.error.Tpn, digits = 1), sep = "")) + 
  annotate("text", x = 75, y = 14, label = paste("St. Dev Residual = ", round(stdev.error.Tpn, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_Tpn.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_Tpn.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.Tpn, 
     result.Tpn, 
     main="Methlyation Age vs Sample Age in ER/PR+ Her2- Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,150), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.Tpn~sample.ages.Tpn), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.Tpn~sample.ages.Tpn))$r.squared
text(90,30, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# Tnp PLOTS 
##########################################################################
residual.Tnp <- result.Tnp - sample.ages.Tnp 
mean.error.Tnp <- mean(residual.Tnp)
stdev.error.Tnp <- sd(residual.Tnp)

p <- ggplot(data.frame(res = residual.Tnp), aes(x=res)) +
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
  labs(title = "Prediction Residuals in ER/PR- Her2+ Tumor Tissue") + 
  annotate("text", x = 75, y = 15, label = paste("Mean Residual = ", round(mean.error.Tnp, digits = 1), sep = "")) + 
  annotate("text", x = 75, y = 14, label = paste("St. Dev Residual = ", round(stdev.error.Tnp, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_Tnp.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_Tnp.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.Tnp, 
     result.Tnp, 
     main="Methlyation Age vs Sample Age in ER/PR- Her2+ Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,150), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.Tnp~sample.ages.Tnp), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.Tnp~sample.ages.Tnp))$r.squared
text(90,30, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# Tnn PLOTS 
##########################################################################
residual.Tnn <- result.Tnn - sample.ages.Tnn 
mean.error.Tnn <- mean(residual.Tnn)
stdev.error.Tnn <- sd(residual.Tnn)

p <- ggplot(data.frame(res = residual.Tnn), aes(x=res)) +
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
  labs(title = "Prediction Residuals in ER/PR- Her2- Tumor Tissue") + 
  annotate("text", x = 75, y = 15, label = paste("Mean Residual = ", round(mean.error.Tnn, digits = 1), sep = "")) + 
  annotate("text", x = 75, y = 14, label = paste("St. Dev Residual = ", round(stdev.error.Tnn, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_Tnn.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_Tnn.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.Tnn, 
     result.Tnn, 
     main="Methlyation Age vs Sample Age in ER/PR- Her2- Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,150), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.Tnn~sample.ages.Tnn), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.Tnn~sample.ages.Tnn))$r.squared
text(90,30, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# g3.5 PLOTS 
##########################################################################
residual.g3.5 <- result.g3.5 - sample.ages.g3.5 
mean.error.g3.5 <- mean(residual.g3.5)
stdev.error.g3.5 <- sd(residual.g3.5)

p <- ggplot(data.frame(res = residual.g3.5), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 9)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Grade 3-5 Tumor Tissue") + 
  annotate("text", x = 15, y = 8, label = paste("Mean Residual = ", round(mean.error.g3.5, digits = 1), sep = "")) + 
  annotate("text", x = 15, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.g3.5, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_g3.5.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_g3.5.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.g3.5, 
     result.g3.5, 
     main="Methlyation Age vs Sample Age in Grade 3-5 Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.g3.5~sample.ages.g3.5), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.g3.5~sample.ages.g3.5))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# g6.7 PLOTS 
##########################################################################
residual.g6.7 <- result.g6.7 - sample.ages.g6.7 
mean.error.g6.7 <- mean(residual.g6.7)
stdev.error.g6.7 <- sd(residual.g6.7)

p <- ggplot(data.frame(res = residual.g6.7), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 9)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Grade 6-7 Tumor Tissue") + 
  annotate("text", x = 80, y = 8, label = paste("Mean Residual = ", round(mean.error.g6.7, digits = 1), sep = "")) + 
  annotate("text", x = 80, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.g6.7, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_g6.7.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_g6.7.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.g6.7, 
     result.g6.7, 
     main="Methlyation Age vs Sample Age in Grade 6-7 Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.g6.7~sample.ages.g6.7), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.g6.7~sample.ages.g6.7))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# g8.9 PLOTS 
##########################################################################
residual.g8.9 <- result.g8.9 - sample.ages.g8.9 
mean.error.g8.9 <- mean(residual.g8.9)
stdev.error.g8.9 <- sd(residual.g8.9)

p <- ggplot(data.frame(res = residual.g8.9), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 15)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Grade 8-9 Tumor Tissue") + 
  annotate("text", x = 75, y = 8, label = paste("Mean Residual = ", round(mean.error.g8.9, digits = 1), sep = "")) + 
  annotate("text", x = 75, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.g8.9, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_g8.9.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_g8.9.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.g8.9, 
     result.g8.9, 
     main="Methlyation Age vs Sample Age in Grade 8-9 Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.g8.9~sample.ages.g8.9), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.g8.9~sample.ages.g8.9))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# sII PLOTS 
##########################################################################
residual.sII <- result.sII - sample.ages.sII 
mean.error.sII <- mean(residual.sII)
stdev.error.sII <- sd(residual.sII)

p <- ggplot(data.frame(res = residual.sII), aes(x=res)) +
  geom_histogram(binwidth = 15, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 9)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Stage II Tumor Tissue") + 
  annotate("text", x = 65, y = 8, label = paste("Mean Residual = ", round(mean.error.sII, digits = 1), sep = "")) + 
  annotate("text", x = 65, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.sII, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_sII.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_sII.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.sII, 
     result.sII, 
     main="Methlyation Age vs Sample Age in Stage II Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.sII~sample.ages.sII), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.sII~sample.ages.sII))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# sIII.IV PLOTS 
##########################################################################
residual.sIII.IV <- result.sIII.IV - sample.ages.sIII.IV 
mean.error.sIII.IV <- mean(residual.sIII.IV)
stdev.error.sIII.IV <- sd(residual.sIII.IV)

p <- ggplot(data.frame(res = residual.sIII.IV), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 9)
  ) + 
  labs(x = "Meth Age - Sample Age") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Stages III & IV Tumor Tissue") + 
  annotate("text", x = -25, y = 8, label = paste("Mean Residual = ", round(mean.error.sIII.IV, digits = 1), sep = "")) + 
  annotate("text", x = -25, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.sIII.IV, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_sIII.IV.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_sIII.IV.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.sIII.IV, 
     result.sIII.IV, 
     main="Methlyation Age vs Sample Age in Stages III & IV Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.sIII.IV~sample.ages.sIII.IV), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.sIII.IV~sample.ages.sIII.IV))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# sIV PLOTS 
##########################################################################
#residual.sIV <- result.sIV - sample.ages.sIV 
#mean.error.sIV <- mean(residual.sIV)
#stdev.error.sIV <- sd(residual.sIV)

#p <- ggplot(data.frame(res = residual.sIV), aes(x=res)) +
#  geom_histogram(binwidth = 10, color="black", fill="white") +
#  scale_y_continuous(
#    expand=c(0, 0),
#    limits = c(0, 9)
#  ) + 
#  labs(x = "Meth Age - Sample Age") + 
#  labs(y = "Frequency") + 
#  labs(title = "Prediction Residuals in Stage IV Tumor Tissue") + 
#  annotate("text", x = 18, y = 8, label = paste("Mean Residual = ", round(mean.error.sIV, digits = 1), sep = "")) + 
#  annotate("text", x = 18, y = 7.5, label = paste("St. Dev Residual = ", round(stdev.error.sIV, digits = 1), sep = "")) +
#  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
#  theme_bw() + 
#  theme(
#    axis.ticks.length=unit(-0.25, "cm"), 
#    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
#    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
#    plot.title = element_text(hjust = 0.5)
#  )


#png( paste(model.dir, "CancerStudies/residual_hist_sIV.png", sep = ''), width = 500, height = 500, units = "px" )
#p
#dev.off()

#png( paste(model.dir, "CancerStudies/MethAgevsSampleAge_sIV.png", sep = ''), width = 500, height = 500, units = "px" )
#plot(sample.ages.sIV, 
#     result.sIV, 
#     main="Methlyation Age vs Sample Age in Stage IV Tumor Tissue", 
#     xlab="Sample Age ", 
#     ylab="Methylation Age ", 
#     pch=19,
#     xlim=c(15,87), 
#     xaxs="i",
#     ylim=c(15,87), 
#     yaxs="i",
#     tck = 0.02
#) 
#abline(lm(result.sIV~sample.ages.sIV), col="black") # regression line (y~x) 
#abline(a=0, b=1, col = "red", lty = 2)
#rsq <- summary(lm(result.sIV~sample.ages.sIV))$r.squared
#text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
#dev.off()


##########################################################################
# MERGED SUBTYPE RESIDUAL PLOT 
##########################################################################
median.error.LumA = median(residual.LumA)
median.error.LumB = median(residual.LumB)
median.error.TrpN = median(residual.TrpN)

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


p <- ggplot(df.ABN, aes(x = res, stat(density), color = ttype, linetype = ttype)) +
  geom_freqpoly(size = 1.2, binwidth = 20) +
  scale_linetype_manual(
    name = "Cancer Subtype", 
    values = c("solid", "dashed", "twodash")
  ) + 
  scale_color_manual(
    name = "Cancer Subtype", 
    values = c("deepskyblue4","darkorange3","seagreen4"), 
    breaks = c("Luminal A", "Luminal B", "Triple Negative")
  ) +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(-40, 125)
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 0.05)
  ) + 
  labs(x = "DNAm Age Acceleration") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals for Cancer Subtypes") + 
  annotate(
    "text", x = 90, y = 0.037, 
    label = paste("Median Residual (LumA) = ", round(median.error.LumA,  digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) + 
  annotate(
    "text", x = 90, y = 0.035, 
    label = paste("St. Dev Residual (LumA) = ", round(stdev.error.LumA, digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) +
  annotate(
    "text", x = 90, y = 0.032, 
    label = paste("Median Residual (LumB) = ", round(median.error.LumB,  digits = 1), sep = ""), 
    color = "darkorange3"
  ) + 
  annotate(
    "text", x = 90, y = 0.030, 
    label = paste("St. Dev Residual (LumB) = ", round(stdev.error.LumB, digits = 1), sep = ""), 
    color = "darkorange3"
  ) +
  annotate(
    "text", x = 90, y = 0.027,  
    label = paste("Median Residual (TrpN) = ", round(median.error.TrpN,  digits = 1), sep = ""), 
    color = "seagreen4"
  ) + 
  annotate(
    "text", x = 90, y = 0.025,
    label = paste("St. Dev Residual (TrpN) = ", round(stdev.error.TrpN, digits = 1), sep = ""), 
    color = "seagreen4"
  ) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(legend.key.width = unit(3, "line"), legend.position=c(0.8, 0.9)) +
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/residual_hist_CancerSub.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED SUBTYPE RESIDUAL PLOT 2
##########################################################################
median.error.Tpp = median(residual.Tpp)
median.error.Tpn = median(residual.Tpn)
median.error.Tnp = median(residual.Tnp)
median.error.Tnn = median(residual.Tnn)

tmp.Tpp <- data.frame(residual.Tpp)
tmp.Tpp$type <- "1. ER/PR+ Her2+"
colnames(tmp.Tpp) <- c("res", "ttype")

tmp.Tpn <- data.frame(residual.Tpn)
tmp.Tpn$type <- "2. ER/PR+ Her2-"
colnames(tmp.Tpn) <- c("res", "ttype")

tmp.Tnp <- data.frame(residual.Tnp)
tmp.Tnp$type <- "3. ER/PR- Her2+"
colnames(tmp.Tnp) <- c("res", "ttype")

tmp.Tnn <- data.frame(residual.Tnn)
tmp.Tnn$type <- "4. ER/PR- Her2-"
colnames(tmp.Tnn) <- c("res", "ttype")

df.ABN <- rbind(tmp.Tpp, tmp.Tpn, tmp.Tnp, tmp.Tnn)
rm(tmp.Tpp); rm(tmp.Tpn); rm(tmp.Tnp); rm(tmp.Tnn); gc();

p <- ggplot(df.ABN, aes(x = res, stat(density), color = ttype, linetype = ttype)) +
  geom_freqpoly(size = 1.2, binwidth = 25) +
  scale_linetype_manual(
    name = "Cancer Subtype", 
    values = c("solid", "dashed", "twodash", "longdash"),
    breaks =  c("1. ER/PR+ Her2+", "2. ER/PR+ Her2-", "3. ER/PR- Her2+", "4. ER/PR- Her2-")
  ) + 
  scale_color_manual(
    name = "Cancer Subtype", 
    values = c("deepskyblue4","darkorange3","seagreen4", "firebrick3"), 
    breaks = c("1. ER/PR+ Her2+", "2. ER/PR+ Her2-", "3. ER/PR- Her2+", "4. ER/PR- Her2-")
  ) +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(-40, 125)
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 0.035)
  ) + 
  labs(x = "DNAm Age Acceleration [Years]") + 
  labs(y = "Frequency [Arbitrary Units]") + 
  labs(title = "DNAm Age Acceleration for Cancer Subtypes") + 
  annotate(
    "text", x = 85, y = 0.025, 
    label = paste("Median Accel. (ER/PR+ Her2+) = ", round(median.error.Tpp,  digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) + 
  #annotate(
  #  "text", x = 85, y = 0.0245, 
  #  label = paste("St. Dev Accel. (ER/PR+ Her2+) = ", round(stdev.error.Tpp, digits = 1), sep = ""), 
  #  color = "deepskyblue4"
  #) +
  annotate(
    "text", x = 85, y = 0.023, 
    label = paste("Median Accel. (ER/PR+ Her2-) = ", round(median.error.Tpn,  digits = 1), sep = ""), 
    color = "darkorange3"
  ) + 
  #annotate(
  #  "text", x = 85, y = 0.021, 
  #  label = paste("St. Dev Accel. (ER/PR+ Her2-) = ", round(stdev.error.Tpn, digits = 1), sep = ""), 
  #  color = "darkorange3"
  #) +
  annotate(
    "text", x = 85, y = 0.021,  
    label = paste("Median Accel. (ER/PR- Her2+) = ", round(median.error.Tnp,  digits = 1), sep = ""), 
    color = "seagreen4"
  ) + 
  #annotate(
  #  "text", x = 85, y = 0.0175,
  #  label = paste("St. Dev Accel. (ER/PR- Her2+) = ", round(stdev.error.Tnp, digits = 1), sep = ""), 
  #  color = "seagreen4"
  #) +
  annotate(
    "text", x = 85, y = 0.019,  
    label = paste("Median Accel. (ER/PR- Her2-) = ", round(median.error.Tnn,  digits = 1), sep = ""), 
    color = "firebrick3"
  ) + 
  #annotate(
  #  "text", x = 85, y = 0.014,
  #  label = paste("St. Dev Accel. (ER/PR- Her2-) = ", round(stdev.error.Tnn, digits = 1), sep = ""), 
  #  color = "firebrick3"
  #) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
  theme_bw(base_size = 15) +
  theme(legend.key.width = unit(3, "line"), legend.position=c(0.8, 0.87)) +
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "CancerStudies/1_residual_hist_CancerSub.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED GRADE PLOT
##########################################################################
median.error.g3.5 <- median(residual.g3.5)
median.error.g6.7 <- median(residual.g6.7)
median.error.g8.9 <- median(residual.g8.9)

tmp.g3.5 <- data.frame(residual.g3.5)
tmp.g3.5$type <- "Grade 3-5"
colnames(tmp.g3.5) <- c("res", "ttype")

tmp.g6.7 <- data.frame(residual.g6.7)
tmp.g6.7$type <- "Grade 6-7"
colnames(tmp.g6.7) <- c("res", "ttype")

tmp.g8.9 <- data.frame(residual.g8.9)
tmp.g8.9$type <- "Grade 8-9"
colnames(tmp.g8.9) <- c("res", "ttype")

df.ABN <- rbind(tmp.g3.5, tmp.g6.7, tmp.g8.9)
rm(tmp.g3.5); rm(tmp.g6.7); rm(tmp.g8.9); gc();

p <- ggplot(df.ABN, aes(x = res, stat(density), color = ttype, linetype = ttype)) +
  #geom_histogram(data=subset(df.ABN, ttype == 'Grade 3-5'), aes(fill = "Grade 3-5"), binwidth = 10, alpha = 1) +
  #geom_histogram(data=subset(df.ABN, ttype == 'Grade 6-7'), aes(fill = "Grade 6-7"), binwidth = 10, alpha = 1) +
  #geom_histogram(data=subset(df.ABN, ttype == 'Grade 8-9'), aes(fill = "Grade 8-9"), binwidth = 10, alpha = 1) +
  geom_freqpoly(size = 1.2, binwidth = 20) +
  scale_linetype_manual(
    name = "Tumor Grade", 
    values = c("solid", "dashed", "twodash"),
    labels = c("Grade 3-5", "Grade 6-7", "Grade 8-9")
  ) + 
  scale_color_manual(
    name = "Tumor Grade", 
    values = c("gray16","pink4","paleturquoise4"), 
    labels = c("Grade 3-5", "Grade 6-7", "Grade 8-9")
  ) +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(-40, 125)
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 0.045)
  ) + 
  labs(x = "DNAm Age Acceleration [Years]") + 
  labs(y = "Frequency [Arbitrary Units]") + 
  labs(title = "DNAm Age Acceleration for Tumor Grades") + 
  annotate(
    "text", x = 90, y = 0.032, 
    label = paste("Median Accel. (3-5) = ", round(median.error.g3.5,  digits = 1), sep = ""), 
    color = "gray16"
  ) + 
  #annotate(
  #  "text", x = 90, y = 0.030, 
  #  label = paste("St. Dev Accel. (3-5) = ", round(stdev.error.g3.5, digits = 1), sep = ""), 
  #  color = "gray16"
  #) +
  annotate(
    "text", x = 90, y = 0.030, 
    label = paste("Median Accel. (6-7) = ", round(median.error.g6.7,  digits = 1), sep = ""), 
    color = "pink4"
  ) + 
  #annotate(
  #  "text", x = 90, y = 0.025, 
  #  label = paste("St. Dev Accel. (6-7) = ", round(stdev.error.g6.7, digits = 1), sep = ""), 
  #  color = "pink4"
  #) +
  annotate(
    "text", x = 90, y = 0.028,  
    label = paste("Median Accel. (8-9) = ", round(median.error.g8.9,  digits = 1), sep = ""), 
    color = "paleturquoise4"
  ) + 
  #annotate(
  #  "text", x = 90, y = 0.020,
  #  label = paste("St. Dev Accel. (8-9) = ", round(stdev.error.g8.9, digits = 1), sep = ""), 
  #  color = "paleturquoise4"
  #) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
  theme_bw(base_size = 15) +
  theme(legend.key.width = unit(3, "line"), legend.position=c(0.8, 0.87)) +
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  ); p

png( paste(model.dir, "CancerStudies/2_residual_hist_TumorGrade.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED STAGE PLOT
##########################################################################
median.error.sII <- median(residual.sII)
median.error.sIII.IV <- median(residual.sIII.IV)
#median.error.sIV <- median(residual.sIV)

tmp.sII <- data.frame(residual.sII)
tmp.sII$type <- "Stage II"
colnames(tmp.sII) <- c("res", "ttype")

tmp.sIII.IV <- data.frame(residual.sIII.IV)
tmp.sIII.IV$type <- "Stage III & IV"
colnames(tmp.sIII.IV) <- c("res", "ttype")

#tmp.sIV <- data.frame(residual.sIV)
#tmp.sIV$type <- "Stage IV"
#colnames(tmp.sIV) <- c("res", "ttype")

df.ABN <- rbind(tmp.sII, tmp.sIII.IV)
rm(tmp.sII); rm(tmp.sIII.IV); gc();

p <- ggplot(df.ABN, aes(x = res, stat(density), color = ttype, linetype = ttype)) +
  #geom_histogram(data=subset(df.ABN, ttype == 'Stage II'),  aes(fill = "Stage II"),  binwidth = 20, alpha = 1) +
  #geom_histogram(data=subset(df.ABN, ttype == 'Stage III & IV'), aes(fill = "Stage III & IV"), binwidth = 20, alpha = 1) +
  #geom_histogram(data=subset(df.ABN, ttype == 'Stage IV'),  aes(fill = "Stage IV"),  binwidth = 20, alpha = 1) +
  geom_freqpoly(size = 1.2, binwidth = 20) +
  scale_linetype_manual(
    name = "Tumor Stage", 
    values = c("solid", "dashed"),
    labels = c("Stage II", "Stage III & IV")
  ) +
  scale_color_manual(
    name = "Tumor Stage", 
    values = c("darkslategrey","firebrick3"), 
    labels = c("Stage II", "Stage III & IV")
  ) +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(-40, 125)
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 0.04)
  ) + 
  labs(x = "DNAm Age Acceleration [Years]") + 
  labs(y = "Frequency [Arbitrary Units]") + 
  labs(title = "DNAm Age Acceleration for Tumor Stages") + 
  annotate(
    "text", x = 90, y = 0.030, 
    label = paste("Median Accel. (II) = ", round(median.error.sII,  digits = 1), sep = ""), 
    color = "darkslategrey"
  ) + 
  #annotate(
  #  "text", x = 90, y = 0.028, 
  #  label = paste("St. Dev Accel. (II) = ", round(stdev.error.sII, digits = 1), sep = ""), 
  #  color = "darkslategrey"
  #) +
  annotate(
    "text", x = 90, y = 0.028, 
    label = paste("Median Accel. (III & IV) = ", round(median.error.sIII.IV,  digits = 1), sep = ""), 
    color = "firebrick3"
  ) + 
  #annotate(
  #  "text", x = 90, y = 0.023, 
  #  label = paste("St. Dev Accel. (III & IV) = ", round(stdev.error.sIII.IV, digits = 1), sep = ""), 
  #  color = "firebrick3"
  #) +
  #annotate(
  #  "text", x = 90, y = 0.051,  
  #  label = paste("Median Accel. (IV) = ", round(median.error.sIV,  digits = 1), sep = ""), 
  #  color = "paleturquoise4"
  #) + 
  #annotate(
  #  "text", x = 90, y = 0.047,
  #  label = paste("St. Dev Accel. (IV) = ", round(stdev.error.sIV, digits = 1), sep = ""), 
  #  color = "paleturquoise4"
  #) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
  theme_bw(base_size = 15) + 
  theme(legend.key.width = unit(3, "line"), legend.position=c(0.8, 0.87)) +
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  ); p

png( paste(model.dir, "CancerStudies/3_residual_hist_TumorStage.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# STATISTICAL TESTS
##########################################################################
wilcox.test(residual.LumA, mu=0, conf.int = TRUE)
wilcox.test(residual.LumB, mu=0, conf.int = TRUE)
wilcox.test(residual.TrpN, mu=0, conf.int = TRUE)
wilcox.test(residual.Tpp, mu=0, conf.int = TRUE)
wilcox.test(residual.Tpn, mu=0, conf.int = TRUE)
wilcox.test(residual.Tnp, mu=0, conf.int = TRUE)
wilcox.test(residual.Tnn, mu=0, conf.int = TRUE)
wilcox.test(residual.g3.5, mu=0, conf.int = TRUE)
wilcox.test(residual.g6.7, mu=0, conf.int = TRUE)
wilcox.test(residual.g8.9, mu=0, conf.int = TRUE)
wilcox.test(residual.sII,  mu=0, conf.int = TRUE)
wilcox.test(residual.sIII.IV, mu=0, conf.int = TRUE)
#wilcox.test(residual.sIV,  mu=0, conf.int = TRUE)

