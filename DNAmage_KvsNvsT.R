setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(RColorBrewer)

rm(list=ls()); gc();

seed        <- "123"
model.dir   <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')
meth.file.K <- paste("data/meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs_seed", seed, ".txt", sep = "")
cov.file.K  <- paste("data/cov_K_vali_seed", seed, ".txt", sep = "")
meth.file.N <- paste("data/meth_N_cpgs_in_KNT_imputed_ClockCpGs_seed", seed, ".txt", sep = "")
cov.file.N  <- paste("data/cov_N_seed", seed, ".txt", sep = "")
meth.file.T <- paste("data/meth_T_cpgs_in_KNT_imputed_ClockCpGs_seed", seed, ".txt", sep = "")
cov.file.T  <- paste("data/cov_T_seed", seed, ".txt", sep = "")

model.residual = FALSE

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
# LOAD METH/COV DATA AND PREDICT
###########################################################################################

##### METH ####
meth.data.K <- read.table(meth.file.K, header = TRUE, sep = '\t')
n.samples.K <- dim(meth.data.K)[[2]] - 1
x <- c(1)
x <- rep(x, n.samples.K)
x <- c("Intercept", x)
df <- data.frame(x)
df <- t(df)
colnames(df) <- colnames(meth.data.K)
meth.data.K <- rbind(meth.data.K, df)
row.names(meth.data.K) <- NULL

meth.data.N <- read.table(meth.file.N, header = TRUE, sep = '\t')
n.samples.N <- dim(meth.data.N)[[2]] - 1
x <- c(1)
x <- rep(x, n.samples.N)
x <- c("Intercept", x)
df <- data.frame(x)
df <- t(df)
colnames(df) <- colnames(meth.data.N)
meth.data.N <- rbind(meth.data.N, df)
row.names(meth.data.N) <- NULL

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

meth.data.N <- meth.data.N[ -c(1) ]
meth.data.T <- meth.data.T[ -c(1) ]
meth <- cbind(meth.data.K, meth.data.N, meth.data.T)
rm(meth.data.K, meth.data.N, meth.data.T); gc()

##### COV #####
cov.K <- read.table(cov.file.K, header = TRUE, row.names = 1, sep = '\t')
cov.N <- read.table(cov.file.N, header = TRUE, row.names = 1, sep = '\t')
cov.T <- read.table(cov.file.T, header = TRUE, row.names = 1, sep = '\t')

K.samples <- colnames(cov.K)
N.samples <- colnames(cov.N)
T.samples <- colnames(cov.T)

covariates <- unique(c(rownames(cov.K), rownames(cov.N), rownames(cov.T)))
for(i in covariates){
  if (!(i %in% row.names(cov.K))) {cov.K[i,] <- NA}
  if (!(i %in% row.names(cov.N))) {cov.N[i,] <- NA}
  if (!(i %in% row.names(cov.T))) {cov.T[i,] <- NA}
}

cov.K <- cov.K[ order(row.names(cov.K)), ]
cov.N <- cov.N[ order(row.names(cov.N)), ]
cov.T <- cov.T[ order(row.names(cov.T)), ]

cov <- cbind(cov.K, cov.N, cov.T)
rm(cov.K, cov.N, cov.T); gc()

##### PREDICT #####
clock.cpg.coef <- read.csv(paste(model.dir, "model_coefficients.csv", sep = ''), stringsAsFactors = FALSE)
clock.cpg.coef <- clock.cpg.coef[ c("model.coefficients.name", "model.coefficients.x") ]
clock.cpg.coef[clock.cpg.coef$model.coefficients.name == "(Intercept)", "model.coefficients.name"] <- "Intercept"

meth <- meth[tolower(order(meth$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth$position == clock.cpg.coef$model.coefficients.name

meth$position <- NULL
X    <- data.matrix(meth)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result <- t(X) %*% beta
result <- sapply(result, anti.trafo)

ages <- as.numeric(as.vector(cov["Age",]))
res  <- result - ages
if(model.residual){
  res <- lm(result ~ ages)$residuals
}

##### STORE RESULTS IN COV #####
cov["DNAm Age",] <- result
cov["DNAm Age Residual",] <- res

##### RESIDUALS #####
residual.K <- as.numeric(as.vector(cov["DNAm Age Residual", K.samples]))
residual.N <- as.numeric(as.vector(cov["DNAm Age Residual", N.samples]))
residual.T <- as.numeric(as.vector(cov["DNAm Age Residual", T.samples]))


##########################################################################
# K PLOTS 
##########################################################################
mean.error.K <- mean(residual.K)
stdev.error.K <- sd(residual.K)
cor.test(sample.ages.K,result.K, method = "pearson")
median(residual.K)
p <- ggplot(data.frame(res = residual.K), aes(x=res)) +
  geom_histogram(binwidth = 5, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 33)
  ) + 
  labs(x = "DNAm Age Acceleration") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Normal Tissue") + 
  annotate("text", x = -16, y = 32, label = paste("Median Residual = ", round(mean.error.K, digits = 1), sep = "")) + 
  annotate("text", x = -16, y = 30.5, label = paste("St. Dev Residual = ", round(stdev.error.K, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "TissueStudies/residual_hist_K.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "TissueStudies/MethAgevsSampleAge_K.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.K, 
     result.K, 
     main="Methlyation Age vs Sample Age in Normal Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.K~sample.ages.K), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.K~sample.ages.K))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# N PLOTS 
##########################################################################
mean.error.N <- mean(residual.N)
stdev.error.N <- sd(residual.N)

p <- ggplot(data.frame(res = residual.N), aes(x=res)) +
  geom_histogram(binwidth = 5, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 15)
  ) + 
  labs(x = "DNAm Age Acceleration") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Adjacent Normal Tissue") + 
  annotate("text", x = -21, y = 14, label = paste("Median Residual = ", round(mean.error.N, digits = 1), sep = "")) + 
  annotate("text", x = -21, y = 13, label = paste("St. Dev Residual = ", round(stdev.error.N, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "TissueStudies/residual_hist_N.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "TissueStudies/MethAgevsSampleAge_N.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.N, 
     result.N, 
     main="Methlyation Age vs Sample Age in Adjacent Normal Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.N~sample.ages.N), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.N~sample.ages.N))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()


##########################################################################
# T PLOTS 
##########################################################################
mean.error.T <- mean(residual.T)
stdev.error.T <- sd(residual.T)

p <- ggplot(data.frame(res = residual.T), aes(x=res)) +
  geom_histogram(binwidth = 10, color="black", fill="white") +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 33)
  ) + 
  labs(x = "DNAm Age Acceleration") + 
  labs(y = "Frequency") + 
  labs(title = "Prediction Residuals in Tumor Tissue") + 
  annotate("text", x = 100, y = 32, label = paste("Median Residual = ", round(mean.error.T, digits = 1), sep = "")) + 
  annotate("text", x = 100, y = 30.5, label = paste("St. Dev Residual = ", round(stdev.error.T, digits = 1), sep = "")) +
  geom_vline(xintercept = 0, color = "red", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


png( paste(model.dir, "TissueStudies/residual_hist_T.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "TissueStudies/MethAgevsSampleAge_T.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.T, 
     result.T, 
     main="Methlyation Age vs Sample Age in Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,87), 
     xaxs="i",
     ylim=c(15,87), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(result.T~sample.ages.T), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
rsq <- summary(lm(result.T~sample.ages.T))$r.squared
text(24,80, paste("r^2 = ", round(rsq, digits = 4), sep = ""))
dev.off()



##########################################################################
# MERGED RESIDUAL PLOT 
##########################################################################
median.error.K <- median(residual.K)
median.error.N <- median(residual.N)
median.error.T <- median(residual.T)

tmp.K <- data.frame(residual.K)
tmp.K$type <- "K"
colnames(tmp.K) <- c("res", "ttype")

tmp.N <- data.frame(residual.N)
tmp.N$type <- "N"
colnames(tmp.N) <- c("res", "ttype")

tmp.T <- data.frame(residual.T)
tmp.T$type <- "T"
colnames(tmp.T) <- c("res", "ttype")

df.KT <- rbind(tmp.K, tmp.N, tmp.T)
rm(tmp.K)
rm(tmp.N)
rm(tmp.T)
gc()

p <- ggplot(df.KT, aes(x = res, stat(density), color = ttype, linetype = ttype)) +
  geom_freqpoly(size = 1.2, binwidth = 20) +
  scale_linetype_manual(
    name = "Tissue Type", 
    values = c("solid", "dashed", "twodash"),
    labels = c("K", "N", "T")
  ) + 
  scale_color_manual(
    name = "Tissue Type", 
    values = c("black","darkorange3","deepskyblue4"), 
    labels = c("K", "N", "T")
  ) +
  #scale_x_continuous(
  #  expand=c(0, 0),
  #  limits = c(-60, 200)
  #) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 0.05)
  ) + 
  labs(x = "DNAm Age Acceleration [Years]") + 
  labs(y = "Frequency [Arbitrary Units]") + 
  labs(title = "DNAm Age Acceleration for Tissue Types") + 
  annotate(
    "text", x = 130, y = 0.037, 
    label = paste("Median Residual (K) = ", round(median.error.K, digits = 1), sep = ""), 
    color = "black"
  ) + 
  #annotate(
  #  "text", x = 100, y = 24, 
  #  label = paste("St. Dev Residual (K) = ", round(stdev.error.K, digits = 1), sep = ""), 
  #  color = "darkorange3"
  #) +
  annotate(
    "text", x = 130, y = 0.034, 
    label = paste("Median Residual (N) = ", round(median.error.N, digits = 1), sep = ""), 
    color = "darkorange3"
  ) + 
  #annotate(
  #  "text", x = 100, y = 21, 
  #  label = paste("St. Dev Residual (N) = ", round(stdev.error.N, digits = 1), sep = ""), 
  #  color = "black"
  #) +
  annotate(
    "text", x = 130, y = 0.031, 
    label = paste("Median Residual (T) = ", round(median.error.T, digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) + 
  #annotate(
  #  "text", x = 100, y = 18, 
  #  label = paste("St. Dev Residual (T) = ", round(stdev.error.T, digits = 1), sep = ""), 
  #  color = "deepskyblue4"
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


png( paste(model.dir, "TissueStudies/residual_hist_KNT.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

png( paste(model.dir, "TissueStudies/MethAgevsSampleAge_KNT.png", sep = ''), width = 500, height = 500, units = "px" )
plot(sample.ages.K, 
     result.K, 
     main="Methlyation Age vs Sample Age in Normal and Tumor Tissue", 
     xlab="Sample Age ", 
     ylab="Methylation Age ", 
     pch=19,
     xlim=c(15,100), 
     xaxs="i",
     ylim=c(15,200), 
     yaxs="i",
     tck = 0.02,
     col = "darkorange3"
) 
points(sample.ages.T, result.T, col = "deepskyblue4")
abline(lm(result.K~sample.ages.K), col="darkorange3") # regression line (y~x) 
points(sample.ages.N, result.N, col = "black", pch = 19)
legend(
  20, 190, 
  title = "Tissue Type",
  legend = c("K", "N", "T"), 
  col = c("darkorange3", "black", "deepskyblue4"), 
  lty = c(0,0), 
  pch = c(19,19, 1), 
  bg='white', 
  bty = "n"
)
dev.off()

wilcox.test(residual.K, mu=0, conf.int = TRUE)
wilcox.test(residual.N, mu=0, conf.int = TRUE)
wilcox.test(residual.T, mu=0, conf.int = TRUE)

