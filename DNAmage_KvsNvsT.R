rm(list=ls()); gc();
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(RColorBrewer)
source("plot_functions.R")

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
result.K <- as.numeric(as.vector(cov["DNAm Age", K.samples]))
result.N <- as.numeric(as.vector(cov["DNAm Age", N.samples]))
result.T <- as.numeric(as.vector(cov["DNAm Age", T.samples]))

residual.K <- as.numeric(as.vector(cov["DNAm Age Residual", K.samples]))
residual.N <- as.numeric(as.vector(cov["DNAm Age Residual", N.samples]))
residual.T <- as.numeric(as.vector(cov["DNAm Age Residual", T.samples]))

age.K <- as.numeric(as.vector(cov["Age", K.samples]))
age.N <- as.numeric(as.vector(cov["Age", N.samples]))
age.T <- as.numeric(as.vector(cov["Age", T.samples]))

##########################################################################
# MERGED RESIDUAL PLOT 
##########################################################################
tmp.K <- data.frame(age.K, result.K, residual.K)
tmp.K$type <- "K"
colnames(tmp.K) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.N <- data.frame(age.N, result.N, residual.N)
tmp.N$type <- "N"
colnames(tmp.N) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.T <- data.frame(age.T, result.T, residual.T)
tmp.T$type <- "T"
colnames(tmp.T) <- c("Chrono.age", "DNAm.age", "res", "ttype")

df.KNT <- rbind(tmp.K, tmp.N, tmp.T)
df.KNT$ttype <- factor(df.KNT$ttype, levels = c("K","N","T"))
rm(tmp.K); rm(tmp.N); rm(tmp.T); gc()

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.KNT, 
  bw = 20, 
  legname = "Tissue Type", 
  linetypes =  c("solid", "dashed", "twodash"),
  colors = c("black","darkorange3","deepskyblue4"),
  labels = c("K", "N", "T"),
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Tissue Types",
  annot.x = 110, 
  annot.y = 0.037,
  annot.sep = 0.003, 
  x.min = -50, 
  x.max = 150,
  y.min = 0, 
  y.max = 0.05, 
  leg.x = 0.8, 
  leg.y = 0.87
); p

png( paste(model.dir, "TissueStudies/residual_hist_KNT.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.K
r.list[[2]] <- residual.N
r.list[[3]] <- residual.T
p <- accel.box.plot(
  df.KNT, 
  residuals = r.list, 
  width = 0.75,
  x.label = "Tissue Type", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Accelerati for Tissue Types"
); p

png( paste(model.dir, "TissueStudies/boxplot_KNT.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.KNT, 
  legname = "Tissue Type", 
  colors = c("black","darkorange3","deepskyblue4"), 
  labels = c("K", "N", "T"), 
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.15, 
  leg.y = 0.87
); p

png( paste(model.dir, "TissueStudies/MethAgevsSampleAge_KNT.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()



##########################################################################
# STATISTICAL TESTS
##########################################################################

##### ALLOW OUTLIERS #####
wilcox.test(residual.K, mu=0, conf.int = TRUE)
wilcox.test(residual.N, mu=0, conf.int = TRUE)
wilcox.test(residual.T, mu=0, conf.int = TRUE)


##### REMOVE OUTLIERS #####
wilcox.test(residual.K[ !(residual.K %in% boxplot.stats(residual.K)$out) ], mu=0)
wilcox.test(residual.N[ !(residual.N %in% boxplot.stats(residual.N)$out) ], mu=0)
wilcox.test(residual.T[ !(residual.T %in% boxplot.stats(residual.T)$out) ], mu=0)