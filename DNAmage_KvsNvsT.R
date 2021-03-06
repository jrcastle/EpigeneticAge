rm(list=ls()); gc();
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(RColorBrewer)
library(car)
source("plot_functions.R")

MODEL          <- 4 # 1 = Old Model, 2 = New Model, 3 = Horvath Model
NORMALIZE      <- FALSE
OMIT.SSIMP     <- FALSE
RM.NA          <- TRUE
model.residual <- FALSE

if(MODEL == 1){
  seed        <- "123"
  model.dir   <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')
  meth.file.K <- paste("data/meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs_seed", seed, ".txt", sep = "")
  cov.file.K  <- paste("data/cov_K_vali_seed", seed, ".txt", sep = "")
  meth.file.N <- paste("data/meth_N_cpgs_in_KNT_imputed_ClockCpGs_seed", seed, ".txt", sep = "")
  cov.file.N  <- paste("data/cov_N_seed", seed, ".txt", sep = "")
  meth.file.T <- paste("data/meth_T_cpgs_in_KNT_imputed_ClockCpGs_seed", seed, ".txt", sep = "")
  cov.file.T  <- paste("data/cov_T_seed", seed, ".txt", sep = "")
}
if(MODEL == 2){
  seed        <- "123"
  model.dir   <- paste("gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_seed", seed, "/", sep = '')
  meth.file.K <- paste("data/meth_K_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_vali_ClockCpGs_seed", seed, ".txt", sep = "")
  cov.file.K  <- paste("data/cov_K_vali_seed", seed, ".txt", sep = "")
  meth.file.N <- paste("data/meth_N_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_ClockCpGs_seed", seed, ".txt", sep = "")
  cov.file.N  <- paste("data/cov_N_seed", seed, ".txt", sep = "")
  meth.file.T <- paste("data/meth_T_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_ClockCpGs_seed", seed, ".txt", sep = "")
  cov.file.T  <- paste("data/cov_T_seed", seed, ".txt", sep = "")
}
if(MODEL == 3){
  model.dir <- "HorvathClock/"
  meth.file.K <- "data/meth_K_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_CGHorvathClock.txt"
  cov.file.K  <- "data/cov_K.txt"
  meth.file.N <- "data/meth_N_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_CGHorvathClock.txt"
  cov.file.N  <- "data/cov_N_seed123.txt"
  meth.file.T <- "data/meth_T_gt10R_AddMissHorvCpGs_KNT_KnnImp_SSImpWgtd_FINAL_CGHorvathClock.txt"
  cov.file.T  <- "data/cov_T_seed123.txt"
}
if(MODEL == 4){
  model.dir <- "HorvathClock/"
  meth.file.K <- "data/meth_K_gt10R_AddMissHorvCpGs_KNT_CGHorvClock.txt"
  cov.file.K  <- "data/cov_K.txt"
  meth.file.N <- "data/meth_N_gt10R_AddMissHorvCpGs_KNT_CGHorvClock.txt"
  cov.file.N  <- "data/cov_N_seed123.txt"
  meth.file.T <- "data/meth_T_gt10R_AddMissHorvCpGs_KNT_CGHorvClock.txt"
  cov.file.T  <- "data/cov_T_seed123.txt"
}

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

meth.data.K <- read.table(meth.file.K, header = TRUE, sep = '\t')
meth.data.N <- read.table(meth.file.N, header = TRUE, sep = '\t')
meth.data.T <- read.table(meth.file.T, header = TRUE, sep = '\t')

if((MODEL == 3 | MODEL == 4) & NORMALIZE){
  source("~/Documents/EpigeneticAge/AF24_RCodeForNormalizingDNAMethylationData.txt")
  library(dynamicTreeCut)
  prbAnn21kdatMeth <- read.csv("~/Documents/EpigeneticAge/AF22_AdditionalProbeAnnotationFileForRTutorial.csv")
  prbAnn21kdatMeth <- prbAnn21kdatMeth[which(prbAnn21kdatMeth$Name %in% meth.data.K$position),]
  prbAnn21kdatMeth <- prbAnn21kdatMeth[order(prbAnn21kdatMeth$Name),] 
  
  meth.data.K      <- meth.data.K[order(meth.data.K$position),]
  meth.data.N      <- meth.data.N[order(meth.data.N$position),]
  meth.data.T      <- meth.data.T[order(meth.data.T$position),]
  cpgs             <- as.character(meth.data.K$position)
  
  # K
  df.tmp                <- t(meth.data.K[,-1])
  colnames(df.tmp)      <- cpgs
  meth.data.K           <- data.frame(t(BMIQcalibration(datM = df.tmp, 
                                                        goldstandard.beta = prbAnn21kdatMeth$goldstandard2,
                                                        plots = FALSE
                                                        )))
  meth.data.K$position  <- rownames(meth.data.K)
  rownames(meth.data.K) <- NULL
  meth.data.K           <- meth.data.K[,c(ncol(meth.data.K), 1:ncol(meth.data.K)-1)]
  
  # N
  df.tmp                <- t(meth.data.N[,-1])
  colnames(df.tmp)      <- cpgs
  meth.data.N           <- data.frame(t(BMIQcalibration(datM = df.tmp, 
                                                        goldstandard.beta = prbAnn21kdatMeth$goldstandard2,
                                                        plots = FALSE
                                                        )))
  meth.data.N$position  <- rownames(meth.data.N)
  rownames(meth.data.N) <- NULL
  meth.data.N           <- meth.data.N[,c(ncol(meth.data.N), 1:ncol(meth.data.N)-1)]
  
  # T 
  df.tmp                <- t(meth.data.T[,-1])
  colnames(df.tmp)      <- cpgs
  meth.data.T           <- data.frame(t(BMIQcalibration(datM = df.tmp, 
                                                        goldstandard.beta = prbAnn21kdatMeth$goldstandard2,
                                                        plots = FALSE
                                                        )))
  meth.data.T$position  <- rownames(meth.data.T)
  rownames(meth.data.T) <- NULL
  meth.data.T           <- meth.data.T[,c(ncol(meth.data.T), 1:ncol(meth.data.T)-1)]
  
}

if(MODEL == 3 & OMIT.SSIMP){
  bad.cpgs <- read.table("data/MissingHorvathClockCpGs_CGNumber.txt")
  bad.cpgs <- as.character(bad.cpgs$V1)
  meth.data.K[which(meth.data.K$position %in% bad.cpgs), c(2:ncol(meth.data.K))] <- 0
  meth.data.N[which(meth.data.N$position %in% bad.cpgs), c(2:ncol(meth.data.N))] <- 0
  meth.data.T[which(meth.data.T$position %in% bad.cpgs), c(2:ncol(meth.data.T))] <- 0
}

if(MODEL == 4 & RM.NA){
  meth.data.K[is.na(meth.data.K)] <- 0
  meth.data.N[is.na(meth.data.N)] <- 0
  meth.data.T[is.na(meth.data.T)] <- 0
}

# K
n.samples.K <- dim(meth.data.K)[[2]] - 1
x <- c(1)
x <- rep(x, n.samples.K)
x <- c("Intercept", x)
df <- data.frame(x)
df <- t(df)
colnames(df) <- colnames(meth.data.K)
meth.data.K <- rbind(meth.data.K, df)
row.names(meth.data.K) <- NULL

# N
n.samples.N <- dim(meth.data.N)[[2]] - 1
x <- c(1)
x <- rep(x, n.samples.N)
x <- c("Intercept", x)
df <- data.frame(x)
df <- t(df)
colnames(df) <- colnames(meth.data.N)
meth.data.N <- rbind(meth.data.N, df)
row.names(meth.data.N) <- NULL

# T 
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

cov.K['Tissue.Type',] = "K"
cov.N['Tissue.Type',] = "N"
cov.T['Tissue.Type',] = "T"

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
  res2 <- lm(result ~ ages)$residuals
  cov["DNAm Age Model Residual",] <- res2
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
  x.label = "Epigenetic Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "Tissue-specific Epigenetic Age Acceleration",
  annot.x = 103, 
  annot.y = 0.036,
  annot.sep = 0.003,
  annot.size = 6,
  x.min = -50, 
  x.max = 149,
  y.min = 0, 
  y.max = 0.05, 
  leg.x = 0.8, 
  leg.y = 0.87
); p

tiff( paste(model.dir, "TissueStudies/residual_hist_KNT.tiff", sep = ''),  width = 2100, height = 2100, units = "px",res = 300 )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.K
r.list[[2]] <- residual.N
r.list[[3]] <- residual.T
q <- accel.box.plot(
  df.KNT, 
  residuals = r.list, 
  width = 0.75,
  x.label = "Tissue Type", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "Epigenetic Age Acceleration for Tissue Types",
  leg.x = 0.1, 
  leg.y = 0.95
); q

tiff( paste(model.dir, "TissueStudies/boxplot_KNT.tiff", sep = ''),  width = 2100, height = 2100, units = "px",res = 300 )
q
dev.off()

##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.KNT, 
  legname = "Tissue Type", 
  colors = c("black","darkorange3","deepskyblue4"), 
  symbol.shapes = c(16, 1, 15),
  labels = c("K", "N", "T"), 
  x.label = "Chronological Age [Years]", 
  y.label = "Epigenetic Age [Years]", 
  title = "Epigenetic Age vs Chronological Age",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.15, 
  leg.y = 0.87
); p

if( model.residual ){
  p <- p + 
    geom_abline(
      slope = lm(result ~ ages)$coefficients[[2]], 
      intercept = lm(result ~ ages)$coefficients[[1]],
      linetype = "dotted"
    )
}

tiff( paste(model.dir, "TissueStudies/MethAgevsSampleAge_KNT.tiff", sep = ''),  width = 2100, height = 2100, units = "px",res = 300 )
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

##### LEVINE'S TEST #####
res.K <- residual.K[ !(residual.K %in% boxplot.stats(residual.K)$out) ]
res.N <- residual.N[ !(residual.N %in% boxplot.stats(residual.N)$out) ]
res.T <- residual.T[ !(residual.T %in% boxplot.stats(residual.T)$out) ]

# K-N
y <- c(res.K, res.N)
grp <- as.factor( c(rep("K", length(res.K)), rep("N", length(res.N))) )
leveneTest(y, grp)

# N-T
y <- c(res.N, res.T)
grp <- as.factor( c(rep("N", length(res.N)), rep("T", length(res.T))) )
leveneTest(y, grp)

# K-T
y <- c(res.K, res.T)
grp <- as.factor( c(rep("K", length(res.K)), rep("T", length(res.T))) )
leveneTest(y, grp)

##### ONE-WAY ANOVA #####
ag.K <- age.K[ !(residual.K %in% boxplot.stats(residual.K)$out) ]
ag.N <- age.N[ !(residual.N %in% boxplot.stats(residual.N)$out) ]
ag.T <- age.T[ !(residual.T %in% boxplot.stats(residual.T)$out) ]

tmp.K <- data.frame(res.K, ag.K)
tmp.K$ttype <- "K"
colnames(tmp.K) <- c("res", "age", "ttype")

tmp.N <- data.frame(res.N, ag.N)
tmp.N$ttype <- "N"
colnames(tmp.N) <- c("res", "age", "ttype")

tmp.T <- data.frame(res.T, ag.T)
tmp.T$ttype <- "T"
colnames(tmp.T) <- c("res", "age", "ttype")

df.merge <- rbind(tmp.K, tmp.N, tmp.T)
df.merge$age <- as.factor(df.merge$age)
df.merge$ttype <- as.factor(df.merge$ttype)
rm(tmp.K, tmp.N, tmp.T); gc()

fit <- aov(res ~ age + ttype, data = df.merge)
TukeyHSD(fit, which = 'ttype')


##### TTest #####
if(model.residual){

  res.K <- as.numeric(as.vector(cov["DNAm Age Model Residual", which(cov["Tissue.Type",] == "K")]))
  res.N <- as.numeric(as.vector(cov["DNAm Age Model Residual", which(cov["Tissue.Type",] == "N")]))
  res.T <- as.numeric(as.vector(cov["DNAm Age Model Residual", which(cov["Tissue.Type",] == "T")]))

  res.K <- res.K[ !(res.K %in% boxplot.stats(res.K)$out) ]
  res.N <- res.N[ !(res.N %in% boxplot.stats(res.N)$out) ]
  res.T <- res.T[ !(res.T %in% boxplot.stats(res.T)$out) ]

  t.test(res.N, res.T, var.equal = FALSE)

}

