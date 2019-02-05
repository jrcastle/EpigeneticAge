rm(list=ls()); gc();
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(RColorBrewer)
source("plot_functions.R")

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
# LOAD METH/COV DATA AND PREDICT
###########################################################################################

##### METH ####
meth <- read.table(meth.file.T, header = TRUE, sep = '\t')
n.samples.T <- dim(meth)[[2]] - 1
x <- c(1)
x <- rep(x, n.samples.T)
x <- c("Intercept", x)
df <- data.frame(x)
df <- t(df)
colnames(df) <- colnames(meth)
meth <- rbind(meth, df)
row.names(meth) <- NULL

##### COV #####
cov <- read.table(cov.file.T, header = TRUE, row.names = 1, sep = '\t')

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
#res <- lm(result ~ ages)$residuals
res <- result - ages

##### STORE RESULTS IN COV #####
cov["DNAm Age",] <- result
cov["DNAm Age Residual",] <- res


###########################################################################################
# SPLIT SAMPLES
###########################################################################################
tmp <- data.frame(t(cov))

##### SAMPLES #####
Her2p.White.samples      <- rownames(tmp[ which(tmp$Her2pos == 1 & tmp$RaceWhite == 1 ), ])
ERPRpHer2n.White.samples <- rownames(tmp[ which( (tmp$ERpos == 1 | tmp$PRpos == 1) & tmp$Her2pos == 0 & tmp$RaceWhite == 1), ])
ERPRnHer2n.White.samples <- rownames(tmp[ which(tmp$ERpos == 0 & tmp$PRpos == 0 & tmp$Her2pos == 0 & tmp$RaceWhite == 1), ])

Her2p.Black.samples      <- rownames(tmp[ which(tmp$Her2pos == 1 & tmp$RaceWhite == 0 ), ])
ERPRpHer2n.Black.samples <- rownames(tmp[ which( (tmp$ERpos == 1 | tmp$PRpos == 1) & tmp$Her2pos == 0 & tmp$RaceWhite == 0), ])
ERPRnHer2n.Black.samples <- rownames(tmp[ which(tmp$ERpos == 0 & tmp$PRpos == 0 & tmp$Her2pos == 0 & tmp$RaceWhite == 0), ])

##### AGES #####
age.White.H2p     <- as.numeric(as.vector(cov["Age", Her2p.White.samples]))
age.White.EPpHn   <- as.numeric(as.vector(cov["Age", ERPRpHer2n.White.samples]))
age.White.EPnHn   <- as.numeric(as.vector(cov["Age", ERPRnHer2n.White.samples]))

age.Black.H2p     <- as.numeric(as.vector(cov["Age", Her2p.Black.samples]))
age.Black.EPpHn   <- as.numeric(as.vector(cov["Age", ERPRpHer2n.Black.samples]))
age.Black.EPnHn   <- as.numeric(as.vector(cov["Age", ERPRnHer2n.Black.samples]))

##### RESIDUALS #####
residual.White.H2p     <- as.numeric(as.vector(cov["DNAm Age Residual", Her2p.White.samples]))
residual.White.EPpHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRpHer2n.White.samples]))
residual.White.EPnHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnHer2n.White.samples]))

residual.Black.H2p     <- as.numeric(as.vector(cov["DNAm Age Residual", Her2p.Black.samples]))
residual.Black.EPpHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRpHer2n.Black.samples]))
residual.Black.EPnHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnHer2n.Black.samples]))


###########################################################################################
# REMOVE OUTLIERS
###########################################################################################

##### OUTLIER POSITIONS #####
outlier.White.H2p   <- match(c(boxplot.stats(residual.White.H2p)$out), residual.White.H2p)
outlier.White.EPpHn <- match(c(boxplot.stats(residual.White.EPpHn)$out), residual.White.EPpHn)
outlier.White.EPnHn <- match(c(boxplot.stats(residual.White.EPnHn)$out), residual.White.EPnHn)

outlier.Black.H2p   <- match(c(boxplot.stats(residual.Black.H2p)$out), residual.Black.H2p)
outlier.Black.EPpHn <- match(c(boxplot.stats(residual.Black.EPpHn)$out), residual.Black.EPpHn)
outlier.Black.EPnHn <- match(c(boxplot.stats(residual.Black.EPnHn)$out), residual.Black.EPnHn)

##### AGE #####
if( length(outlier.White.H2p)   > 0) { age.White.H2p   <- age.White.H2p[-outlier.White.H2p] }
if( length(outlier.White.EPpHn) > 0) { age.White.EPpHn <- age.White.EPpHn[-outlier.White.EPpHn] }
if( length(outlier.White.EPnHn) > 0) { age.White.EPnHn <- age.White.EPnHn[-outlier.White.EPnHn] }

if( length(outlier.Black.H2p)   > 0) { age.Black.H2p   <- age.Black.H2p[-outlier.Black.H2p] }
if( length(outlier.Black.EPpHn) > 0) { age.Black.EPpHn <- age.Black.EPpHn[-outlier.Black.EPpHn] }
if( length(outlier.Black.EPnHn) > 0) { age.Black.EPnHn <- age.Black.EPnHn[-outlier.Black.EPnHn] }


##### RESIDUALS #####
if( length(outlier.White.H2p)   > 0) { residual.White.H2p     <- residual.White.H2p[-outlier.White.H2p] }
if( length(outlier.White.EPpHn) > 0) { residual.White.EPpHn   <- residual.White.EPpHn[-outlier.White.EPpHn] }
if( length(outlier.White.EPnHn) > 0) { residual.White.EPnHn   <- residual.White.EPnHn[-outlier.White.EPnHn] }

if( length(outlier.Black.H2p)   > 0) { residual.Black.H2p     <- residual.Black.H2p[-outlier.Black.H2p] }
if( length(outlier.Black.EPpHn) > 0) { residual.Black.EPpHn   <- residual.Black.EPpHn[-outlier.Black.EPpHn] }
if( length(outlier.Black.EPnHn) > 0) { residual.Black.EPnHn   <- residual.Black.EPnHn[-outlier.Black.EPnHn] }


###########################################################################################
# LINEAR MODELS
###########################################################################################

# Her2p
Her2p.age <- c(age.White.H2p, age.Black.H2p)
Her2p.res <- c(residual.White.H2p, residual.Black.H2p)
Her2p.rce <- c( rep(1, length(age.White.H2p)), rep(0, length(age.Black.H2p)) )
Her2p.rce <- as.factor(Her2p.rce)
Her2p.lm  <- lm(Her2p.res ~ Her2p.age + Her2p.rce)

# EPpHn
EPpHn.age <- c(age.White.EPpHn, age.Black.EPpHn)
EPpHn.res <- c(residual.White.EPpHn, residual.Black.EPpHn)
EPpHn.rce <- c( rep(1, length(age.White.EPpHn)), rep(0, length(age.Black.EPpHn)) )
EPpHn.rce <- as.factor(EPpHn.rce)
EPpHn.lm  <- lm(EPpHn.res ~ EPpHn.age + EPpHn.rce)

# EPnHn
EPnHn.age <- c(age.White.EPnHn, age.Black.EPnHn)
EPnHn.res <- c(residual.White.EPnHn, residual.Black.EPnHn)
EPnHn.rce <- c( rep(1, length(age.White.EPnHn)), rep(0, length(age.Black.EPnHn)) )
EPnHn.rce <- as.factor(EPnHn.rce)
EPnHn.lm  <- lm(EPnHn.res ~ EPnHn.age + EPnHn.rce)

# Data frames
df.Her2p.lm <- data.frame(summary(Her2p.lm)$coefficients)
df.EPpHn.lm <- data.frame(summary(EPpHn.lm)$coefficients)
df.EPnHn.lm <- data.frame(summary(EPnHn.lm)$coefficients)


###########################################################################################
# Her2+ Boxplot
###########################################################################################
df1           <- data.frame(residual.White.H2p)
df1$ttype     <- "Caucasian"
colnames(df1) <- c("res", "ttype")
df2           <- data.frame(residual.Black.H2p)
df2$ttype     <- "African American"
colnames(df2) <- c("res", "ttype")
df.tmp.Age    <- rbind(df1, df2)

r.list <- list()
r.list[[1]] <- residual.Black.H2p
r.list[[2]] <- residual.White.H2p
p <- accel.box.plot(
  df.tmp.Age, 
  residuals = r.list, 
  width = 0.6,
  x.label = "Race", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "Her2+",
  leg.x = 0.12,
  leg.y = 0.93
) + 
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(-25, 50)
  ) +
  annotate(
    "text", x = 1, y = -22, size = 6,
    label = paste("n = ", length(residual.Black.H2p), sep = '')
  ) + 
  annotate(
    "text", x = 2, y = -22, size = 6,
    label = paste("n = ", length(residual.White.H2p), sep = '')
  ) + 
  annotate(
    "text", x = 0.5, y = 48, 
    label = "a)", 
    size = 6,
    fontface = 2
  ) 
#p <- p +
#  annotate(
#    "text", x = 0.75, y = 18, 
#    label = "Her2+ Samples"
#  ) + 
#  annotate(
#    "text", x = 0.75, y = 15, 
#    label = paste("beta = ", round(df.Her2p.lm[3,1],  digits = 3), sep = "") 
#  ) + 
#  annotate(
#    "text", x = 0.75, y = 12, 
#    label = paste("p = ", round(df.Her2p.lm[3,4],  digits = 3), sep = "")
#  ); p

tiff( paste(model.dir, "CancerStudies/WhiteBlack_Her2p.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()
  

###########################################################################################
#  EPpHn Boxplot
###########################################################################################
df1           <- data.frame(residual.White.EPpHn)
df1$ttype     <- "Caucasian"
colnames(df1) <- c("res", "ttype")
df2           <- data.frame(residual.Black.EPpHn)
df2$ttype     <- "African American"
colnames(df2) <- c("res", "ttype")
df.tmp.Age    <- rbind(df1, df2)

r.list <- list()
r.list[[1]] <- residual.Black.EPpHn
r.list[[2]] <- residual.White.EPpHn
p <- accel.box.plot(
  df.tmp.Age, 
  residuals = r.list, 
  width = 0.6,
  x.label = "Race", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "HR+, Her2-",
  leg.x = 0.24,
  leg.y = 0.93,
  show.leg = FALSE
) + 
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(-30, 90)
  ) +
  annotate(
    "text", x = 1, y = -25, size = 6,
    label = paste("n = ", length(residual.Black.EPpHn), sep = '')
  ) + 
  annotate(
    "text", x = 2, y = -25, size = 6,
    label = paste("n = ", length(residual.White.EPpHn), sep = '')
  ) +
  annotate(
    "text", x = 0.5, y = 87, 
    label = "b)", 
    size = 6,
    fontface = 2
  ) 
#p <- p +
#  annotate(
#    "text", x = 2.25, y = 70, 
#    label = "ER/PR+ Her2- Samples"
#  ) + 
#  annotate(
#    "text", x = 2.25, y = 66, 
#    label = paste("beta = ", round(df.EPpHn.lm[3,1],  digits = 3), sep = "") 
#  ) + 
#  annotate(
#    "text", x = 2.25, y = 62, 
#    label = paste("p = ", round(df.EPpHn.lm[3,4],  digits = 3), sep = "")
#  ); p

tiff( paste(model.dir, "CancerStudies/WhiteBlack_EPpHn.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()


###########################################################################################
#  EPnHn Boxplot
###########################################################################################
df1           <- data.frame(residual.White.EPnHn)
df1$ttype     <- "Caucasian"
colnames(df1) <- c("res", "ttype")
df2           <- data.frame(residual.Black.EPnHn)
df2$ttype     <- "African American"
colnames(df2) <- c("res", "ttype")
df.tmp.Age    <- rbind(df1, df2)

r.list <- list()
r.list[[1]] <- residual.Black.EPnHn
r.list[[2]] <- residual.White.EPnHn
p <- accel.box.plot(
  df.tmp.Age, 
  residuals = r.list, 
  width = 0.6,
  x.label = "Race", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "Triple Negative",
  leg.x = 0.24,
  leg.y = 0.93,
  show.leg = FALSE
) + 
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(-15, 45)
  ) +
  annotate(
    "text", x = 1, y = -12, size = 6,
    label = paste("n = ", length(residual.Black.EPnHn), sep = '')
  ) + 
  annotate(
    "text", x = 2, y = -12, size = 6,
    label = paste("n = ", length(residual.White.EPnHn), sep = '')
  ) +
  annotate(
    "text", x = 0.5, y = 43, 
    label = "c)", 
    size = 6,
    fontface = 2
  ) 

#p <- p + 
#  annotate(
#    "text", x = 2.25, y = 70, 
#    label = "ER- PR- Her2- Samples"
#  ) + 
#  annotate(
#    "text", x = 2.25, y = 66, 
#    label = paste("beta = ", round(df.EPnHn.lm[3,1],  digits = 3), sep = "") 
#  ) + 
#  annotate(
#    "text", x = 2.25, y = 62, 
#    label = paste("p = ", round(df.EPnHn.lm[3,4],  digits = 3), sep = "")
#  ); p

tiff( paste(model.dir, "CancerStudies/WhiteBlack_EPnHn.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()


###########################################################################################
# MODEL SUMMARIES
###########################################################################################
summary(Her2p.lm)
summary(EPpHn.lm)
summary(EPnHn.lm)
