rm(list=ls()); gc();
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(RColorBrewer)
source("plot_functions.R")

model.dir     <- "HorvathClock/"
result.file.T <- paste(model.dir, "Output_T.csv", sep = '')
cov.file.T    <- "data/cov_T_seed123.txt"
 

###########################################################################################
# LOAD RESULTS/COV DATA 
###########################################################################################

##### METH T ####
result <- read.csv(result.file.T, header = TRUE)

##### COV T #####
cov <- read.table(cov.file.T, header = TRUE, row.names = 1, sep = '\t')
DNAmAge <- as.numeric(as.vector(result$DNAmAge))
DNAmAgeResidual <- as.numeric(as.vector(cov["Age",]))
DNAmAgeResidual <- DNAmAge - DNAmAgeResidual
cov["DNAm Age",] <- DNAmAge
cov["DNAm Age Residual",] <- DNAmAgeResidual


###########################################################################################
# SPLIT BY CANCER SUBTYPE
###########################################################################################
tmp <- data.frame(t(cov))

##### DIVIDE SAMPLES BY TYPE (V1) #####
LuminalA.samples       <-rownames(tmp[ which(tmp$Cancer.subtypeLuminal.A == 1), ])
LuminalB.samples       <-rownames(tmp[ which(tmp$Cancer.subtypeLuminal.B == 1), ])
Triplenegative.samples <-rownames(tmp[ which(tmp$Cancer.subtypeTriplenegative == 1), ])

##### DIVIDE SAMPLES BY TYPE (V2) #####
ERPRposHer2pos.samples <- rownames(tmp[ which(tmp$ERPRposHer2pos == 1), ])
ERPRposHer2neg.samples <- rownames(tmp[ which(tmp$ERPRposHer2neg == 1), ])
ERPRnegHer2pos.samples <- rownames(tmp[ which(tmp$ERPRnegHer2pos == 1), ])
ERPRnegHer2neg.samples <- rownames(tmp[ which(tmp$ERPRnegHer2neg == 1), ])

##### DIVIDE SAMPLES BY TYPE (V3) #####
Her2p.samples      <- rownames(tmp[ which(tmp$Her2pos == 1), ])
ERPRpHer2n.samples <- rownames(tmp[ which( (tmp$ERpos == 1 | tmp$PRpos == 1) & tmp$Her2pos == 0), ])
ERPRnHer2n.samples <- rownames(tmp[ which(tmp$ERpos == 0 & tmp$PRpos == 0 & tmp$Her2pos == 0), ])

Her2p.White.samples      <- rownames(tmp[ which(tmp$Her2pos == 1 & tmp$RaceWhite == 1 ), ])
ERPRpHer2n.White.samples <- rownames(tmp[ which( (tmp$ERpos == 1 | tmp$PRpos == 1) & tmp$Her2pos == 0 & tmp$RaceWhite == 1), ])
ERPRnHer2n.White.samples <- rownames(tmp[ which(tmp$ERpos == 0 & tmp$PRpos == 0 & tmp$Her2pos == 0 & tmp$RaceWhite == 1), ])

Her2p.Black.samples      <- rownames(tmp[ which(tmp$Her2pos == 1 & tmp$RaceWhite == 0 ), ])
ERPRpHer2n.Black.samples <- rownames(tmp[ which( (tmp$ERpos == 1 | tmp$PRpos == 1) & tmp$Her2pos == 0 & tmp$RaceWhite == 0), ])
ERPRnHer2n.Black.samples <- rownames(tmp[ which(tmp$ERpos == 0 & tmp$PRpos == 0 & tmp$Her2pos == 0 & tmp$RaceWhite == 0), ])


##### DIVIDE SAMPLES BY TUMOR GRADE #####
grade3.5.samples <- rownames(tmp[ which(tmp$Cancer.gradeBin6to7 == 0 & tmp$Cancer.gradeBin8to9 == 0), ])
grade6.7.samples <- rownames(tmp[ which(tmp$Cancer.gradeBin6to7 == 1), ])
grade8.9.samples <- rownames(tmp[ which(tmp$Cancer.gradeBin8to9 == 1), ])
grade6.samples <- rownames(tmp[ which(tmp$Cancer.grade6 == 1), ])
grade7.samples <- rownames(tmp[ which(tmp$Cancer.grade7 == 1), ])
grade8.samples <- rownames(tmp[ which(tmp$Cancer.grade8 == 1), ])
grade9.samples <- rownames(tmp[ which(tmp$Cancer.grade9 == 1), ])

##### DIVIDE SAMPLES BY CANCER STAGE #####
stageII.samples <- rownames(tmp[ which(tmp$Cancer.stageIII == 0 & tmp$Cancer.stageIV == 0), ])
stageIII.IV.samples <- rownames(tmp[ which(tmp$Cancer.stageIII == 1 | tmp$Cancer.stageIV == 1), ])
#stageIV.samples <- rownames(tmp[ which(tmp$Cancer.stageIV == 1), ])
rm(tmp); gc()

##### AGES #####
age.LumA    <- as.numeric(as.vector(cov["Age", LuminalA.samples]))
age.LumB    <- as.numeric(as.vector(cov["Age", LuminalB.samples]))
age.TrpN    <- as.numeric(as.vector(cov["Age", Triplenegative.samples]))

age.Tpp     <- as.numeric(as.vector(cov["Age", ERPRposHer2pos.samples]))
age.Tpn     <- as.numeric(as.vector(cov["Age", ERPRposHer2neg.samples]))
age.Tnp     <- as.numeric(as.vector(cov["Age", ERPRnegHer2pos.samples]))
age.Tnn     <- as.numeric(as.vector(cov["Age", ERPRnegHer2neg.samples]))

age.H2p     <- as.numeric(as.vector(cov["Age", Her2p.samples]))
age.EPpHn   <- as.numeric(as.vector(cov["Age", ERPRpHer2n.samples]))
age.EPnHn   <- as.numeric(as.vector(cov["Age", ERPRnHer2n.samples]))

age.White.H2p     <- as.numeric(as.vector(cov["Age", Her2p.White.samples]))
age.White.EPpHn   <- as.numeric(as.vector(cov["Age", ERPRpHer2n.White.samples]))
age.White.EPnHn   <- as.numeric(as.vector(cov["Age", ERPRnHer2n.White.samples]))

age.Black.H2p     <- as.numeric(as.vector(cov["Age", Her2p.Black.samples]))
age.Black.EPpHn   <- as.numeric(as.vector(cov["Age", ERPRpHer2n.Black.samples]))
age.Black.EPnHn   <- as.numeric(as.vector(cov["Age", ERPRnHer2n.Black.samples]))

age.g3.5    <- as.numeric(as.vector(cov["Age", grade3.5.samples]))
age.g6.7    <- as.numeric(as.vector(cov["Age", grade6.7.samples]))
age.g8.9    <- as.numeric(as.vector(cov["Age", grade8.9.samples]))
age.g6      <- as.numeric(as.vector(cov["Age", grade6.samples]))
age.g7      <- as.numeric(as.vector(cov["Age", grade7.samples]))
age.g8      <- as.numeric(as.vector(cov["Age", grade8.samples]))
age.g9      <- as.numeric(as.vector(cov["Age", grade9.samples]))

age.sII     <- as.numeric(as.vector(cov["Age", stageII.samples]))
age.sIII.IV <- as.numeric(as.vector(cov["Age", stageIII.IV.samples]))

##### DNAm AGES #####
result.LumA    <- as.numeric(as.vector(cov["DNAm Age", LuminalA.samples]))
result.LumB    <- as.numeric(as.vector(cov["DNAm Age", LuminalB.samples]))
result.TrpN    <- as.numeric(as.vector(cov["DNAm Age", Triplenegative.samples]))

result.Tpp     <- as.numeric(as.vector(cov["DNAm Age", ERPRposHer2pos.samples]))
result.Tpn     <- as.numeric(as.vector(cov["DNAm Age", ERPRposHer2neg.samples]))
result.Tnp     <- as.numeric(as.vector(cov["DNAm Age", ERPRnegHer2pos.samples]))
result.Tnn     <- as.numeric(as.vector(cov["DNAm Age", ERPRnegHer2neg.samples]))

result.H2p     <- as.numeric(as.vector(cov["DNAm Age", Her2p.samples]))
result.EPpHn   <- as.numeric(as.vector(cov["DNAm Age", ERPRpHer2n.samples]))
result.EPnHn   <- as.numeric(as.vector(cov["DNAm Age", ERPRnHer2n.samples]))

result.White.H2p     <- as.numeric(as.vector(cov["DNAm Age", Her2p.White.samples]))
result.White.EPpHn   <- as.numeric(as.vector(cov["DNAm Age", ERPRpHer2n.White.samples]))
result.White.EPnHn   <- as.numeric(as.vector(cov["DNAm Age", ERPRnHer2n.White.samples]))

result.Black.H2p     <- as.numeric(as.vector(cov["DNAm Age", Her2p.Black.samples]))
result.Black.EPpHn   <- as.numeric(as.vector(cov["DNAm Age", ERPRpHer2n.Black.samples]))
result.Black.EPnHn   <- as.numeric(as.vector(cov["DNAm Age", ERPRnHer2n.Black.samples]))

result.g3.5    <- as.numeric(as.vector(cov["DNAm Age", grade3.5.samples]))
result.g6.7    <- as.numeric(as.vector(cov["DNAm Age", grade6.7.samples]))
result.g8.9    <- as.numeric(as.vector(cov["DNAm Age", grade8.9.samples]))
result.g6      <- as.numeric(as.vector(cov["DNAm Age", grade6.samples]))
result.g7      <- as.numeric(as.vector(cov["DNAm Age", grade7.samples]))
result.g8      <- as.numeric(as.vector(cov["DNAm Age", grade8.samples]))
result.g9      <- as.numeric(as.vector(cov["DNAm Age", grade9.samples]))

result.sII     <- as.numeric(as.vector(cov["DNAm Age", stageII.samples]))
result.sIII.IV <- as.numeric(as.vector(cov["DNAm Age", stageIII.IV.samples]))


##### RESIDUALS #####
residual.LumA    <- as.numeric(as.vector(cov["DNAm Age Residual", LuminalA.samples]))
residual.LumB    <- as.numeric(as.vector(cov["DNAm Age Residual", LuminalB.samples]))
residual.TrpN    <- as.numeric(as.vector(cov["DNAm Age Residual", Triplenegative.samples]))

residual.Tpp     <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRposHer2pos.samples]))
residual.Tpn     <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRposHer2neg.samples]))
residual.Tnp     <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnegHer2pos.samples]))
residual.Tnn     <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnegHer2neg.samples]))

residual.H2p     <- as.numeric(as.vector(cov["DNAm Age Residual", Her2p.samples]))
residual.EPpHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRpHer2n.samples]))
residual.EPnHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnHer2n.samples]))

residual.White.H2p     <- as.numeric(as.vector(cov["DNAm Age Residual", Her2p.White.samples]))
residual.White.EPpHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRpHer2n.White.samples]))
residual.White.EPnHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnHer2n.White.samples]))

residual.Black.H2p     <- as.numeric(as.vector(cov["DNAm Age Residual", Her2p.Black.samples]))
residual.Black.EPpHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRpHer2n.Black.samples]))
residual.Black.EPnHn   <- as.numeric(as.vector(cov["DNAm Age Residual", ERPRnHer2n.Black.samples]))

residual.g3.5    <- as.numeric(as.vector(cov["DNAm Age Residual", grade3.5.samples]))
residual.g6.7    <- as.numeric(as.vector(cov["DNAm Age Residual", grade6.7.samples]))
residual.g8.9    <- as.numeric(as.vector(cov["DNAm Age Residual", grade8.9.samples]))
residual.g6      <- as.numeric(as.vector(cov["DNAm Age Residual", grade6.samples]))
residual.g7      <- as.numeric(as.vector(cov["DNAm Age Residual", grade7.samples]))
residual.g8      <- as.numeric(as.vector(cov["DNAm Age Residual", grade8.samples]))
residual.g9      <- as.numeric(as.vector(cov["DNAm Age Residual", grade9.samples]))

residual.sII     <- as.numeric(as.vector(cov["DNAm Age Residual", stageII.samples]))
residual.sIII.IV <- as.numeric(as.vector(cov["DNAm Age Residual", stageIII.IV.samples]))
#residual.sIV  <- as.numeric(as.vector(cov["DNAm Age Residual", stageIV.samples]))


##########################################################################
# MERGED SUBTYPE RESIDUAL PLOT
##########################################################################
tmp.H2p <- data.frame(age.H2p, result.H2p, residual.H2p)
tmp.H2p$type <- "Her2+"
colnames(tmp.H2p) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.EPpHn <- data.frame(age.EPpHn, result.EPpHn, residual.EPpHn)
tmp.EPpHn$type <- "ER/PR+ Her2-"
colnames(tmp.EPpHn) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.EPnHn <- data.frame(age.EPnHn, result.EPnHn, residual.EPnHn)
tmp.EPnHn$type <- "ER- PR- Her2-"
colnames(tmp.EPnHn) <- c("Chrono.age", "DNAm.age", "res", "ttype")

df.ABN <- rbind(tmp.H2p, tmp.EPpHn, tmp.EPnHn)
rm(tmp.H2p); rm(tmp.EPpHn); rm(tmp.EPnHn); gc();
df.ABN$ttype <- factor(df.ABN$ttype, levels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"))

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.ABN, 
  bw = 25, 
  legname = "Cancer Subtype", 
  linetypes = c("solid", "dashed", "twodash"), 
  colors = c("deepskyblue4","darkorange3","seagreen4"), 
  labels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"), 
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Cancer Subtypes",
  annot.x = 85, 
  annot.y = 0.025,
  annot.sep = 0.002, 
  x.min = -45, 
  x.max = 125,
  y.min = 0, 
  y.max = 0.035, 
  leg.x = 0.8, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/subtype_histogram.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.H2p
r.list[[2]] <- residual.EPpHn
r.list[[3]] <- residual.EPnHn
p <- accel.box.plot(
  df.ABN, 
  residuals = r.list,
  x.label = "Subtype", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Acceleration for Cancer Subtypes",
  width = 0.6,
  leg.x = 0.74,
  leg.y = 0.92
); p

png( paste(model.dir, "CancerStudies/subtype_boxplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.ABN, 
  legname = "Cancer Subtype", 
  colors = c("deepskyblue4","darkorange3","seagreen4"), 
  labels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"), 
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.2, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/subtype_scatterplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED SUBTYPE RESIDUAL PLOT (White/Black)
##########################################################################
tmp.H2p <- data.frame(age.White.H2p, result.White.H2p, residual.White.H2p)
tmp.H2p$type <- "Her2+"
colnames(tmp.H2p) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.EPpHn <- data.frame(age.White.EPpHn, result.White.EPpHn, residual.White.EPpHn)
tmp.EPpHn$type <- "ER/PR+ Her2-"
colnames(tmp.EPpHn) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.EPnHn <- data.frame(age.White.EPnHn, result.White.EPnHn, residual.White.EPnHn)
tmp.EPnHn$type <- "ER- PR- Her2-"
colnames(tmp.EPnHn) <- c("Chrono.age", "DNAm.age", "res", "ttype")

df.ABN <- rbind(tmp.H2p, tmp.EPpHn, tmp.EPnHn)
rm(tmp.H2p); rm(tmp.EPpHn); rm(tmp.EPnHn); gc();
df.ABN$ttype <- factor(df.ABN$ttype, levels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"))

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.ABN, 
  bw = 25, 
  legname = "Cancer Subtype", 
  linetypes = c("solid", "dashed", "twodash"), 
  colors = c("deepskyblue4","darkorange3","seagreen4"), 
  labels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"), 
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Cancer Subtypes in White Samples",
  annot.x = 85, 
  annot.y = 0.025,
  annot.sep = 0.002, 
  x.min = -45, 
  x.max = 125,
  y.min = 0, 
  y.max = 0.035, 
  leg.x = 0.8, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/subtype_white_histogram.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.White.H2p
r.list[[2]] <- residual.White.EPpHn
r.list[[3]] <- residual.White.EPnHn
p <- accel.box.plot(
  df.ABN, 
  residuals = r.list,
  x.label = "Subtype", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Acceleration for Cancer Subtypes in White Samples",
  width = 0.6,
  leg.x = 0.74,
  leg.y = 0.92
); p

png( paste(model.dir, "CancerStudies/subtype_white_boxplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.ABN, 
  legname = "Cancer Subtype", 
  colors = c("deepskyblue4","darkorange3","seagreen4"), 
  labels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"), 
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age in White Samples",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.2, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/subtype_white_scatterplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED SUBTYPE RESIDUAL PLOT (White/Black)
##########################################################################
tmp.H2p <- data.frame(age.Black.H2p, result.Black.H2p, residual.Black.H2p)
tmp.H2p$type <- "Her2+"
colnames(tmp.H2p) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.EPpHn <- data.frame(age.Black.EPpHn, result.Black.EPpHn, residual.Black.EPpHn)
tmp.EPpHn$type <- "ER/PR+ Her2-"
colnames(tmp.EPpHn) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.EPnHn <- data.frame(age.Black.EPnHn, result.Black.EPnHn, residual.Black.EPnHn)
tmp.EPnHn$type <- "ER- PR- Her2-"
colnames(tmp.EPnHn) <- c("Chrono.age", "DNAm.age", "res", "ttype")

df.ABN <- rbind(tmp.H2p, tmp.EPpHn, tmp.EPnHn)
rm(tmp.H2p); rm(tmp.EPpHn); rm(tmp.EPnHn); gc();
df.ABN$ttype <- factor(df.ABN$ttype, levels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"))

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.ABN, 
  bw = 25, 
  legname = "Cancer Subtype", 
  linetypes = c("solid", "dashed", "twodash"), 
  colors = c("deepskyblue4","darkorange3","seagreen4"), 
  labels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"), 
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Cancer Subtypes in Black Samples",
  annot.x = 85, 
  annot.y = 0.025,
  annot.sep = 0.002, 
  x.min = -45, 
  x.max = 125,
  y.min = 0, 
  y.max = 0.035, 
  leg.x = 0.8, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/subtype_black_histogram.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.Black.H2p
r.list[[2]] <- residual.Black.EPpHn
r.list[[3]] <- residual.Black.EPnHn
p <- accel.box.plot(
  df.ABN, 
  residuals = r.list,
  x.label = "Subtype", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Acceleration for Cancer Subtypes in Black Samples",
  width = 0.6,
  leg.x = 0.74,
  leg.y = 0.92
); p

png( paste(model.dir, "CancerStudies/subtype_black_boxplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.ABN, 
  legname = "Cancer Subtype", 
  colors = c("deepskyblue4","darkorange3","seagreen4"), 
  labels = c("Her2+", "ER/PR+ Her2-", "ER- PR- Her2-"), 
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age in Black Samples",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.2, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/subtype_black_scatterplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED GRADE PLOT
##########################################################################
tmp.g3.5 <- data.frame(age.g3.5, result.g3.5, residual.g3.5)
tmp.g3.5$type <- "Grade 3-5"
colnames(tmp.g3.5) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.g6.7 <- data.frame(age.g6.7, result.g6.7, residual.g6.7)
tmp.g6.7$type <- "Grade 6-7"
colnames(tmp.g6.7) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.g8.9 <- data.frame(age.g8.9, result.g8.9, residual.g8.9)
tmp.g8.9$type <- "Grade 8-9"
colnames(tmp.g8.9) <- c("Chrono.age", "DNAm.age", "res", "ttype")

df.ABN <- rbind(tmp.g3.5, tmp.g6.7, tmp.g8.9)
rm(tmp.g3.5); rm(tmp.g6.7); rm(tmp.g8.9); gc();

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.ABN, 
  bw = 20, 
  legname = "Tumor Grade", 
  linetypes = c("solid", "dashed", "twodash"), 
  colors = c("gray16", "pink4", "paleturquoise4"), 
  labels = c("Grade 3-5", "Grade 6-7", "Grade 8-9"), 
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Tumor Grades",
  annot.x = 90, 
  annot.y = 0.032,
  annot.sep = 0.002, 
  x.min = -40, 
  x.max = 125,
  y.min = 0, 
  y.max = 0.045, 
  leg.x = 0.8, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/grade_histogram.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.g3.5
r.list[[2]] <- residual.g6.7
r.list[[3]] <- residual.g8.9
p <- accel.box.plot(
  df.ABN, 
  residuals = r.list, 
  width = 0.6,
  x.label = "Grade", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Acceleration for Tumor Grades"
); p

png( paste(model.dir, "CancerStudies/grade_boxplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.ABN, 
  legname = "Tumor Grade", 
  colors = c("gray16", "pink4", "paleturquoise4"), 
  labels = c("Grade 3-5", "Grade 6-7", "Grade 8-9"), 
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.2, 
  leg.y = 0.87
); p

png( paste(model.dir, "CancerStudies/grade_scatterplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED GRADE PLOT VERSION 2
##########################################################################
tmp.g3.5 <- data.frame(age.g3.5, result.g3.5, residual.g3.5)
tmp.g3.5$type <- "Grade 3-5"
colnames(tmp.g3.5) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.g6 <- data.frame(age.g6, result.g6, residual.g6)
tmp.g6$type <- "Grade 6"
colnames(tmp.g6) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.g7 <- data.frame(age.g7, result.g7, residual.g7)
tmp.g7$type <- "Grade 7"
colnames(tmp.g7) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.g8 <- data.frame(age.g8, result.g8, residual.g8)
tmp.g8$type <- "Grade 8"
colnames(tmp.g8) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.g9 <- data.frame(age.g9, result.g9, residual.g9)
tmp.g9$type <- "Grade 9"
colnames(tmp.g9) <- c("Chrono.age", "DNAm.age", "res", "ttype")

df.ABN <- rbind(tmp.g3.5, tmp.g6, tmp.g7, tmp.g8, tmp.g9)
rm(tmp.g3.5); rm(tmp.g6); rm(tmp.g7); rm(tmp.g8); rm(tmp.g9); gc();

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.ABN, 
  bw = 20, 
  legname = "Tumor Grade", 
  linetypes = c("solid", "dashed", "twodash", "longdash", "dotted"), 
  colors = c("gray16", "pink4", "paleturquoise4", "red", "blue"), 
  labels = c("Grade 3-5", "Grade 6", "Grade 7", "Grade 8", "Grade 9"), 
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Tumor Grades",
  annot.x = 90, 
  annot.y = 0.03,
  annot.sep = 0.002, 
  x.min = -40, 
  x.max = 125,
  y.min = 0, 
  y.max = 0.045, 
  leg.x = 0.8, 
  leg.y = 0.85
)

png( paste(model.dir, "CancerStudies/gradeV2_histogram.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.g3.5
r.list[[2]] <- residual.g6
r.list[[3]] <- residual.g7
r.list[[4]] <- residual.g8
r.list[[5]] <- residual.g9
p <- accel.box.plot(
  df.ABN, 
  residuals = r.list, 
  width = 0.6,
  x.label = "Grade", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Acceleration for Tumor Grades",
  leg.y = 0.9
)

png( paste(model.dir, "CancerStudies/gradeV2_boxplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.ABN, 
  legname = "Tumor Grade", 
  colors = c("gray16", "pink4", "paleturquoise4", "red", "blue"), 
  labels = c("Grade 3-5", "Grade 6", "Grade 7", "Grade 8", "Grade 9"), 
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.2, 
  leg.y = 0.85
); p

png( paste(model.dir, "CancerStudies/gradeV2_scatterplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


##########################################################################
# MERGED STAGE PLOT
##########################################################################
tmp.sII <- data.frame(age.sII, result.sII, residual.sII)
tmp.sII$type <- "Stage II"
colnames(tmp.sII) <- c("Chrono.age", "DNAm.age", "res", "ttype")

tmp.sIII.IV <- data.frame(age.sIII.IV, result.sIII.IV, residual.sIII.IV)
tmp.sIII.IV$type <- "Stage III & IV"
colnames(tmp.sIII.IV) <- c("Chrono.age", "DNAm.age", "res", "ttype")

#tmp.sIV <- data.frame(residual.sIV)
#tmp.sIV$type <- "Stage IV"
#colnames(tmp.sIV) <- c("res", "ttype")

df.ABN <- rbind(tmp.sII, tmp.sIII.IV)
rm(tmp.sII); rm(tmp.sIII.IV); gc();

##### HISTOGRAM #####
p <- accel.hist.plot(
  df.ABN, 
  bw = 20, 
  legname = "Tumor Stage", 
  linetypes = c("solid", "dashed"), 
  colors = c("darkslategrey","firebrick3"), 
  labels = c("Stage II", "Stage III & IV"), 
  x.label = "DNAm Age Acceleration [Years]", 
  y.label = "Frequency [Arbitrary Units]", 
  title = "DNAm Age Acceleration for Tumor Stages",
  annot.x = 85, 
  annot.y = 0.03,
  annot.sep = 0.002, 
  x.min = -40, 
  x.max = 125,
  y.min = 0, 
  y.max = 0.04, 
  leg.x = 0.8, 
  leg.y = 0.87
)

png( paste(model.dir, "CancerStudies/stage_histogram.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### BAR PLOT #####
r.list <- list()
r.list[[1]] <- residual.sII
r.list[[2]] <- residual.sIII.IV
p <- accel.box.plot(
  df.ABN, 
  residuals = r.list, 
  width = 0.6,
  x.label = "Stage", 
  y.label = "DNAm Age Acceleration [Years]",
  title = "DNAm Age Acceleration for Cancer Stages",
  leg.x = 0.75
)

png( paste(model.dir, "CancerStudies/stage_boxplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##### DNAm AGE vs CHRONOLOGICAL AGE #####
p <- DNAmAge.ChronoAge.plot(
  df.ABN, 
  legname = "Tumor Stage", 
  colors = c("darkslategrey","firebrick3"), 
  labels = c("Stage II", "Stage III & IV"),
  x.label = "Chronological Age [Years]", 
  y.label = "DNAm Age [Years]", 
  title = "DNAm Age vs Chronological Age",
  x.min = 0, 
  x.max = 90,
  y.min = 0, 
  y.max = 185, 
  leg.x = 0.2, 
  leg.y = 0.85
); p

png( paste(model.dir, "CancerStudies/stage_scatterplot.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()

##########################################################################
# EXPLORE STAGE-GRADE RELATION
##########################################################################
stages.with.grades <- c(4, 2, 4, 2, 2, 3, 2, 4, 2, 3, 2, 3, 2, 2, 2, 4, 4)
grades.with.stages <- c(8, 9, 7, 7, 6, 9, 6, 9, 7, 8, 6, 8, 6, 7, 9, 9, 8)

png( paste(model.dir, "CancerStudies/grade-stage.png", sep = ''), width = 500, height = 500, units = "px" )
plot(stages.with.grades, 
     grades.with.stages, 
     main="Tumor Grade-Stage Association", 
     xlab="Stage", 
     ylab="Grade", 
     pch=19,
     xlim=c(0,5), 
     xaxs="i",
     ylim=c(0,10), 
     yaxs="i",
     tck = 0.02
) 
abline(lm(grades.with.stages~stages.with.grades), col="black") # regression line (y~x) 
abline(a=0, b=1, col = "red", lty = 2)
r <- sqrt(summary(lm(grades.with.stages~stages.with.grades))$r.squared)
text(0.5,9.5, paste("r = ", round(r, digits = 4), sep = ""))
dev.off()


##########################################################################
# STATISTICAL TESTS
##########################################################################

##### ALLOW OUTLIERS #####
wilcox.test(residual.LumA,    mu=0, conf.int = TRUE)
wilcox.test(residual.LumB,    mu=0, conf.int = TRUE)
wilcox.test(residual.TrpN,    mu=0, conf.int = TRUE)
wilcox.test(residual.Tpp,     mu=0, conf.int = TRUE)
wilcox.test(residual.Tpn,     mu=0, conf.int = TRUE)
wilcox.test(residual.Tnp,     mu=0, conf.int = TRUE)
wilcox.test(residual.Tnn,     mu=0, conf.int = TRUE)
wilcox.test(residual.H2p,     mu=0, conf.int = TRUE)
wilcox.test(residual.EPpHn,   mu=0, conf.int = TRUE)
wilcox.test(residual.EPnHn,   mu=0, conf.int = TRUE)
wilcox.test(residual.White.H2p,     mu=0, conf.int = TRUE)
wilcox.test(residual.White.EPpHn,   mu=0, conf.int = TRUE)
wilcox.test(residual.White.EPnHn,   mu=0, conf.int = TRUE)
wilcox.test(residual.Black.H2p,     mu=0, conf.int = TRUE)
wilcox.test(residual.Black.EPpHn,   mu=0, conf.int = TRUE)
wilcox.test(residual.Black.EPnHn,   mu=0, conf.int = TRUE)
wilcox.test(residual.g3.5,    mu=0, conf.int = TRUE)
wilcox.test(residual.g6.7,    mu=0, conf.int = TRUE)
wilcox.test(residual.g8.9,    mu=0, conf.int = TRUE)
wilcox.test(residual.g6,      mu=0, conf.int = TRUE)
wilcox.test(residual.g7,      mu=0, conf.int = TRUE)
wilcox.test(residual.g8,      mu=0, conf.int = TRUE)
wilcox.test(residual.g9,      mu=0, conf.int = TRUE)
wilcox.test(residual.sII,     mu=0, conf.int = TRUE)
wilcox.test(residual.sIII.IV, mu=0, conf.int = TRUE)


#wilcox.test(residual.sIV,  mu=0, conf.int = TRUE)
median(residual.sIII.IV)
mean(residual.sII[ !(residual.sII %in% boxplot.stats(residual.sII)$out) ] )

##### REMOVE OUTLIERS #####
wilcox.test(residual.LumA[ !(residual.LumA %in% boxplot.stats(residual.LumA)$out) ],          mu=0, conf.int = TRUE)
wilcox.test(residual.LumB[ !(residual.LumB %in% boxplot.stats(residual.LumB)$out) ],          mu=0, conf.int = TRUE)
wilcox.test(residual.TrpN[ !(residual.TrpN %in% boxplot.stats(residual.TrpN)$out) ],          mu=0, conf.int = TRUE)
wilcox.test(residual.Tpp[ !(residual.Tpp %in% boxplot.stats(residual.Tpp)$out) ],             mu=0, conf.int = TRUE)
wilcox.test(residual.Tpn[ !(residual.Tpn %in% boxplot.stats(residual.Tpn)$out) ],             mu=0, conf.int = TRUE)
wilcox.test(residual.Tnp[ !(residual.Tnp %in% boxplot.stats(residual.Tnp)$out) ],             mu=0, conf.int = TRUE)
wilcox.test(residual.Tnn[ !(residual.Tnn %in% boxplot.stats(residual.Tnn)$out) ],             mu=0, conf.int = TRUE)
wilcox.test(residual.H2p[ !(residual.H2p %in% boxplot.stats(residual.H2p)$out) ],             mu=0, conf.int = TRUE)
wilcox.test(residual.EPpHn[ !(residual.EPpHn %in% boxplot.stats(residual.EPpHn)$out) ],       mu=0, conf.int = TRUE)
wilcox.test(residual.EPnHn[ !(residual.EPnHn %in% boxplot.stats(residual.EPnHn)$out) ],       mu=0, conf.int = TRUE)
wilcox.test(residual.g3.5[ !(residual.g3.5 %in% boxplot.stats(residual.g3.5)$out) ],          mu=0, conf.int = TRUE)
wilcox.test(residual.g6.7[ !(residual.g6.7 %in% boxplot.stats(residual.g6.7)$out) ],          mu=0, conf.int = TRUE)
wilcox.test(residual.g8.9[ !(residual.g8.9 %in% boxplot.stats(residual.g8.9)$out) ],          mu=0, conf.int = TRUE)
wilcox.test(residual.g6[ !(residual.g6 %in% boxplot.stats(residual.g6)$out) ],                mu=0, conf.int = TRUE)
wilcox.test(residual.g7[ !(residual.g7 %in% boxplot.stats(residual.g7)$out) ],                mu=0, conf.int = TRUE)
wilcox.test(residual.g8[ !(residual.g8 %in% boxplot.stats(residual.g8)$out) ],                mu=0, conf.int = TRUE)
wilcox.test(residual.g9[ !(residual.g9 %in% boxplot.stats(residual.g9)$out) ],                mu=0, conf.int = TRUE)
wilcox.test(residual.sII[ !(residual.sII %in% boxplot.stats(residual.sII)$out) ],             mu=0, conf.int = TRUE)
wilcox.test(residual.sIII.IV[ !(residual.sIII.IV %in% boxplot.stats(residual.sIII.IV)$out) ], mu=0, conf.int = TRUE)
