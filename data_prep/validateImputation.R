rm(list=ls()); gc();
setwd( '/Users/jrca253/Documents/EpigeneticAge/test_code/ImputeCheck_850k')
library(ggplot2)
library(matrixStats)

meth.file.850       <- "Meth_850k_HorvathClock.csv"
meth.file.K.truseq  <- "meth_K_CGNumber_HorvathClockCpGs_Imputed.csv"
meth.file.T.truseq  <- "meth_T_CGNumber_HorvathClockCpGs_Imputed.csv"
imputed.cpgs.file.K <- "ImputedHorvathClockCpGs_K.txt"
imputed.cpgs.file.T <- "ImputedHorvathClockCpGs_T.txt"

###########################################################################################
# LOAD/PREP DATA
###########################################################################################

##### IMPUTED CpG List #####
imputed.cpgs.K <- read.table(imputed.cpgs.file.K)
imputed.cpgs.K <- as.character(as.vector(imputed.cpgs.K[,1]))

imputed.cpgs.T <- read.table(imputed.cpgs.file.T)
imputed.cpgs.T <- as.character(as.vector(imputed.cpgs.T[,1]))

##### 850 #####
df.850   <- read.csv(meth.file.850, header = TRUE)
df.850$X <- as.character(df.850$X)
df.850.K <- df.850[, c("X", "K104496", "K100875", "K102645")]
df.850.T <- df.850[, c("X", "T060", "T040", "T004", "T042")]
rm(df.850); gc()

df.850.K <- df.850.K[which( df.850.K$X %in% imputed.cpgs.K),]
df.850.T <- df.850.T[which( df.850.T$X %in% imputed.cpgs.T),]

##### TruSeq #####
df.truseq.K <- read.csv(meth.file.K.truseq, header = TRUE)
df.truseq.K$position <- as.character(df.truseq.K$position)
df.truseq.K <- df.truseq.K[which( df.truseq.K$position %in% imputed.cpgs.K),]
df.truseq.K <- df.truseq.K[which( df.truseq.K$position %in% df.850.K$X),]

df.truseq.T <- read.csv(meth.file.T.truseq, header = TRUE)
df.truseq.T$position <- as.character(df.truseq.T$position)
df.truseq.T <- df.truseq.T[which( df.truseq.T$position %in% imputed.cpgs.T),]
df.truseq.T <- df.truseq.T[which( df.truseq.T$position %in% df.850.T$X),]

##### ORDER ROWS #####
df.850.K   <- df.850.K[ order(df.850.K$X),]
df.850.T   <- df.850.T[ order(df.850.T$X),]

rownames(df.850.K) <- as.character(as.vector(df.850.K$X))
rownames(df.850.T) <- as.character(as.vector(df.850.T$X))

df.truseq.K <- df.truseq.K[ order(df.truseq.K$position),]
df.truseq.T <- df.truseq.T[ order(df.truseq.T$position),]

rownames(df.truseq.K) <- as.character(as.vector(df.truseq.K$position))
rownames(df.truseq.T) <- as.character(as.vector(df.truseq.T$position))

df.850.K$X          <- NULL
df.850.T$X          <- NULL
df.truseq.K$position <- NULL
df.truseq.T$position <- NULL


###########################################################################################
# COMPARE METHYLATION VALUES PER IMPUTED CpG
###########################################################################################

##### MEDIAN #####
df.850.K.Median   <- rowMedians(as.matrix(df.850.K))
df.850.T.Median   <- rowMedians(as.matrix(df.850.T))
df.truseq.K.Median <- rowMedians(as.matrix(df.truseq.K))
df.truseq.T.Median <- rowMedians(as.matrix(df.truseq.T))

#### ADD TO DATA.FRAMES #####
df.850.K$Median   <- df.850.K.Median
df.850.T$Median   <- df.850.T.Median
df.truseq.K$Median <- df.truseq.K.Median
df.truseq.T$Median <- df.truseq.T.Median

###########################################################################################
# MEAN PLOT
###########################################################################################
truseq.850.median.ratio.K <- df.truseq.K$Median / df.850.K$Median
ttype <- rep("K", length(truseq.850.median.ratio.K))
df.compare.medians.K <- data.frame(truseq.850.median.ratio.K, ttype)
colnames(df.compare.medians.K) <- c("Ratio", "ttype")

truseq.850.median.ratio.T <- df.truseq.T$Median / df.850.T$Median
ttype <- rep("T", length(truseq.850.median.ratio.T))
df.compare.medians.T <- data.frame(truseq.850.median.ratio.T, ttype)
colnames(df.compare.medians.T) <- c("Ratio", "ttype")

df.compare <- rbind(df.compare.medians.K, df.compare.medians.T)
rm(df.compare.medians.K); rm(df.compare.medians.T); gc()

p <- ggplot(df.compare, aes(x=ttype, y=Ratio)) +
  geom_boxplot(aes(x=ttype, y=Ratio), data = df.compare, width = 0.5) +
  geom_hline(yintercept = 1, color = "gray", size = 1, linetype = "dashed") + 
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 5)
  ) + 
  annotate(
    "text", x = 1, y = median(truseq.850.median.ratio.K)+0.1, 
    label = paste("Median = ", round(median(truseq.850.median.ratio.K), digits = 2), sep = ""),
    color = "blue"
  ) + 
  annotate(
    "text", x = 2, y = median(truseq.850.median.ratio.T)+0.1, 
    label = paste("Median = ", round(median(truseq.850.median.ratio.T), digits = 2), sep = ""),
    color = "blue"
  ) + 
  labs(x = "Tissue Type") + 
  labs(y = "(Median Imputed TruSeq beta) / (Median 850 beta)") + 
  labs(title = "Methylation Values Per CpG Comparison") + 
  theme_bw(base_size = 15) +
  theme(
    legend.position=c(0.2, 0.85),
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png("ValidateImputation.png", width = 500, height = 500, units = 'px')
p
dev.off()