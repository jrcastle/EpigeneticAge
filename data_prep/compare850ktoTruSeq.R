rm(list=ls()); gc();
setwd( '/Users/jrca253/Documents/EpigeneticAge/test_code/ImputeCheck_850k')
library(ggplot2)
library(matrixStats)

meth.file.850k     <- "M-value_11_sample.csv"
meth.file.K.truseq <- "meth_K_850k_CGNumber.txt"
meth.file.T.truseq <- "meth_T_850k_CGNumber.txt"

###########################################################################################
# LOAD/PREP DATA
###########################################################################################

##### TruSeq #####
df.truseq.K <- read.csv(meth.file.K.truseq, header = TRUE)
df.truseq.T <- read.csv(meth.file.T.truseq, header = TRUE)

##### 850k #####
df.850k   <- read.csv(meth.file.850k, header = TRUE)
df.850k.K <- df.850k[, c("X", "K104496", "K100875", "K102645")]
df.850k.T <- df.850k[, c("X", "T060", "T040", "T004", "T042")]
rm(df.850k); gc()

##### COMPLETE CASES ONLY #####
df.truseq.K <- df.truseq.K[rowSums(is.na(df.truseq.K)) <= 0,]
df.truseq.T <- df.truseq.T[rowSums(is.na(df.truseq.T)) <= 0,]
df.850k.K   <- df.850k.K[rowSums(is.na(df.850k.K)) <= 0,]
df.850k.T   <- df.850k.T[rowSums(is.na(df.850k.T)) <= 0,]

##### OVERLAPS #####
df.truseq.K$position <- as.character(df.truseq.K$position)
df.truseq.T$position <- as.character(df.truseq.T$position)
df.850k.K$X          <- as.character(df.850k.K$X)
df.850k.T$X          <- as.character(df.850k.T$X)

df.truseq.K <- df.truseq.K[which(df.truseq.K$position %in% df.850k.K$X),]
df.truseq.T <- df.truseq.T[which(df.truseq.T$position %in% df.850k.T$X),]
df.850k.K   <- df.850k.K[which(df.850k.K$X %in% df.truseq.K$position),]
df.850k.T   <- df.850k.T[which(df.850k.T$X %in% df.truseq.T$position),]

##### ORDER ROWS #####
df.truseq.K<- df.truseq.K[order(df.truseq.K$position),]
df.850k.K   <- df.850k.K[order(df.850k.K$X),]
df.850k.K$X == df.truseq.K$position

df.truseq.T<- df.truseq.T[order(df.truseq.T$position),]
df.850k.T   <- df.850k.T[order(df.850k.T$X),]
df.850k.T$X == df.truseq.T$position

rownames(df.truseq.K) <- as.character(as.vector(df.truseq.K$position))
rownames(df.truseq.T) <- as.character(as.vector(df.truseq.T$position))
rownames(df.850k.K) <- as.character(as.vector(df.850k.K$X))
rownames(df.850k.T) <- as.character(as.vector(df.850k.T$X))

df.truseq.K$position <- NULL
df.truseq.T$position <- NULL
df.850k.K$X          <- NULL
df.850k.T$X          <- NULL

###########################################################################################
# COMPARE METHYLATION VALUES PER CpG
###########################################################################################

##### MEDIAN #####
df.truseq.K.Median <- rowMedians(as.matrix(df.truseq.K))
df.truseq.T.Median <- rowMedians(as.matrix(df.truseq.T))
df.850k.K.Median   <- rowMedians(as.matrix(df.850k.K))
df.850k.T.Median   <- rowMedians(as.matrix(df.850k.T))

#### ADD TO DATA.FRAMES #####
df.truseq.K$Median <- df.truseq.K.Median
df.truseq.T$Median <- df.truseq.T.Median
df.850k.K$Median   <- df.850k.K.Median
df.850k.T$Median   <- df.850k.T.Median


###########################################################################################
# MEAN PLOT
###########################################################################################
truseq.850K.median.ratio.K <- df.truseq.K$Median / df.850k.K$Median
ttype <- rep("K", length(truseq.850K.median.ratio.K))
df.compare.medians.K <- data.frame(truseq.850K.median.ratio.K, ttype)
colnames(df.compare.medians.K) <- c("Ratio", "ttype")

truseq.850K.median.ratio.T <- df.truseq.T$Median / df.850k.T$Median
ttype <- rep("T", length(truseq.850K.median.ratio.T))
df.compare.medians.T <- data.frame(truseq.850K.median.ratio.T, ttype)
colnames(df.compare.medians.T) <- c("Ratio", "ttype")

df.compare <- rbind(df.compare.medians.K, df.compare.medians.T)
rm(df.compare.medians.K); rm(df.compare.medians.T); gc()

p <- ggplot(df.compare, aes(x=ttype, y=Ratio)) +
  geom_boxplot(aes(x=ttype, y=Ratio), data = df.compare, width = 0.5) +
  geom_hline(yintercept = 1, color = "gray", size = 1, linetype = "dashed") + 
  scale_color_manual(
    name = "Measure", 
    values = c("red", "blue"), 
    labels = c("Mean", "Median")
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 2)
  ) + 
  annotate(
    "text", x = 1, y = median(truseq.850K.median.ratio.K)-0.05, 
    label = paste("Median = ", round(median(truseq.850K.median.ratio.K), digits = 2), sep = ""),
    color = "blue"
  ) + 
  annotate(
    "text", x = 2, y = median(truseq.850K.median.ratio.T)-0.1, 
    label = paste("Median = ", round(median(truseq.850K.median.ratio.T), digits = 2), sep = ""),
    color = "blue"
  ) + 
  labs(x = "Tissue Type") + 
  labs(y = "(Median TruSeq beta) / (Median 850k beta)") + 
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

png("Truseq_850k_compare.png",  width = 500, height = 500, units = 'px')
p
dev.off()
