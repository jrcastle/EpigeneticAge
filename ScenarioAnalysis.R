setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)

rm(list=ls()); gc();

model.dir   <- "HorvathClock/"
result.file <- "ScenarioResults.csv"
cov.file.K  <- "cov_K.txt"


###########################################################################################
# LOAD DNAm AGES AND SAMPLE AGES
###########################################################################################

##### DNAm AGE SCENARIOS #####
scenarios <- read.csv( paste(model.dir, result.file, sep = '') )
DNAmAge.21k.Imputed.Normalized   <- as.numeric( scenarios$DNAmAge.21k.Imputed.Normalized )
DNAmAge.21k.Imputed.Unnormalized <- as.numeric( scenarios$DNAmAge.21k.Imputed.Unnormalized )
DNAmAge.HorvathClock.Raw         <- as.numeric( scenarios$DNAmAge.HorvathClock.Raw )

##### SAMPLE AGES #####
sample.ages <- read.table(paste(model.dir, cov.file.K, sep = ''), header = TRUE, row.names = 1, sep = '\t')
sample.ages <- as.numeric(as.vector(sample.ages["Age",]))


###########################################################################################
# RESIDUALS
###########################################################################################
Res.DNAmAge.21k.Imputed.Normalized   <- DNAmAge.21k.Imputed.Normalized   - sample.ages
Res.DNAmAge.21k.Imputed.Unnormalized <- DNAmAge.21k.Imputed.Unnormalized - sample.ages
Res.DNAmAge.HorvathClock.Raw         <- DNAmAge.HorvathClock.Raw         - sample.ages


###########################################################################################
# PLOT
###########################################################################################
tmp.1 <- data.frame(Res.DNAmAge.21k.Imputed.Normalized)
tmp.1$Scenario = "s1"
colnames(tmp.1) <- c("res", "Scenario")

tmp.2 <- data.frame(Res.DNAmAge.21k.Imputed.Unnormalized)
tmp.2$Scenario = "s2"
colnames(tmp.2) <- c("res", "Scenario")

tmp.3 <- data.frame(Res.DNAmAge.HorvathClock.Raw)
tmp.3$Scenario = "s3"
colnames(tmp.3) <- c("res", "Scenario")

df.123 <- rbind(tmp.1, tmp.2, tmp.3)
rm(tmp.1)
rm(tmp.2)
rm(tmp.3)
gc()

median.error.1 <- median(Res.DNAmAge.21k.Imputed.Normalized)
stddev.error.1 <- sd(Res.DNAmAge.21k.Imputed.Normalized)
median.error.2 <- median(Res.DNAmAge.21k.Imputed.Unnormalized)
stddev.error.2 <- sd(Res.DNAmAge.21k.Imputed.Unnormalized)
median.error.3 <- median(Res.DNAmAge.HorvathClock.Raw)
stddev.error.3 <- sd(Res.DNAmAge.HorvathClock.Raw)

#df.123 <- data.frame(Res.DNAmAge.21k.Imputed.Normalized, Res.DNAmAge.21k.Imputed.Unnormalized, Res.DNAmAge.HorvathClock.Raw)
#colnames(df.123) <- c("s1", "s2", "s3")

#p <- ggplot(df.123, aes(x = res)) +
  #geom_histogram(binwidth = 5, alpha = .5) +
  #scale_fill_manual(name = "ttype", values = c("red","blue"), labels=c("K", "T")) + 
#  geom_freqpoly(data=subset(df.123, Scenario == "s1"),  binwidth = 5, alpha = 1) +
#  geom_freqpoly(data=subset(df.123, Scenario == "s3"),  binwidth = 5, alpha = 0.8) +
#  geom_freqpoly(data=subset(df.123, Scenario == "s2"),  binwidth = 5, alpha = 0.8) +
#    scale_fill_manual(
#    name = "Scenario", 
#    values = c("darkorange3", "black", "deepskyblue4"), 
#    labels=c("1) 21k.Imputed.Normalized", "2) 21k.Imputed.Unnormalized", "3) HorvathClock.Raw")
#  ) +
p <- ggplot(df.123, aes(res, color = Scenario)) +
  geom_freqpoly(size = 1.2, binwidth = 5) + 
  scale_color_manual(
    name = "Scenario", 
    values = c("darkorange3", "black", "deepskyblue4"), 
    #values = c("red", "black", "blue")
    labels=c("1) 21k.Imputed.Normalized", "2) 21k.Imputed.Unnormalized", "3) HorvathClock.Raw")
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 100)
  ) + 
  labs(x = "DNAm Age Acceleration") + 
  labs(y = "Frequency") + 
  labs(title = "Scenario Results") + 
  annotate(
    "text", x = 50, y = 65, 
    label = paste("Med Accel (1) = ", round(median.error.1, digits = 1), sep = ""), 
    color = "darkorange3"
  ) + 
  annotate(
    "text", x = 50, y = 60, 
    label = paste("SD Accel (1) = ", round(stddev.error.1, digits = 1), sep = ""), 
    color = "darkorange3"
  ) +
  annotate(
    "text", x = 50, y = 55, 
    label = paste("Med Accel (2) = ", round(median.error.2, digits = 1), sep = ""), 
    color = "black"
  ) + 
  annotate(
    "text", x = 50, y = 50, 
    label = paste("SD Accel (2) = ", round(stddev.error.2, digits = 1), sep = ""), 
    color = "black"
  ) +
  annotate(
    "text", x = 50, y = 45, 
    label = paste("Med Accel (3) = ", round(median.error.3, digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) + 
  annotate(
    "text", x = 50, y = 40, 
    label = paste("SD Accel (3) = ", round(stddev.error.3, digits = 1), sep = ""), 
    color = "deepskyblue4"
  ) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
  theme_bw() + 
  theme(legend.position=c(0.8, 0.85)) +
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  ); p

png( paste(model.dir, "ScenarioCompare.png", sep = ''), width = 500, height = 500, units = "px" )
p
dev.off()


