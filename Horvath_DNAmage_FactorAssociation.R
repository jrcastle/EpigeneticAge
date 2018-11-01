setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(moments)

rm(list=ls()); gc();

model.dir     <- "HorvathClock/"
result.file.K <- paste(model.dir, "ScenarioResults.csv", sep = "")
cov.file.K    <- paste(model.dir, "cov_K.txt", sep = "")


###########################################################################################
# LOAD DNAm AGES AND SAMPLE AGES
###########################################################################################

##### DNAm AGES ####
scenarios <- read.csv( result.file.K )
result.K  <- as.numeric( scenarios$DNAmAge.21k.Imputed.Normalized )

##### SAMPLE AGES #####
df.cov.K      <- read.table(cov.file.K, header = TRUE, row.names = 1, sep = '\t')
sample.ages.K <- as.numeric(as.vector(df.cov.K["Age",]))
residual.K    <- result.K - sample.ages.K

##### Add result and residual to cov dataframe #####
df.cov.K["DNAm.Age",] <- result.K
df.cov.K["DNAm.Age.Residual",] <- residual.K

for (i in 2:dim(df.cov.K)[[2]] ){
  df.cov.K[,i]=as.numeric(as.character(df.cov.K[,i]))
}

##########################################################################
# Multivariate Correlation Analysis
##########################################################################

DNAm.Age          <- as.numeric(as.vector(df.cov.K["DNAm.Age",]))
DNAm.Age.Residual <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",]))
Race              <- as.numeric(as.vector(df.cov.K["RaceWhite",]))
BMI               <- as.numeric(as.vector(df.cov.K["LocationUrban",]))
Cig.Pack.Years    <- as.numeric(as.vector(df.cov.K["Cig.Pack.Years",]))
Smoking           <- as.numeric(as.vector(df.cov.K["SmokingYes",]))
Drinking          <- as.numeric(as.vector(df.cov.K["DrinkingYes",]))
Menarche          <- as.numeric(as.vector(df.cov.K["Menarche",]))
Been.PregnantYes  <- as.numeric(as.vector(df.cov.K["Been.PregnantYes",]))
Times.Pregnant    <- as.numeric(as.vector(df.cov.K["Times.Pregnant",]))
Parity            <- as.numeric(as.vector(df.cov.K["Parity",]))
Parity[Parity == 41] <- NA #Outlier in the data...
Age.FB            <- as.numeric(as.vector(df.cov.K["Age.FB",]))
Menopause.Age     <- as.numeric(as.vector(df.cov.K["Menopause.Age",]))
VDYes             <- as.numeric(as.vector(df.cov.K["VDYes",]))

Pre.Menopause  <- as.numeric(as.vector(df.cov.K["MenopausePre-menopausal",]))
Post.Menopause <- as.numeric(as.vector(df.cov.K["MenopausePost-menopausal",]))
Location.Urban <- as.numeric(as.vector(df.cov.K["LocationUrban",]))
Location.Rural <- as.numeric(as.vector(df.cov.K["LocationRural",]))

Menopause <- c()
Location <- c()

for(i in 1:length(Location.Rural) ){
  
  pre <- Pre.Menopause[i]
  post <- Post.Menopause[i]
  rur <- Location.Rural[i]
  urb <- Location.Urban[i]
  
  post.yes <- ifelse(post == 1 && pre == 0, 1, ifelse(post == 0 && pre == 1, 0, NA) )
  urban.yes <- ifelse(urb == 1 && rur == 0, 1, ifelse(urb == 0 && rur == 1, 0, NA) )
  
  Location <- as.vector(as.vector(c( Location, urban.yes )))
  Menopause <- as.vector(as.vector(c( Menopause, post.yes )))
}

Location  <- as.numeric(Location)
Menopause <- as.numeric(Menopause)

lm.DNAm.Age <- lm(DNAm.Age ~ 
                    Race + 
                    BMI + 
                    Cig.Pack.Years + 
                    Smoking + 
                    Drinking + 
                    Menarche + 
                    Been.PregnantYes + 
                    Times.Pregnant + 
                    Parity + 
                    Age.FB + 
                    #Menopause.Age + 
                    VDYes + 
                    Menopause + 
                    Location
)

lm.DNAm.Age.Residual <- lm(DNAm.Age.Residual ~ 
                             Race + 
                             BMI + 
                             Cig.Pack.Years + 
                             Smoking + 
                             Drinking + 
                             Menarche + 
                             Been.PregnantYes + 
                             Times.Pregnant + 
                             Parity + 
                             Age.FB + 
                             #Menopause.Age + 
                             VDYes + 
                             Menopause + 
                             Location
)


##########################################################################
# Univariate Analysis - Race
##########################################################################
df.univar <- data.frame()

##### Correlate DNAm Age and Residual to Race #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["RaceWhite",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["RaceWhite",]) ]))
Race                <- as.numeric(as.vector(df.cov.K["RaceWhite", !is.na(df.cov.K["RaceWhite",]) ]))

DNAm.Age.Rho  <- cor.test(Race, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Race, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Race, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Race, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Race", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Race", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Race", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Race", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
DNAm.Age.White            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 1)]))
DNAm.Age.AfrAmer          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 0)]))
DNAm.Age.Residual.White   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 1)]))
DNAm.Age.Residual.AfrAmer <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 0)]))

race                   <- c("African American", "White")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.AfrAmer), mean(DNAm.Age.White)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.AfrAmer), sd(DNAm.Age.White)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.AfrAmer), mean(DNAm.Age.Residual.White)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.AfrAmer), sd(DNAm.Age.Residual.White)) )

df.tmp <- data.frame(race, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=race, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Race") + 
  labs(y = "Mean DNAm Age") + 
  annotate("text", x = 0.75, y = 74.5, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 71.5, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.AfrAmer), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.White), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Race-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = race, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Race") + 
  labs(y = "Mean DNAm Age Acceleration") + 
  annotate("text", x = 0.75, y = 21, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 20, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.AfrAmer), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.White), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Race-DNAmAgeResidual.png", sep = '') )
p
dev.off()

##########################################################################
# Univariate Analysis - Location
##########################################################################

##### Plot Correlations #####
Location.Rural          <- as.numeric(as.vector(df.cov.K["LocationRural", (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
Loacation.Urban         <- as.numeric(as.vector(df.cov.K["LocationUrban", (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
DNAm.Age.Rural          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
DNAm.Age.Urban          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
DNAm.Age.Residual.Rural <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
DNAm.Age.Residual.Urban <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))

location               <- c("Rural", "Urban")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.Rural), mean(DNAm.Age.Urban)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.Rural), sd(DNAm.Age.Urban)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.Rural), mean(DNAm.Age.Residual.Urban)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.Rural), sd(DNAm.Age.Residual.Urban)) )

##### Correlate DNAm Age and Residual to Location #####
Location.Rural[Location.Rural == 1] <- 0
DNAm.Age            <- as.numeric( c(DNAm.Age.Rural, DNAm.Age.Urban) )
DNAm.Age.Residual   <- as.numeric( c(DNAm.Age.Residual.Rural, DNAm.Age.Residual.Urban) )
Location            <- as.numeric( c(Location.Rural, Loacation.Urban) )

DNAm.Age.Rho  <- cor.test(Location, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Location, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Location, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Location, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Location", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Location", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Location", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Location", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

df.tmp <- data.frame(location, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=location, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Location") + 
  labs(y = "Mean DNAm Age") + 
  annotate("text", x = 0.75, y = 68, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 65, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.Rural), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.Urban), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Location-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = location, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Location") + 
  labs(y = "Mean DNAm Age Acceleration") + 
  annotate("text", x = 2.25, y = 23, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 2.25, y = 22, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Rural), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.Urban), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Location-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - BMI
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["BMI",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["BMI",]) ]))
BMI                 <- as.numeric(as.vector(df.cov.K["BMI", !is.na(df.cov.K["BMI",]) ]))

DNAm.Age.Rho  <- cor.test(BMI, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(BMI, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(BMI, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(BMI, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["BMI", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["BMI", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["BMI", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["BMI", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(BMI, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=BMI, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "BMI") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 62, y = 80, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 62, y = 78, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 62, y = 76, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/BMI-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = BMI, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "BMI") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 62, y = 50, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 62, y = 47, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 62, y = 44, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/BMI-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Cigarette Pack Years
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Cig.Pack.Years",] >= 0 & !is.na(df.cov.K["Cig.Pack.Years",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Cig.Pack.Years",] >= 0 & !is.na(df.cov.K["Cig.Pack.Years",])) ]))
Cig.Pack.Years      <- as.numeric(as.vector(df.cov.K["Cig.Pack.Years", (df.cov.K["Cig.Pack.Years",] >= 0 & !is.na(df.cov.K["Cig.Pack.Years",])) ]))

DNAm.Age.Rho  <- cor.test(Cig.Pack.Years, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Cig.Pack.Years, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Cig.Pack.Years, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Cig.Pack.Years, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Cig.Pack.Years", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Cig.Pack.Years", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Cig.Pack.Years", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Cig.Pack.Years", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(Cig.Pack.Years, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Cig.Pack.Years, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Cigarette Pack Years") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 27, y = 35, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 27, y = 33, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 27, y = 31, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/CigPackYears-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Cig.Pack.Years, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Cigarette Pack Years") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 27, y = 50, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 27, y = 47, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 27, y = 44, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/CigPackYears-DNAmAgeResidual.png", sep = '') )
p
dev.off()

##########################################################################
# Univariate Analysis - Current Smoker
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age           <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["SmokingYes",]) ]))
DNAm.Age.Residual  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["SmokingYes",]) ]))
Smoking            <- as.numeric(as.vector(df.cov.K["SmokingYes", !is.na(df.cov.K["SmokingYes",]) ]))

DNAm.Age.Rho  <- cor.test(Smoking, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Smoking, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Smoking, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Smoking, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Smoking", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Smoking", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Smoking", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Smoking", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
DNAm.Age.SmokingYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 1)]))
DNAm.Age.SmokingNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 0)]))
DNAm.Age.Residual.SmokingYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 1)]))
DNAm.Age.Residual.SmokingNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 0)]))


smoking               <- c("Current Smoker - No", "Current Smoker - Yes")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.SmokingNo), mean(DNAm.Age.SmokingYes)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.SmokingNo), sd(DNAm.Age.SmokingYes)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.SmokingNo), mean(DNAm.Age.Residual.SmokingYes)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.SmokingNo), sd(DNAm.Age.Residual.SmokingYes)) )

df.tmp <- data.frame(smoking, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=smoking, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "") + 
  labs(y = "Mean DNAm Age") + 
  annotate("text", x = 0.75, y = 71, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 69, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.SmokingNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.SmokingYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Smoking-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = smoking, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "") + 
  labs(y = "Mean DNAm Age Acceleration") + 
  annotate("text", x = 0.75, y = -4.4, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = -5.5, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.SmokingNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.Residual.SmokingYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Smoking-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Current Drinker
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["DrinkingYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["DrinkingYes",]) ]))
Drinking            <- as.numeric(as.vector(df.cov.K["DrinkingYes", !is.na(df.cov.K["DrinkingYes",]) ]))

DNAm.Age.Rho  <- cor.test(Drinking, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Drinking, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Drinking, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Drinking, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Drinking", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Drinking", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Drinking", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Drinking", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
DNAm.Age.DrinkingYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 1)]))
DNAm.Age.DrinkingNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 0)]))
DNAm.Age.Residual.DrinkingYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 1)]))
DNAm.Age.Residual.DrinkingNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 0)]))


drinking               <- c("Alcohol Consumption - No", "Alcohol Consumption - Yes")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.DrinkingNo), mean(DNAm.Age.DrinkingYes)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.DrinkingNo), sd(DNAm.Age.DrinkingYes)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.DrinkingNo), mean(DNAm.Age.Residual.DrinkingYes)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.DrinkingNo), sd(DNAm.Age.Residual.DrinkingYes)) )

df.tmp <- data.frame(drinking, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=drinking, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "") + 
  labs(y = "Mean DNAm Age") + 
  annotate("text", x = 0.75, y = 71, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 69, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.DrinkingNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.DrinkingYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Alcohol-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = drinking, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "") + 
  labs(y = "Mean DNAm Age Acceleration") + 
  annotate("text", x = 0.75, y = 22, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 21, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.DrinkingNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.Residual.DrinkingYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Alcohol-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Age at Menarche
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))
Menarche            <- as.numeric(as.vector(df.cov.K["Menarche", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))

DNAm.Age.Rho  <- cor.test(Menarche, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Menarche, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Menarche, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Menarche, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Menarche", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Menarche", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Menarche", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Menarche", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(Menarche, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Menarche, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at Menarche") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 8.5, y = 30, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 8.5, y = 28, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 8.5, y = 26, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/MenarcheAge-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Menarche, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at Menarche") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 16, y = 49, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 16, y = 46, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 16, y = 43, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/MenarcheAge-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Have you been pregnant
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["Been.PregnantYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Been.PregnantYes",]) ]))
Been.PregnantYes    <- as.numeric(as.vector(df.cov.K["Been.PregnantYes", !is.na(df.cov.K["Been.PregnantYes",]) ]))

DNAm.Age.Rho  <- cor.test(Been.PregnantYes, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Been.PregnantYes, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Been.PregnantYes, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Been.PregnantYes, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Been.PregnantYes", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Been.PregnantYes", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Been.PregnantYes", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Been.PregnantYes", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
DNAm.Age.Been.PregnantYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Been.PregnantYes",] == 1 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Been.PregnantNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Been.PregnantYes",] == 0 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Residual.Been.PregnantYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Been.PregnantYes",] == 1 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Residual.Been.PregnantNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Been.PregnantYes",] == 0 & !is.na(df.cov.K["Been.PregnantYes",]))]))

pregnant               <- c("Have you been pregnant? - No", "Have you been pregnant? - Yes")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.Been.PregnantNo), mean(DNAm.Age.Been.PregnantYes)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.Been.PregnantNo), sd(DNAm.Age.Been.PregnantYes)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.Been.PregnantNo), mean(DNAm.Age.Residual.Been.PregnantYes)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.Been.PregnantNo), sd(DNAm.Age.Residual.Been.PregnantYes)) )

df.tmp <- data.frame(pregnant, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=pregnant, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "") + 
  labs(y = "Mean DNAm Age") + 
  annotate("text", x = 0.75, y = 64.5, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 61.5, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.Been.PregnantNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.Been.PregnantYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/HaveYouBeenPregnant-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = pregnant, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "") + 
  labs(y = "Mean DNAm Age Acceleration") + 
  annotate("text", x = 2.25, y = 25, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 2.25, y = 24, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.Been.PregnantNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.Residual.Been.PregnantYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/HaveYouBeenPregnant-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Times Pregnant
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))
Times.Pregnant      <- as.numeric(as.vector(df.cov.K["Times.Pregnant", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))

DNAm.Age.Rho  <- cor.test(Times.Pregnant, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Times.Pregnant, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Times.Pregnant, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Times.Pregnant, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Times.Pregnant", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Times.Pregnant", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Times.Pregnant", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Times.Pregnant", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(Times.Pregnant, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Times.Pregnant, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Times Pregnant") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 6, y = 36, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 6, y = 34, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 6, y = 32, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/TimesPregnant-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Times.Pregnant, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Times Pregnant") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 6, y = 51, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 6, y = 48, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 6, y = 45, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/TimesPregnant-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Parity
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Parity",] >= 0 & df.cov.K["Parity",] != 41 & !is.na(df.cov.K["Parity",]) ) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Parity",] >= 0 & df.cov.K["Parity",] != 41 & !is.na(df.cov.K["Parity",])) ]))
Parity              <- as.numeric(as.vector(df.cov.K["Parity", (df.cov.K["Parity",] >= 0 & df.cov.K["Parity",] != 41 & !is.na(df.cov.K["Parity",])) ]))

DNAm.Age.Rho  <- cor.test(Parity, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Parity, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Parity, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Parity, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Parity", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Parity", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Parity", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Parity", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(Parity, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Parity, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Parity") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 6, y = 35, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 6, y = 33, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 6, y = 31, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Parity-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Parity, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Parity") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 6, y = 50, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 6, y = 47, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 6, y = 44, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/Parity-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Age at first birth
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))
Age.FB              <- as.numeric(as.vector(df.cov.K["Age.FB", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))

DNAm.Age.Rho  <- cor.test(Age.FB, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Age.FB, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Age.FB, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Age.FB, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Age.FB", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Age.FB", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Age.FB", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Age.FB", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(Age.FB, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Age.FB, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at First Birth") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 40, y = 80, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 40, y = 78, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 40, y = 76, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/AgeAtFirstBirth-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Age.FB, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at First Birth") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 40, y = 50, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 40, y = 47, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 40, y = 44, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/AgeAtFirstBirth-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Menstural Status
##########################################################################

##### Plot Correlations #####
MenopausePost                   <- as.numeric(as.vector(df.cov.K["MenopausePost-menopausal", (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
MenopausePre                    <- as.numeric(as.vector(df.cov.K["MenopausePre-menopausal", (df.cov.K["MenopausePre-menopausal",] == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
DNAm.Age.MenopausePost          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
DNAm.Age.MenopausePre           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["MenopausePre-menopausal",] == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
DNAm.Age.Residual.MenopausePost <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
DNAm.Age.Residual.MenopausePre  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["MenopausePre-menopausal",] == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))

menopause              <- c("0.Pre-menopausal", "1.Post-menopausal")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.MenopausePre), mean(DNAm.Age.MenopausePost)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.MenopausePre), sd(DNAm.Age.MenopausePost)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.MenopausePre), mean(DNAm.Age.Residual.MenopausePost)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.MenopausePre), sd(DNAm.Age.Residual.MenopausePost)) )

##### Correlate DNAm Age and Residual to Menopause #####
MenopausePre[MenopausePre == 1] <- 0
DNAm.Age            <- as.numeric( c(DNAm.Age.MenopausePre, DNAm.Age.MenopausePost) )
DNAm.Age.Residual   <- as.numeric( c(DNAm.Age.Residual.MenopausePre, DNAm.Age.Residual.MenopausePost) )
Menopause           <- as.numeric( c(MenopausePre, MenopausePost) )

DNAm.Age.Rho  <- cor.test(Menopause, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Menopause, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Menopause, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Menopause, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Menopause", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Menopause", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Menopause", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Menopause", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

df.tmp <- data.frame(menopause, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=menopause, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Menopause Status") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.75, y = 69, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 66, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.MenopausePre), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.MenopausePost), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/MenopauseStatus-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = menopause, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Menopause Status") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.75, y = -5, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = -6, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.MenopausePre), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.Residual.MenopausePost), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/MenopauseStatus-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Age at Menopause
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Menopause.Age",] >= 0 & !is.na(df.cov.K["Menopause.Age",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Menopause.Age",] >= 0 & !is.na(df.cov.K["Menopause.Age",])) ]))
Age.Menopause       <- as.numeric(as.vector(df.cov.K["Menopause.Age", (df.cov.K["Menopause.Age",] >= 0 & !is.na(df.cov.K["Menopause.Age",])) ]))

DNAm.Age.Rho  <- cor.test(Age.Menopause, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Age.Menopause, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Age.Menopause, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Age.Menopause, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Age.Menopause", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Age.Menopause", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Age.Menopause", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Age.Menopause", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
df.tmp <- data.frame(Age.Menopause, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Age.Menopause, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at Menopause") +
  labs(y = "DNAm Age") + 
  annotate("text", x = 28, y = 35, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 28, y = 33, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 28, y = 31, label = paste("n = ", length(DNAm.Age), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/MenopauseAge-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Age.Menopause, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at Menopause") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 28, y = -22, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 28, y = -24, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 28, y = -26, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/MenopauseAge-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate Analysis - Vitamin Use
##########################################################################

##### Correlate DNAm Age and Residual to Race #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["VDYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["VDYes",]) ]))
VD                  <- as.numeric(as.vector(df.cov.K["VDYes", complete.cases(df.cov.K["VDYes",]) ]))

DNAm.Age.Rho  <- cor.test(VD, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(VD, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(VD, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(VD, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["VD", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["VD", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["VD", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["VD", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

##### Plot Correlations #####
DNAm.Age.VDYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["VDYes",] == 1]))
DNAm.Age.VDNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["VDYes",] == 0]))
DNAm.Age.Residual.VDYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["VDYes",] == 1]))
DNAm.Age.Residual.VDNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["VDYes",] == 0]))

vd                     <- c("Multivitamin Use - No", "Multivitamin Use - Yes")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.VDNo), mean(DNAm.Age.VDYes)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.VDNo), sd(DNAm.Age.VDYes)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.VDNo), mean(DNAm.Age.Residual.VDYes)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.VDNo), sd(DNAm.Age.Residual.VDYes)) )

df.tmp <- data.frame(vd, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=vd, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.75, y = 70, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = 67, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age[1]+2, label = paste("n = ", length(DNAm.Age.VDNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age[2]+2, label = paste("n = ", length(DNAm.Age.VDYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/VDUse-DNAmAge.png", sep = '') )
p
dev.off()

p <- ggplot(df.tmp, aes(x = vd, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.75, y = -5, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.75, y = -6.2, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.VDNo), sep = ""), col = "blue") +
  annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]+0.5, label = paste("n = ", length(DNAm.Age.Residual.VDYes), sep = ""), col = "blue") +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

png( paste(model.dir, "FactorAssociation/VDUse-DNAmAgeResidual.png", sep = '') )
p
dev.off()


##########################################################################
# Univariate/Multivariate Analysis Summary
##########################################################################

##### Uni #####
df.univar[which(df.univar$DNAm.Age.pval < 0.05),1:2]
df.univar[which(df.univar$DNAm.Age.Residual.pval < 0.05),3:4]

df.univar

##### Multi #####
summary(lm.DNAm.Age)
summary(lm.DNAm.Age.Residual)

