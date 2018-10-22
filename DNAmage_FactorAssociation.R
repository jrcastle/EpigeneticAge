setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)

rm(list=ls()); gc();

model.dir   <- "cpgs_in_KNT_imputed/"
cov.file.K  <- "data/cov_K_vali.txt"
meth.file.K <- "data/meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs.txt"


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

##### METH K ####
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

##### SAMPLE AGES K #####
df.cov.K      <- read.table(cov.file.K, header = TRUE, row.names = 1, sep = '\t')
sample.ages.K <- as.numeric(as.vector(df.cov.K["Age",]))


###########################################################################################
# MODEL COEFFICIENTS 
###########################################################################################
clock.cpg.coef <- read.csv(paste(model.dir, "model_coefficients.csv", sep = ''), stringsAsFactors = FALSE)
clock.cpg.coef <- clock.cpg.coef[ c("model.coefficients.name", "model.coefficients.x") ]
clock.cpg.coef[clock.cpg.coef$model.coefficients.name == "(Intercept)", "model.coefficients.name"] <- "Intercept"

##### SORT K #####
meth.data.K <- meth.data.K[tolower(order(meth.data.K$position)),]
clock.cpg.coef <- clock.cpg.coef[order(clock.cpg.coef$model.coefficients.name),]
meth.data.K$position == clock.cpg.coef$model.coefficients.name


###########################################################################################
# PREDICT
###########################################################################################

##### K #####
meth.data.K$position <- NULL
X <- data.matrix(meth.data.K)
beta <- data.matrix(clock.cpg.coef$model.coefficients.x)

result.K <- t(X) %*% beta
result.K <- sapply(result.K, anti.trafo)

residual.K <- result.K - sample.ages.K 
mean.error.K <- mean(residual.K)
stdev.error.K <- sd(residual.K)

##### Add result and residual to cov dataframe #####
df.cov.K["DNAm.Age",] <- result.K
df.cov.K["DNAm.Age.Residual",] <- residual.K


##########################################################################
# Race
##########################################################################

##### Correlate DNAm Age and Residual to Race #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["RaceWhite",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["RaceWhite",]) ]))
Race                <- as.numeric(as.vector(df.cov.K["RaceWhite", complete.cases(df.cov.K["RaceWhite",]) ]))

DNAm.Age.Rho  <- cor.test(Race, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Race, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Race, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Race, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
DNAm.Age.White            <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["RaceWhite",] == 1]))
DNAm.Age.AfrAmer          <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["RaceWhite",] == 0]))
DNAm.Age.Residual.White   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["RaceWhite",] == 1]))
DNAm.Age.Residual.AfrAmer <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["RaceWhite",] == 0]))

race                   <- c("African American", "White")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.AfrAmer), mean(DNAm.Age.White)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.AfrAmer), sd(DNAm.Age.White)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.AfrAmer), mean(DNAm.Age.Residual.White)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.AfrAmer), sd(DNAm.Age.Residual.White)) )

df.tmp <- data.frame(race, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

ggplot(df.tmp, aes(x=race, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Race") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.65, y = 56, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 53, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = race, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Race") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.65, y = 6, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 5.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Location
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

df.tmp <- data.frame(location, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

ggplot(df.tmp, aes(x=location, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Location") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.65, y = 56, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 53, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = location, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Location") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.65, y = 9, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 8.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# BMI
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["BMI",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["BMI",]) ]))
BMI                <- as.numeric(as.vector(df.cov.K["BMI", complete.cases(df.cov.K["BMI",]) ]))

DNAm.Age.Rho  <- cor.test(BMI, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(BMI, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(BMI, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(BMI, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
df.tmp <- data.frame(BMI, DNAm.Age, DNAm.Age.Residual)

ggplot(df.tmp, aes(x=BMI, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "BMI") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 15, y = 72, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 15, y = 70, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = BMI, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "BMI") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 16, y = 15, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 16, y = 13, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Cigarette Pack Years
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Cig.Pack.Years",] >= 0 & !is.na(df.cov.K["Cig.Pack.Years",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Cig.Pack.Years",] >= 0 & !is.na(df.cov.K["Cig.Pack.Years",])) ]))
Cig.Pack.Years      <- as.numeric(as.vector(df.cov.K["Cig.Pack.Years", (df.cov.K["Cig.Pack.Years",] >= 0 & !is.na(df.cov.K["Cig.Pack.Years",])) ]))

DNAm.Age.Rho  <- cor.test(Cig.Pack.Years, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Cig.Pack.Years, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Cig.Pack.Years, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Cig.Pack.Years, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
df.tmp <- data.frame(Cig.Pack.Years, DNAm.Age, DNAm.Age.Residual)

ggplot(df.tmp, aes(x=Cig.Pack.Years, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Cigarette Pack Years") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 29, y = 70, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 29, y = 68, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = Cig.Pack.Years, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Cigarette Pack Years") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 27, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 27, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Current Drinker
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["DrinkingYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["DrinkingYes",]) ]))
Drinking            <- as.numeric(as.vector(df.cov.K["DrinkingYes", complete.cases(df.cov.K["DrinkingYes",]) ]))

DNAm.Age.Rho  <- cor.test(Drinking, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Drinking, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Drinking, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Drinking, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
DNAm.Age.DrinkingYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["DrinkingYes",] == 1]))
DNAm.Age.DrinkingNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["DrinkingYes",] == 0]))
DNAm.Age.Residual.DrinkingYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["DrinkingYes",] == 1]))
DNAm.Age.Residual.DrinkingNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["DrinkingYes",] == 0]))

drinking               <- c("Alcohol Consumption - No", "Alcohol Consumption - Yes")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.DrinkingNo), mean(DNAm.Age.DrinkingYes)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.DrinkingNo), sd(DNAm.Age.DrinkingYes)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.DrinkingNo), mean(DNAm.Age.Residual.DrinkingYes)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.DrinkingNo), sd(DNAm.Age.Residual.DrinkingYes)) )

df.tmp <- data.frame(drinking, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

ggplot(df.tmp, aes(x=drinking, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Race") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.65, y = 56, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 53, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = drinking, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Race") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.65, y = 6, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 5.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Age at Menarche
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))
Menarche            <- as.numeric(as.vector(df.cov.K["Menarche", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))

DNAm.Age.Rho  <- cor.test(Menarche, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Menarche, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Menarche, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Menarche, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
df.tmp <- data.frame(Menarche, DNAm.Age, DNAm.Age.Residual)

ggplot(df.tmp, aes(x=Menarche, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at Menarche") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 16.5, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 16.5, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = Menarche, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at Menarche") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 16.5, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 16.5, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Have you been pregnant
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["Been.PregnantYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["Been.PregnantYes",]) ]))
Been.PregnantYes    <- as.numeric(as.vector(df.cov.K["DrinkingYes", complete.cases(df.cov.K["Been.PregnantYes",]) ]))

DNAm.Age.Rho  <- cor.test(Been.PregnantYes, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Been.PregnantYes, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Been.PregnantYes, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Been.PregnantYes, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
DNAm.Age.Been.PregnantYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Been.PregnantYes",] == 1 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Been.PregnantNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Been.PregnantYes",] == 0 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Residual.Been.PregnantYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Been.PregnantYes",] == 1 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Residual.Been.PregnantNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Been.PregnantYes",] == 0 & !is.na(df.cov.K["Been.PregnantYes",]))]))

pregnant               <- c("Have you been pregnant? - No", "Have you been pregnant? - Yes")
Mean.DNAm.Age          <- as.numeric( c(mean(DNAm.Age.DrinkingNo), mean(DNAm.Age.DrinkingYes)) ) 
std.DNAm.Age           <- as.numeric( c(sd(DNAm.Age.DrinkingNo), sd(DNAm.Age.DrinkingYes)) ) 
Mean.DNAm.Age.Residual <- as.numeric( c(mean(DNAm.Age.Residual.DrinkingNo), mean(DNAm.Age.Residual.DrinkingYes)) ) 
std.DNAm.Age.Residual  <- as.numeric( c(sd(DNAm.Age.Residual.DrinkingNo), sd(DNAm.Age.Residual.DrinkingYes)) )

df.tmp <- data.frame(pregnant, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

ggplot(df.tmp, aes(x=pregnant, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.65, y = 56, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 53, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = pregnant, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.65, y = 6, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 5.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Times Pregnant
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))
Times.Pregnant      <- as.numeric(as.vector(df.cov.K["Times.Pregnant", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))

DNAm.Age.Rho  <- cor.test(Times.Pregnant, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Times.Pregnant, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Times.Pregnant, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Times.Pregnant, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
df.tmp <- data.frame(Times.Pregnant, DNAm.Age, DNAm.Age.Residual)

ggplot(df.tmp, aes(x=Times.Pregnant, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Times Pregnant") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 6.5, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 6.5, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = Times.Pregnant, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Times Pregnant") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 6.5, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 6.5, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Parity
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))
Parity              <- as.numeric(as.vector(df.cov.K["Parity", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))

DNAm.Age.Rho  <- cor.test(Parity, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Parity, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Parity, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Parity, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
df.tmp <- data.frame(Parity, DNAm.Age, DNAm.Age.Residual)

ggplot(df.tmp, aes(x=Parity, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Parity") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 3.5, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 3.5, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = Parity, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Parity") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 3.5, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 3.5, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Age at first birth
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))
Age.FB              <- as.numeric(as.vector(df.cov.K["Age.FB", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))

DNAm.Age.Rho  <- cor.test(Age.FB, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Age.FB, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Age.FB, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Age.FB, DNAm.Age.Residual, method = "spearman" )$p.value

##### Plot Correlations #####
df.tmp <- data.frame(Age.FB, DNAm.Age, DNAm.Age.Residual)

ggplot(df.tmp, aes(x=Age.FB, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at First Birth") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 40, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 40, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = Age.FB, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at First Birth") +
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 40, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 40, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )


##########################################################################
# Menstural Status
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

df.tmp <- data.frame(menopause, Mean.DNAm.Age, std.DNAm.Age, Mean.DNAm.Age.Residual, std.DNAm.Age.Residual)

ggplot(df.tmp, aes(x=menopause, y=Mean.DNAm.Age)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
  labs(x = "Menopause Status") + 
  labs(y = "DNAm Age") + 
  annotate("text", x = 0.65, y = 56, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 53, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df.tmp, aes(x = menopause, y = Mean.DNAm.Age.Residual)) + 
  geom_bar(stat="identity", width = 0.5) + 
  geom_errorbar(aes(ymin = Mean.DNAm.Age.Residual - std.DNAm.Age.Residual, ymax = Mean.DNAm.Age.Residual + std.DNAm.Age.Residual), width = 0.2) + 
  labs(x = "Menopause Status") + 
  labs(y = "DNAm Age Acceleration") + 
  annotate("text", x = 0.65, y = 9, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
  annotate("text", x = 0.65, y = 8.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
  theme_bw() + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5)
  )



