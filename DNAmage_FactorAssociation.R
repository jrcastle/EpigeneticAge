setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(moments)

rm(list=ls()); gc();

seed        <- "123"
model.dir   <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')
meth.file.K <- paste("data/meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs_seed", seed, ".txt", sep = "")
cov.file.K  <- paste("data/cov_K_vali_seed", seed, ".txt", sep = "")


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
# Multivariate Correlation Analysis
##########################################################################

DNAm.Age          <- as.numeric(as.vector(df.cov.K["DNAm.Age",]))
DNAm.Age.Residual <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",]))
Race              <- as.numeric(as.vector(df.cov.K["RaceWhite",]))
BMI               <- as.numeric(as.vector(df.cov.K["LocationUrban",]))
Cig.Pack.Years    <- as.numeric(as.vector(df.cov.K["Cig.Pack.Years",]))
Drinking          <- as.numeric(as.vector(df.cov.K["DrinkingYes",]))
Menarche          <- as.numeric(as.vector(df.cov.K["Menarche",]))
Been.PregnantYes  <- as.numeric(as.vector(df.cov.K["Been.PregnantYes",]))
Times.Pregnant    <- as.numeric(as.vector(df.cov.K["Times.Pregnant",]))
Parity            <- as.numeric(as.vector(df.cov.K["Parity",]))
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
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["RaceWhite",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["RaceWhite",]) ]))
Race                <- as.numeric(as.vector(df.cov.K["RaceWhite", complete.cases(df.cov.K["RaceWhite",]) ]))

DNAm.Age.Rho  <- cor.test(Race, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Race, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Race, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Race, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Race", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Race", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Race", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Race", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

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

p <- ggplot(df.tmp, aes(x=race, y=Mean.DNAm.Age)) + 
     geom_bar(stat="identity", width = 0.5) + 
     geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
     labs(x = "Race") + 
     labs(y = "Mean DNAm Age") + 
     annotate("text", x = 0.75, y = 64.5, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 61.5, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 0.75, y = 9, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 8.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 0.75, y = 66, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 63, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 0.75, y = -4.5, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = -5.2, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
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
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["BMI",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["BMI",]) ]))
BMI                 <- as.numeric(as.vector(df.cov.K["BMI", complete.cases(df.cov.K["BMI",]) ]))

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
     annotate("text", x = 17, y = 72, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 17, y = 70, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 17, y = 68, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 17, y = 15, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 17, y = 13, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 17, y = 11, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
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
     annotate("text", x = 27, y = 70, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 27, y = 68, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 27, y = 66, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 27, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 27, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 27, y = 16, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
# Univariate Analysis - Current Drinker
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["DrinkingYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["DrinkingYes",]) ]))
Drinking            <- as.numeric(as.vector(df.cov.K["DrinkingYes", complete.cases(df.cov.K["DrinkingYes",]) ]))

DNAm.Age.Rho  <- cor.test(Drinking, DNAm.Age, method = "spearman" )$estimate
DNAm.Age.pval <- cor.test(Drinking, DNAm.Age, method = "spearman" )$p.value

DNAm.Age.Residual.Rho  <- cor.test(Drinking, DNAm.Age.Residual, method = "spearman" )$estimate
DNAm.Age.Residual.pval <- cor.test(Drinking, DNAm.Age.Residual, method = "spearman" )$p.value

df.univar["Drinking", "DNAm.Age.Rho"]           <- DNAm.Age.Rho
df.univar["Drinking", "DNAm.Age.pval"]          <- DNAm.Age.pval
df.univar["Drinking", "DNAm.Age.Residual.Rho"]  <- DNAm.Age.Residual.Rho
df.univar["Drinking", "DNAm.Age.Residual.pval"] <- DNAm.Age.Residual.pval

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

p <- ggplot(df.tmp, aes(x=drinking, y=Mean.DNAm.Age)) + 
     geom_bar(stat="identity", width = 0.5) + 
     geom_errorbar(aes(ymin = Mean.DNAm.Age - std.DNAm.Age, ymax = Mean.DNAm.Age + std.DNAm.Age), width = 0.2) + 
     labs(x = "") + 
     labs(y = "Mean DNAm Age") + 
     annotate("text", x = 0.75, y = 64.5, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 61.5, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 0.75, y = -5.4, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = -6.1, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 16, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 16, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 16, y = 71, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 16, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 16, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 16, y = 16, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
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
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["Been.PregnantYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["Been.PregnantYes",]) ]))
Been.PregnantYes    <- as.numeric(as.vector(df.cov.K["Been.PregnantYes", complete.cases(df.cov.K["Been.PregnantYes",]) ]))

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
     annotate("text", x = 0.75, y = 8, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 7.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 6, y = 76, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 6, y = 74, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 6, y = 72, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 6, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 6, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 6, y = 16, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
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
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))
Parity              <- as.numeric(as.vector(df.cov.K["Parity", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))

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
     annotate("text", x = 3.5, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 3.5, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 3.5, y = 71, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 3.5, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 3.5, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 3.5, y = 16, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
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
     annotate("text", x = 40, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 40, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 40, y = 71, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 40, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 40, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 40, y = 16, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
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
     annotate("text", x = 0.75, y = 59, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 56, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 0.75, y = 10, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 9.3, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.MenopausePre), sep = ""), col = "blue") +
     annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]-0.5, label = paste("n = ", length(DNAm.Age.Residual.MenopausePost), sep = ""), col = "blue") +
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
     annotate("text", x = 33, y = 75, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 33, y = 73, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
     annotate("text", x = 33, y = 71, label = paste("n = ", length(DNAm.Age), sep = "")) +
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
     annotate("text", x = 33, y = 20, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 33, y = 18, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 33, y = 16, label = paste("n = ", length(DNAm.Age.Residual), sep = "")) +
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
# Univariate Analysis - VD Use
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

vd                     <- c("Vitamin D Use - No", "Vitamin D Use - Yes")
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
     annotate("text", x = 0.75, y = 63, label = paste("Spearman Rho: ", round(DNAm.Age.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = 60, label = paste("p value: ", round(DNAm.Age.pval, digits = 4), sep = "")) +
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
     annotate("text", x = 0.75, y = -5.5, label = paste("Spearman Rho: ", round(DNAm.Age.Residual.Rho, digits = 3), sep = "")) +
     annotate("text", x = 0.75, y = -6.2, label = paste("p value: ", round(DNAm.Age.Residual.pval, digits = 4), sep = "")) +
     annotate("text", x = 1, y = Mean.DNAm.Age.Residual[1]+0.5, label = paste("n = ", length(DNAm.Age.Residual.VDNo), sep = ""), col = "blue") +
     annotate("text", x = 2, y = Mean.DNAm.Age.Residual[2]-0.5, label = paste("n = ", length(DNAm.Age.Residual.VDYes), sep = ""), col = "blue") +
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
