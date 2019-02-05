rm(list=ls()); gc();
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")
library(ggplot2)
library(moments)
library(stargazer)
source("debug_contr_error.R")
source("plot_functions.R")

seed        <- "123"
model.dir   <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')
meth.file.K <- paste("data/meth_K_cpgs_in_KNT_imputed_vali_ClockCpGs_seed", seed, ".txt", sep = "")
cov.file.K  <- paste("data/cov_K_vali_seed", seed, ".txt", sep = "")

#seed        <- ""
#model.dir   <- "WhiteBlack_DiffMethAnalysis/ElasticNet/"
#meth.file.K <- "data/meth_K_WhiteBlackDiffMethCpGs_imputed_vali_ClockCpGs.txt"
#cov.file.K  <- "data/cov_K_vali_seed123.txt"

age.covariate  = TRUE
model.residual = TRUE

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
if(model.residual){
  residual.K <- as.numeric(as.vector(lm(result.K ~ sample.ages.K)$residual))
}
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
BMI               <- as.numeric(as.vector(df.cov.K["BMI",]))
Cig.Pack.Years    <- as.numeric(as.vector(df.cov.K["Cig.Pack.Years",]))
Drinking          <- as.numeric(as.vector(df.cov.K["DrinkingYes",]))
Menarche          <- as.numeric(as.vector(df.cov.K["Menarche",]))
Been.PregnantYes  <- as.numeric(as.vector(df.cov.K["Been.PregnantYes",]))
Times.Pregnant    <- as.numeric(as.vector(df.cov.K["Times.Pregnant",]))
Parity            <- as.numeric(as.vector(df.cov.K["Parity",]))
Age.FB            <- as.numeric(as.vector(df.cov.K["Age.FB",]))
Menopause.Age     <- as.numeric(as.vector(df.cov.K["Menopause.Age",]))
VDYes             <- as.numeric(as.vector(df.cov.K["VDYes",]))
SmokingYes        <- as.numeric(as.vector(df.cov.K["SmokingYes",]))

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

# Factorize categorical variables
Race[which(Race == 0)] <- "African American"
Race[which(Race == 1)] <- "Caucasian"
Race <- as.factor(Race)
#if(anyNA(Race)) { Race <- addNA(Race) }

Drinking[which(Drinking == 0)] <- "Current Alcohol Use - No"
Drinking[which(Drinking == 1)] <- "Current Alcohol Use - Yes"
Drinking <- as.factor(Drinking)
#if(anyNA(Drinking)) { Drinking <- addNA(Drinking) }

SmokingYes[which(SmokingYes == 0)] <- "Current Tobacco Use - No"
SmokingYes[which(SmokingYes == 1)] <- "Current Tobacco Use - Yes"
SmokingYes <- as.factor(SmokingYes)
#if(anyNA(Smoking)) { Smoking <- addNA(Smoking) }

Been.PregnantYes[which(Been.PregnantYes == 0)] <- "Been Pregnant - No"
Been.PregnantYes[which(Been.PregnantYes == 1)] <- "Been Pregnant - Yes"
Been.PregnantYes <- as.factor(Been.PregnantYes)
#if(anyNA(Been.PregnantYes)) { Been.PregnantYes <- addNA(Been.PregnantYes) }

VDYes[which(VDYes == 0)] <- "Multivitamin Use - No"
VDYes[which(VDYes == 1)] <- "Multivitamin Use - Yes"
VDYes <- as.factor(VDYes)
#if(anyNA(VDYes)) { VDYes <- addNA(VDYes) }

Menopause[which(Menopause == 0)] <- "Pre-menopausal"
Menopause[which(Menopause == 1)] <- "Post-menopausal"
Menopause <- as.factor(Menopause)
#if(anyNA(Menopause)) { Menopause <- addNA(Menopause) }

Location[which(Location == 0)] <- "Rural"
Location[which(Location == 1)] <- "Urban"
Location <- as.factor(Location)
#if(anyNA(Location)) { Location <- addNA(Location) }

#df <- data.frame(DNAm.Age, Race, BMI, Cig.Pack.Years, Drinking, SmokingYes,
#                 Menarche, Been.PregnantYes, Times.Pregnant, Parity, Age.FB, 
#                 Menopause.Age, VDYes, Menopause, Location, sample.ages.K)
#colnames(df) <- c("DNAm.Age", "Race", "BMI", "Cig.Pack.Years", "Drinking", "SmokingYes", 
#                  "Menarche", "Been.PregnantYes", "Times.Pregnant", "Parity", "Age.FB", 
#                  "Menopause.Age", "VDYes", "Menopause", "Location", "sample.ages.K")
#df<-NA_preproc(df)

if( age.covariate ){
  lm.DNAm.Age <- glm(DNAm.Age ~ 
                     Race +
                     BMI +
                     #Cig.Pack.Years +
                     Drinking + 
                     SmokingYes + 
                     Menarche +
                     #Been.PregnantYes +
                     Times.Pregnant + 
                     Parity + 
                     Age.FB +
                     #Menopause.Age +
                     VDYes +
                     Menopause +
                     Location + 
                     sample.ages.K
                    )
}else{
  lm.DNAm.Age <- glm(DNAm.Age ~ 
                     Race + 
                     BMI + 
                     #Cig.Pack.Years + 
                     Drinking + 
                     SmokingYes + 
                     Menarche + 
                     #Been.PregnantYes + 
                     Times.Pregnant + 
                     Parity + 
                     Age.FB + 
                     #Menopause.Age + 
                     VDYes + 
                     Menopause + 
                     Location
                    )
}

if( age.covariate ){
  lm.DNAm.Age.Residual <- lm(DNAm.Age.Residual ~ 
                             Race + 
                             BMI + 
                             #Cig.Pack.Years + 
                             Drinking + 
                             SmokingYes + 
                             Menarche + 
                             #Been.PregnantYes + 
                             Times.Pregnant + 
                             Parity + 
                             Age.FB + 
                             #Menopause.Age + 
                             VDYes + 
                             Menopause + 
                             Location 
                             sample.ages.K
                            )
} else{
  lm.DNAm.Age.Residual <- lm(DNAm.Age.Residual ~ 
                             Race + 
                             BMI + 
                             #Cig.Pack.Years + 
                             Drinking + 
                             SmokingYes + 
                             Menarche + 
                             #Been.PregnantYes + 
                             Times.Pregnant + 
                             Parity + 
                             Age.FB + 
                             #Menopause.Age + 
                             VDYes + 
                             Menopause + 
                             Location
                            )
}


##########################################################################
# DNAm Age Univariate Analysis
##########################################################################

##### Race DNAm Age #####
DNAm.Age.tmp     <- as.numeric(as.vector(df.cov.K["DNAm.Age",  !is.na(df.cov.K["RaceWhite",]) ]))
Race             <- as.numeric(as.vector(df.cov.K["RaceWhite", !is.na(df.cov.K["RaceWhite",]) ]))
sample.age.tmp   <- as.numeric(as.vector(df.cov.K["Age",       !is.na(df.cov.K["RaceWhite",]) ]))
if( age.covariate ){
  lm.Race.DNAm.Age <- glm(DNAm.Age.tmp ~ Race + sample.age.tmp)
}else{
  lm.Race.DNAm.Age <- glm(DNAm.Age.tmp ~ Race)
}

##### BMI DNAm Age #####
DNAm.Age.tmp    <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["BMI",]) ]))
BMI             <- as.numeric(as.vector(df.cov.K["BMI",      !is.na(df.cov.K["BMI",]) ]))
sample.age.tmp  <- as.numeric(as.vector(df.cov.K["Age",      !is.na(df.cov.K["BMI",]) ]))
if( age.covariate ){
  lm.BMI.DNAm.Age <- glm(DNAm.Age.tmp ~ BMI + sample.age.tmp)
}else{
  lm.BMI.DNAm.Age <- glm(DNAm.Age.tmp ~ BMI)
}

##### Smoking DNAm Age #####
DNAm.Age.tmp         <- as.numeric(as.vector(df.cov.K["DNAm.Age",    !is.na(df.cov.K["SmokingYes",]) ]))
SmokingYes           <- as.numeric(as.vector(df.cov.K["SmokingYes", !is.na(df.cov.K["SmokingYes",]) ]))
sample.age.tmp       <- as.numeric(as.vector(df.cov.K["Age",         !is.na(df.cov.K["SmokingYes",]) ]))
if( age.covariate ){
  lm.SmokingYes.DNAm.Age <- glm(DNAm.Age.tmp ~ SmokingYes + sample.age.tmp)
}else{
  lm.SmokingYes.DNAm.Age <- glm(DNAm.Age.tmp ~ SmokingYes)
}

##### Drinking DNAm Age #####
DNAm.Age.tmp         <- as.numeric(as.vector(df.cov.K["DNAm.Age",    !is.na(df.cov.K["DrinkingYes",]) ]))
Drinking             <- as.numeric(as.vector(df.cov.K["DrinkingYes", !is.na(df.cov.K["DrinkingYes",]) ]))
sample.age.tmp       <- as.numeric(as.vector(df.cov.K["Age",         !is.na(df.cov.K["DrinkingYes",]) ]))
if( age.covariate ){
  lm.Drinking.DNAm.Age <- glm(DNAm.Age.tmp ~ Drinking + sample.age.tmp)
}else{
  lm.Drinking.DNAm.Age <- glm(DNAm.Age.tmp ~ Drinking)
}

##### Menarche DNAm Age #####
DNAm.Age.tmp         <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["Menarche",]) ]))
Menarche             <- as.numeric(as.vector(df.cov.K["Menarche", !is.na(df.cov.K["Menarche",]) ]))
sample.age.tmp       <- as.numeric(as.vector(df.cov.K["Age",      !is.na(df.cov.K["Menarche",]) ]))
if( age.covariate ){
  lm.Menarche.DNAm.Age <- glm(DNAm.Age.tmp ~ Menarche + sample.age.tmp)
}else{
  lm.Menarche.DNAm.Age <- glm(DNAm.Age.tmp ~ Menarche)
}

##### Been.PregnantYes DNAm Age #####
DNAm.Age.tmp                 <- as.numeric(as.vector(df.cov.K["DNAm.Age",         !is.na(df.cov.K["Been.PregnantYes",]) ]))
Been.PregnantYes             <- as.numeric(as.vector(df.cov.K["Been.PregnantYes", !is.na(df.cov.K["Been.PregnantYes",]) ]))
sample.age.tmp               <- as.numeric(as.vector(df.cov.K["Age",              !is.na(df.cov.K["Been.PregnantYes",]) ]))
if( age.covariate ){
  lm.Been.PregnantYes.DNAm.Age <- glm(DNAm.Age.tmp ~ Been.PregnantYes + sample.age.tmp)
}else{
  lm.Been.PregnantYes.DNAm.Age <- glm(DNAm.Age.tmp ~ Been.PregnantYes)
}

##### Times.Pregnant DNAm Age #####
DNAm.Age.tmp               <- as.numeric(as.vector(df.cov.K["DNAm.Age",       !is.na(df.cov.K["Times.Pregnant",]) ]))
Times.Pregnant             <- as.numeric(as.vector(df.cov.K["Times.Pregnant", !is.na(df.cov.K["Times.Pregnant",]) ]))
sample.age.tmp             <- as.numeric(as.vector(df.cov.K["Age",            !is.na(df.cov.K["Times.Pregnant",]) ]))
if( age.covariate ){
  lm.Times.Pregnant.DNAm.Age <- glm(DNAm.Age.tmp ~ Times.Pregnant + sample.age.tmp)
}else{
  lm.Times.Pregnant.DNAm.Age <- glm(DNAm.Age.tmp ~ Times.Pregnant)
}

##### Parity DNAm Age #####
DNAm.Age.tmp       <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["Parity",]) ]))
Parity             <- as.numeric(as.vector(df.cov.K["Parity",   !is.na(df.cov.K["Parity",]) ]))
sample.age.tmp     <- as.numeric(as.vector(df.cov.K["Age",      !is.na(df.cov.K["Parity",]) ]))
if( age.covariate ){
  lm.Parity.DNAm.Age <- glm(DNAm.Age.tmp ~ Parity + sample.age.tmp)
}else{
  lm.Parity.DNAm.Age <- glm(DNAm.Age.tmp ~ Parity)
}

##### Age.FB DNAm Age #####
DNAm.Age.tmp       <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["Age.FB",]) ]))
Age.FB             <- as.numeric(as.vector(df.cov.K["Age.FB",   !is.na(df.cov.K["Age.FB",]) ]))
sample.age.tmp     <- as.numeric(as.vector(df.cov.K["Age",      !is.na(df.cov.K["Age.FB",]) ]))
if( age.covariate ){
  lm.Age.FB.DNAm.Age <- glm(DNAm.Age.tmp ~ Age.FB + sample.age.tmp)
}else{
  lm.Age.FB.DNAm.Age <- glm(DNAm.Age.tmp ~ Age.FB)
}

##### Menopause.Age DNAm Age #####
DNAm.Age.tmp              <- as.numeric(as.vector(df.cov.K["DNAm.Age",      !is.na(df.cov.K["Menopause.Age",]) ]))
Menopause.Age             <- as.numeric(as.vector(df.cov.K["Menopause.Age", !is.na(df.cov.K["Menopause.Age",]) ]))
sample.age.tmp            <- as.numeric(as.vector(df.cov.K["Age",           !is.na(df.cov.K["Menopause.Age",]) ]))
if( age.covariate ){
  lm.Menopause.Age.DNAm.Age <- glm(DNAm.Age.tmp ~ Menopause.Age + sample.age.tmp)
}else{
  lm.Menopause.Age.DNAm.Age <- glm(DNAm.Age.tmp ~ Menopause.Age)
}

##### VDYes DNAm Age #####
DNAm.Age.tmp      <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["VDYes",]) ]))
VDYes             <- as.numeric(as.vector(df.cov.K["VDYes",    !is.na(df.cov.K["VDYes",]) ]))
sample.age.tmp    <- as.numeric(as.vector(df.cov.K["Age",      !is.na(df.cov.K["VDYes",]) ]))
if( age.covariate ){
  lm.VDYes.DNAm.Age <- glm(DNAm.Age.tmp ~ VDYes + sample.age.tmp)
}else{
  lm.VDYes.DNAm.Age <- glm(DNAm.Age.tmp ~ VDYes)
}

##### Location DNAm Age #####
Location.Rural   <- as.numeric(as.vector(df.cov.K["LocationRural", (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
Loacation.Urban  <- as.numeric(as.vector(df.cov.K["LocationUrban", (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
DNAm.Age.Rural   <- as.numeric(as.vector(df.cov.K["DNAm.Age",      (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
DNAm.Age.Urban   <- as.numeric(as.vector(df.cov.K["DNAm.Age",      (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
sample.age.Rural <- as.numeric(as.vector(df.cov.K["Age",           (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
sample.age.Urban <- as.numeric(as.vector(df.cov.K["Age",           (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
Location.Rural[Location.Rural == 1] <- 0

DNAm.Age.tmp         <- c(DNAm.Age.Rural, DNAm.Age.Urban)
Location             <- c(Location.Rural, Loacation.Urban)         
sample.age.tmp       <- c(sample.age.Rural, sample.age.Urban) 
if( age.covariate ){
  lm.Location.DNAm.Age <- glm(DNAm.Age.tmp ~ Location + sample.age.tmp)
}else{
  lm.Location.DNAm.Age <- glm(DNAm.Age.tmp ~ Location)
}

##### Menopause DNAm Age #####
MenopausePost            <- as.numeric(as.vector(df.cov.K["MenopausePost-menopausal", (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
MenopausePre             <- as.numeric(as.vector(df.cov.K["MenopausePre-menopausal",  (df.cov.K["MenopausePre-menopausal",]  == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
DNAm.Age.MenopausePost   <- as.numeric(as.vector(df.cov.K["DNAm.Age",                 (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
DNAm.Age.MenopausePre    <- as.numeric(as.vector(df.cov.K["DNAm.Age",                 (df.cov.K["MenopausePre-menopausal",]  == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
sample.age.MenopausePost <- as.numeric(as.vector(df.cov.K["Age",                      (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
sample.age.MenopausePre  <- as.numeric(as.vector(df.cov.K["Age",                      (df.cov.K["MenopausePre-menopausal",]  == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
MenopausePre[MenopausePre == 1] <- 0

DNAm.Age.tmp   <- as.numeric( c(DNAm.Age.MenopausePre, DNAm.Age.MenopausePost) )
Menopause      <- as.numeric( c(MenopausePre, MenopausePost) )
sample.age.tmp <- as.numeric( c(sample.age.MenopausePre, sample.age.MenopausePost) )
if( age.covariate ){
  lm.Menopause.DNAm.Age <- glm(DNAm.Age.tmp ~ Menopause + sample.age.tmp)
}else{
  lm.Menopause.DNAm.Age <- glm(DNAm.Age.tmp ~ Menopause)
}

##########################################################################
# DNAm Age Residual Univariate Analysis
##########################################################################

##### Race DNAm Age Residual #####
DNAm.Age.Residual.tmp     <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",  !is.na(df.cov.K["RaceWhite",]) ]))
Race                      <- as.numeric(as.vector(df.cov.K["RaceWhite",          !is.na(df.cov.K["RaceWhite",]) ]))
sample.age.tmp            <- as.numeric(as.vector(df.cov.K["Age",                !is.na(df.cov.K["RaceWhite",]) ]))
if( age.covariate ){
  lm.Race.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Race + sample.age.tmp)
}else{
  lm.Race.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Race)
}

##### BMI DNAm Age Residual #####
DNAm.Age.Residual.tmp    <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["BMI",]) ]))
BMI                      <- as.numeric(as.vector(df.cov.K["BMI",               !is.na(df.cov.K["BMI",]) ]))
sample.age.tmp           <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["BMI",]) ]))
if( age.covariate ){
  lm.BMI.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ BMI + sample.age.tmp)
}else{
  lm.BMI.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ BMI)
}

##### Smoking DNAm Age Residual #####
DNAm.Age.Residual.tmp         <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["SmokingYes",]) ]))
SmokingYes                    <- as.numeric(as.vector(df.cov.K["SmokingYes",        !is.na(df.cov.K["SmokingYes",]) ]))
sample.age.tmp                <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["SmokingYes",]) ]))
if( age.covariate ){
  lm.SmokingYes.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ SmokingYes + sample.age.tmp)
}else{
  lm.SmokingYes.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ SmokingYes)
}

##### Drinking DNAm Age Residual #####
DNAm.Age.Residual.tmp         <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["DrinkingYes",]) ]))
Drinking                      <- as.numeric(as.vector(df.cov.K["DrinkingYes",       !is.na(df.cov.K["DrinkingYes",]) ]))
sample.age.tmp                <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["DrinkingYes",]) ]))
if( age.covariate ){
  lm.Drinking.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Drinking + sample.age.tmp)
}else{
  lm.Drinking.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Drinking)
}

##### Menarche DNAm Age Residual #####
DNAm.Age.Residual.tmp         <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Menarche",]) ]))
Menarche                      <- as.numeric(as.vector(df.cov.K["Menarche",          !is.na(df.cov.K["Menarche",]) ]))
sample.age.tmp                <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["Menarche",]) ]))
if( age.covariate ){
  lm.Menarche.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Menarche + sample.age.tmp)
}else{
  lm.Menarche.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Menarche)
}

##### Been.PregnantYes DNAm Age Residual #####
DNAm.Age.Residual.tmp                 <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Been.PregnantYes",]) ]))
Been.PregnantYes                      <- as.numeric(as.vector(df.cov.K["Been.PregnantYes",  !is.na(df.cov.K["Been.PregnantYes",]) ]))
sample.age.tmp                        <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["Been.PregnantYes",]) ]))
if( age.covariate ){
  lm.Been.PregnantYes.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Been.PregnantYes + sample.age.tmp)
}else{
  lm.Been.PregnantYes.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Been.PregnantYes)
}

##### Times.Pregnant DNAm Age Residual #####
DNAm.Age.Residual.tmp               <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",  !is.na(df.cov.K["Times.Pregnant",]) ]))
Times.Pregnant                      <- as.numeric(as.vector(df.cov.K["Times.Pregnant",     !is.na(df.cov.K["Times.Pregnant",]) ]))
sample.age.tmp                      <- as.numeric(as.vector(df.cov.K["Age",                !is.na(df.cov.K["Times.Pregnant",]) ]))
if( age.covariate ){
  lm.Times.Pregnant.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Times.Pregnant + sample.age.tmp)
}else{
  lm.Times.Pregnant.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Times.Pregnant)
}

##### Parity DNAm Age Residual #####
DNAm.Age.Residual.tmp       <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Parity",]) ]))
Parity                      <- as.numeric(as.vector(df.cov.K["Parity",            !is.na(df.cov.K["Parity",]) ]))
sample.age.tmp              <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["Parity",]) ]))
if( age.covariate ){
  lm.Parity.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Parity + sample.age.tmp)
}else{
  lm.Parity.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Parity)
}

##### Age.FB DNAm Age Residual #####
DNAm.Age.Residual.tmp       <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Age.FB",]) ]))
Age.FB                      <- as.numeric(as.vector(df.cov.K["Age.FB",            !is.na(df.cov.K["Age.FB",]) ]))
sample.age.tmp              <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["Age.FB",]) ]))
if( age.covariate ){
  lm.Age.FB.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Age.FB + sample.age.tmp)
}else{
  lm.Age.FB.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Age.FB)
}

##### Menopause.Age DNAm Age Residual #####
DNAm.Age.Residual.tmp              <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Menopause.Age",]) ]))
Menopause.Age                      <- as.numeric(as.vector(df.cov.K["Menopause.Age",     !is.na(df.cov.K["Menopause.Age",]) ]))
sample.age.tmp                     <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["Menopause.Age",]) ]))
if( age.covariate ){
  lm.Menopause.Age.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Menopause.Age + sample.age.tmp)
}else{
  lm.Menopause.Age.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Menopause.Age)
}

##### VDYes DNAm Age Residual #####
DNAm.Age.Residual.tmp      <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["VDYes",]) ]))
VDYes                      <- as.numeric(as.vector(df.cov.K["VDYes",             !is.na(df.cov.K["VDYes",]) ]))
sample.age.tmp             <- as.numeric(as.vector(df.cov.K["Age",               !is.na(df.cov.K["VDYes",]) ]))
if( age.covariate ){
  lm.VDYes.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ VDYes + sample.age.tmp)
}else{
  lm.VDYes.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ VDYes)
}

##### Location DNAm Age Residual #####
Location.Rural            <- as.numeric(as.vector(df.cov.K["LocationRural",     (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
Loacation.Urban           <- as.numeric(as.vector(df.cov.K["LocationUrban",     (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
DNAm.Age.Residual.Rural   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
DNAm.Age.Residual.Urban   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
sample.age.Rural          <- as.numeric(as.vector(df.cov.K["Age",               (df.cov.K["LocationRural",] == 1 & !is.na(df.cov.K["LocationRural",]))]))
sample.age.Urban          <- as.numeric(as.vector(df.cov.K["Age",               (df.cov.K["LocationUrban",] == 1 & !is.na(df.cov.K["LocationUrban",]))]))
Location.Rural[Location.Rural == 1] <- 0

DNAm.Age.Residual.tmp         <- c(DNAm.Age.Residual.Rural, DNAm.Age.Residual.Urban)
Location                      <- c(Location.Rural, Loacation.Urban)         
sample.age.tmp                <- c(sample.age.Rural, sample.age.Urban) 
if( age.covariate ){
  lm.Location.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Location + sample.age.tmp)
}else{
  lm.Location.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Location)
}

##### Menopause DNAm Age Residual #####
MenopausePost                     <- as.numeric(as.vector(df.cov.K["MenopausePost-menopausal", (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
MenopausePre                      <- as.numeric(as.vector(df.cov.K["MenopausePre-menopausal",  (df.cov.K["MenopausePre-menopausal",]  == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
DNAm.Age.Residual.MenopausePost   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",        (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
DNAm.Age.Residual.MenopausePre    <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",        (df.cov.K["MenopausePre-menopausal",]  == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
sample.age.MenopausePost          <- as.numeric(as.vector(df.cov.K["Age",                      (df.cov.K["MenopausePost-menopausal",] == 1 & !is.na(df.cov.K["MenopausePost-menopausal",]))]))
sample.age.MenopausePre           <- as.numeric(as.vector(df.cov.K["Age",                      (df.cov.K["MenopausePre-menopausal",]  == 1 & !is.na(df.cov.K["MenopausePre-menopausal",]))]))
MenopausePre[MenopausePre == 1] <- 0

DNAm.Age.Residual.tmp          <- as.numeric( c(DNAm.Age.Residual.MenopausePre, DNAm.Age.Residual.MenopausePost) )
Menopause                      <- as.numeric( c(MenopausePre, MenopausePost) )
sample.age.tmp                 <- as.numeric( c(sample.age.MenopausePre, sample.age.MenopausePost) )
if( age.covariate ){
  lm.Menopause.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Menopause + sample.age.tmp)
}else{
  lm.Menopause.DNAm.Age.Residual <- glm(DNAm.Age.Residual.tmp ~ Menopause)
}


##########################################################################
# Univariate Analysis - Race
##########################################################################
residual.K <- result.K - sample.ages.K 
df.cov.K["DNAm.Age.Residual",] <- residual.K

##### Correlate DNAm Age and Residual to Race #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["RaceWhite",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["RaceWhite",]) ]))
Race                <- as.numeric(as.vector(df.cov.K["RaceWhite", !is.na(df.cov.K["RaceWhite",]) ]))

##### Plot Correlations #####
DNAm.Age.White            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 1)]))
DNAm.Age.AfrAmer          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 0)]))
DNAm.Age.Residual.White   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 1)]))
DNAm.Age.Residual.AfrAmer <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["RaceWhite",]) & df.cov.K["RaceWhite",] == 0)]))

Median.DNAm.Age.0          <- median(DNAm.Age.AfrAmer)
Median.DNAm.Age.1          <- median(DNAm.Age.White)
Median.DNAm.Age.Residual.0 <- median(DNAm.Age.Residual.AfrAmer)
Median.DNAm.Age.Residual.1 <- median(DNAm.Age.Residual.White)

df1           <- data.frame(DNAm.Age.White)
df1$ttype     <- "White"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.AfrAmer)
df2$ttype     <- "African American"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.White)
df1$ttype     <- "White"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.AfrAmer)
df2$ttype     <- "African American"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)
                          
p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Race", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Race-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.AfrAmer
r.list[[2]] <- DNAm.Age.Residual.White
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Race", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Race-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

tiff( paste(model.dir, "FactorAssociation/Race-DNAmAgeResidual.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
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

##### Correlate DNAm Age and Residual to Location #####
Location.Rural[Location.Rural == 1] <- 0
DNAm.Age            <- as.numeric( c(DNAm.Age.Rural, DNAm.Age.Urban) )
DNAm.Age.Residual   <- as.numeric( c(DNAm.Age.Residual.Rural, DNAm.Age.Residual.Urban) )
Location            <- as.numeric( c(Location.Rural, Loacation.Urban) )

df1           <- data.frame(DNAm.Age.Rural)
df1$ttype     <- "Rural"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.Urban)
df2$ttype     <- "Urban"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.Rural)
df1$ttype     <- "Rural"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.Urban)
df2$ttype     <- "Urban"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)

p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Location", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Location-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.Rural
r.list[[2]] <- DNAm.Age.Residual.Urban
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Location", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Location-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

tiff( paste(model.dir, "FactorAssociation/Location-DNAmAgeResidual.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()


##########################################################################
# Univariate Analysis - BMI
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["BMI",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["BMI",]) ]))
BMI                 <- as.numeric(as.vector(df.cov.K["BMI", !is.na(df.cov.K["BMI",]) ]))

##### Plot Correlations #####
df.tmp <- data.frame(BMI, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=BMI, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "BMI") + 
  labs(y = "Epigenetic Age [Years]") + 
  geom_abline(
    slope = lm(DNAm.Age ~ BMI)$coefficients[[2]], 
    intercept = lm(DNAm.Age ~ BMI)$coefficients[[1]],
    linetype = "dotted",
    color = "red"
  ) +
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

png( paste(model.dir, "FactorAssociation/BMI-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

p <- ggplot(df.tmp, aes(x = BMI, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "BMI") + 
  labs(y = "Epigenetic Age Acceleration [Years]") + 
  geom_abline(
    slope = lm(DNAm.Age.Residual ~ BMI)$coefficients[[2]], 
    intercept = lm(DNAm.Age.Residual ~ BMI)$coefficients[[1]],
    linetype = "dotted",
    color = "red"
  ) + 
  theme_bw(base_size = 18) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/BMI-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

tiff( paste(model.dir, "FactorAssociation/BMI-DNAmAgeResidual.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()

##########################################################################
# Univariate Analysis - Current Drinker
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["DrinkingYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["DrinkingYes",]) ]))
Drinking            <- as.numeric(as.vector(df.cov.K["DrinkingYes", !is.na(df.cov.K["DrinkingYes",]) ]))

##### Plot Correlations #####
DNAm.Age.DrinkingYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 1)]))
DNAm.Age.DrinkingNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 0)]))
DNAm.Age.Residual.DrinkingYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 1)]))
DNAm.Age.Residual.DrinkingNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["DrinkingYes",]) & df.cov.K["DrinkingYes",] == 0)]))

df1           <- data.frame(DNAm.Age.DrinkingNo)
df1$ttype     <- "Current Alcohol Use - No"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.DrinkingYes)
df2$ttype     <- "Current Alcohol Use - Yes"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.DrinkingNo)
df1$ttype     <- "Current Alcohol Use - No"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.DrinkingYes)
df2$ttype     <- "Current Alcohol Use - Yes"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)

p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Alcohol Use", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Drinking-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.DrinkingNo
r.list[[2]] <- DNAm.Age.Residual.DrinkingYes
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Alcohol Use", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 17,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Drinking-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

##########################################################################
# Univariate Analysis - Current Drinker
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["SmokingYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["SmokingYes",]) ]))
Smoking             <- as.numeric(as.vector(df.cov.K["SmokingYes", !is.na(df.cov.K["SmokingYes",]) ]))

##### Plot Correlations #####
DNAm.Age.SmokingYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 1)]))
DNAm.Age.SmokingNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 0)]))
DNAm.Age.Residual.SmokingYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 1)]))
DNAm.Age.Residual.SmokingNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (!is.na(df.cov.K["SmokingYes",]) & df.cov.K["SmokingYes",] == 0)]))

df1           <- data.frame(DNAm.Age.SmokingNo)
df1$ttype     <- "Current Tobacco Use - No"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.SmokingYes)
df2$ttype     <- "Current Tobacco Use - Yes"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.SmokingNo)
df1$ttype     <- "Current Tobacco Use - No"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.SmokingYes)
df2$ttype     <- "Current Tobacco Use - Yes"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)

p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Tobacco Use", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Smoking-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.SmokingNo
r.list[[2]] <- DNAm.Age.Residual.SmokingYes
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Alcohol Use", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 17,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/Smoking-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()


##########################################################################
# Univariate Analysis - Age at Menarche
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))
Menarche            <- as.numeric(as.vector(df.cov.K["Menarche", (df.cov.K["Menarche",] >= 0 & !is.na(df.cov.K["Menarche",])) ]))

##### Plot Correlations #####
df.tmp <- data.frame(Menarche, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Menarche, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at Menarche") + 
  labs(y = "DNAm Age") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/MenarcheAge-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Menarche, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at Menarche") +
  labs(y = "Epigenetic Age Acceleration [Years]") + 
  geom_abline(
    slope = lm(DNAm.Age.Residual ~ Menarche)$coefficients[[2]], 
    intercept = lm(DNAm.Age.Residual ~ Menarche)$coefficients[[1]],
    linetype = "dotted",
    color = "red"
  ) + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/MenarcheAge-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

tiff( paste(model.dir, "FactorAssociation/MenarcheAge-DNAmAgeResidual.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()


##########################################################################
# Univariate Analysis - Have you been pregnant
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", !is.na(df.cov.K["Been.PregnantYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", !is.na(df.cov.K["Been.PregnantYes",]) ]))
Been.PregnantYes    <- as.numeric(as.vector(df.cov.K["Been.PregnantYes", !is.na(df.cov.K["Been.PregnantYes",]) ]))

##### Plot Correlations #####
DNAm.Age.Been.PregnantYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Been.PregnantYes",] == 1 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Been.PregnantNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Been.PregnantYes",] == 0 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Residual.Been.PregnantYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Been.PregnantYes",] == 1 & !is.na(df.cov.K["Been.PregnantYes",]))]))
DNAm.Age.Residual.Been.PregnantNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Been.PregnantYes",] == 0 & !is.na(df.cov.K["Been.PregnantYes",]))]))

df1           <- data.frame(DNAm.Age.Been.PregnantNo)
df1$ttype     <- "No"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.Been.PregnantYes)
df2$ttype     <- "Yes"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.Been.PregnantNo)
df1$ttype     <- "No"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.Been.PregnantYes)
df2$ttype     <- "Yes"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)

p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Have You Been Pregnant?", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/HaveYouBeenPregnant-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.Been.PregnantNo
r.list[[2]] <- DNAm.Age.Residual.Been.PregnantYes
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Have You Been Pregnant?", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 17,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/HaveYouBeenPregnant-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()


##########################################################################
# Univariate Analysis - Times Pregnant
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))
Times.Pregnant      <- as.numeric(as.vector(df.cov.K["Times.Pregnant", (df.cov.K["Times.Pregnant",] >= 0 & !is.na(df.cov.K["Times.Pregnant",])) ]))

##### Plot Correlations #####
df.tmp <- data.frame(Times.Pregnant, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Times.Pregnant, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Times Pregnant") + 
  labs(y = "DNAm Age") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/TimesPregnant-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Times.Pregnant, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Times Pregnant") +
  labs(y = "DNAm Age Acceleration") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/TimesPregnant-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()


##########################################################################
# Univariate Analysis - Parity
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))
Parity              <- as.numeric(as.vector(df.cov.K["Parity", (df.cov.K["Parity",] >= 0 & !is.na(df.cov.K["Parity",])) ]))

##### Plot Correlations #####
df.tmp <- data.frame(Parity, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Parity, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Parity") + 
  labs(y = "DNAm Age") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/Parity-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Parity, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Parity") +
  labs(y = "DNAm Age Acceleration") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/Parity-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()


##########################################################################
# Univariate Analysis - Age at first birth
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))
Age.FB              <- as.numeric(as.vector(df.cov.K["Age.FB", (df.cov.K["Age.FB",] >= 0 & !is.na(df.cov.K["Age.FB",])) ]))

##### Plot Correlations #####
df.tmp <- data.frame(Age.FB, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Age.FB, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at First Birth [Years]") + 
  labs(y = "Epigenetic Age [Years]") + 
  geom_abline(
    slope = lm(DNAm.Age ~ Age.FB)$coefficients[[2]], 
    intercept = lm(DNAm.Age ~ Age.FB)$coefficients[[1]],
    linetype = "dotted",
    color = "red"
  ) +
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/AgeAtFirstBirth-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Age.FB, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at First Birth [Years]") +
  labs(y = "Epigenetic Age Acceleration [Years]") + 
  geom_abline(
    slope = lm(DNAm.Age.Residual ~ Age.FB)$coefficients[[2]], 
    intercept = lm(DNAm.Age.Residual ~ Age.FB)$coefficients[[1]],
    linetype = "dotted",
    color = "red"
  ) +
  theme_bw(base_size = 18) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/AgeAtFirstBirth-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

tiff( paste(model.dir, "FactorAssociation/AgeAtFirstBirth-DNAmAgeResidual.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
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

##### Correlate DNAm Age and Residual to Menopause #####
MenopausePre[MenopausePre == 1] <- 0
DNAm.Age            <- as.numeric( c(DNAm.Age.MenopausePre, DNAm.Age.MenopausePost) )
DNAm.Age.Residual   <- as.numeric( c(DNAm.Age.Residual.MenopausePre, DNAm.Age.Residual.MenopausePost) )
Menopause           <- as.numeric( c(MenopausePre, MenopausePost) )

df1           <- data.frame(DNAm.Age.MenopausePre)
df1$ttype     <- "Pre-menopausal"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.MenopausePost)
df2$ttype     <- "Post-menopausal"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.MenopausePre)
df1$ttype     <- "Pre-menopausal"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.MenopausePost)
df2$ttype     <- "Post-menopausal"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)

df.tmp.Age$ttype <- factor( df.tmp.Age$ttype, levels = unique(as.character(df.tmp.Age$ttype)) )
df.tmp.Res$ttype <- factor( df.tmp.Res$ttype, levels = unique(as.character(df.tmp.Res$ttype)) )

p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Menopause Status", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/MenopauseStatus-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.MenopausePre
r.list[[2]] <- DNAm.Age.Residual.MenopausePost
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Menopause Status", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 17,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/MenopauseStatus-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()

tiff( paste(model.dir, "FactorAssociation/MenopauseStatus-DNAmAgeResidual.tiff", sep = ''), width = 2100, height = 2100, units = "px", res = 300)
p
dev.off()


##########################################################################
# Univariate Analysis - Age at Menopause
##########################################################################

##### Correlate DNAm Age and Residual to Factor #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", (df.cov.K["Menopause.Age",] >= 0 & !is.na(df.cov.K["Menopause.Age",])) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", (df.cov.K["Menopause.Age",] >= 0 & !is.na(df.cov.K["Menopause.Age",])) ]))
Age.Menopause       <- as.numeric(as.vector(df.cov.K["Menopause.Age", (df.cov.K["Menopause.Age",] >= 0 & !is.na(df.cov.K["Menopause.Age",])) ]))

##### Plot Correlations #####
df.tmp <- data.frame(Age.Menopause, DNAm.Age, DNAm.Age.Residual)

p <- ggplot(df.tmp, aes(x=Age.Menopause, y=DNAm.Age)) + 
  geom_point() + 
  labs(x = "Age at Menopause") +
  labs(y = "DNAm Age") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/MenopauseAge-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

p <- ggplot(df.tmp, aes(x = Age.Menopause, y = DNAm.Age.Residual)) + 
  geom_point() + 
  labs(x = "Age at Menopause") +
  labs(y = "DNAm Age Acceleration") + 
  theme_bw(base_size = 22) + 
  theme(
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ); p

png( paste(model.dir, "FactorAssociation/MenopauseAge-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()


##########################################################################
# Univariate Analysis - Vitamin Use
##########################################################################

##### Correlate DNAm Age and Residual to Race #####
DNAm.Age            <- as.numeric(as.vector(df.cov.K["DNAm.Age", complete.cases(df.cov.K["VDYes",]) ]))
DNAm.Age.Residual   <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", complete.cases(df.cov.K["VDYes",]) ]))
VD                  <- as.numeric(as.vector(df.cov.K["VDYes", complete.cases(df.cov.K["VDYes",]) ]))

##### Plot Correlations #####
DNAm.Age.VDYes          <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["VDYes",] == 1]))
DNAm.Age.VDNo           <- as.numeric(as.vector(df.cov.K["DNAm.Age", df.cov.K["VDYes",] == 0]))
DNAm.Age.Residual.VDYes <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["VDYes",] == 1]))
DNAm.Age.Residual.VDNo  <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual", df.cov.K["VDYes",] == 0]))

df1           <- data.frame(DNAm.Age.VDNo)
df1$ttype     <- "No"
colnames(df1) <- c("Age", "ttype")
df2           <- data.frame(DNAm.Age.VDYes)
df2$ttype     <- "Yes"
colnames(df2) <- c("Age", "ttype")
df.tmp.Age    <- rbind(df1, df2)

df1           <- data.frame(DNAm.Age.Residual.VDNo)
df1$ttype     <- "No"
colnames(df1) <- c("Res", "ttype")
df2           <- data.frame(DNAm.Age.Residual.VDYes)
df2$ttype     <- "Yes"
colnames(df2) <- c("Res", "ttype")
df.tmp.Res    <- rbind(df1, df2)

p <- DNAmAgeRF.box.plot(
  df.tmp.Age, 
  x.label = "Multivitamin Use", 
  y.label = "Epigenetic Age [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 18,
  show.leg = FALSE
); p

png( paste(model.dir, "FactorAssociation/VDUse-DNAmAge.png", sep = ''), width = 500, height = 500 )
p
dev.off()

r.list <- list()
r.list[[1]] <- DNAm.Age.Residual.VDNo
r.list[[2]] <- DNAm.Age.Residual.VDYes
p <- DNAmAgeResRF.box.plot(
  df.tmp.Res, 
  residuals = r.list,
  x.label = "Multivitamin Use", 
  y.label = "Epigenetic Age Acceleration [Years]",
  title = "",
  width = 0.6,
  leg.x = 0.77,
  leg.y = 0.93,
  text.size = 17,
  show.leg = FALSE
); p
png( paste(model.dir, "FactorAssociation/VDUse-DNAmAgeResidual.png", sep = ''), width = 500, height = 500 )
p
dev.off()


##########################################################################
# Univariate Analysis Summary
##########################################################################
options(digits=4)
residual.K <- result.K - sample.ages.K 
if(model.residual){
  residual.K <- as.numeric(as.vector(lm(result.K ~ sample.ages.K)$residual))
}
df.cov.K["DNAm.Age.Residual",] <- residual.K

##### DNAm Age Univariate #####
df.Race             <- data.frame(summary(lm.Race.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Location         <- data.frame(summary(lm.Location.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.BMI              <- data.frame(summary(lm.BMI.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.SmokingYes       <- data.frame(summary(lm.SmokingYes.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Drinking         <- data.frame(summary(lm.Drinking.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Menarche         <- data.frame(summary(lm.Menarche.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Been.PregnantYes <- data.frame(summary(lm.Been.PregnantYes.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Times.Pregnant   <- data.frame(summary(lm.Times.Pregnant.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Parity           <- data.frame(summary(lm.Parity.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Age.FB           <- data.frame(summary(lm.Age.FB.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Menopause        <- data.frame(summary(lm.Menopause.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Menopause.Age    <- data.frame(summary(lm.Menopause.Age.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.VDYes            <- data.frame(summary(lm.VDYes.DNAm.Age)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])

df.Race             <- df.Race[-c(1),]
df.Location         <- df.Location[-c(1),]
df.BMI              <- df.BMI[-c(1),]
df.SmokingYes       <- df.SmokingYes[-c(1),]
df.Drinking         <- df.Drinking[-c(1),]
df.Menarche         <- df.Menarche[-c(1),]
df.Been.PregnantYes <- df.Been.PregnantYes[-c(1),]
df.Times.Pregnant   <- df.Times.Pregnant[-c(1),]
df.Parity           <- df.Parity[-c(1),]
df.Age.FB           <- df.Age.FB[-c(1),]
df.Menopause        <- df.Menopause[-c(1),]
df.Menopause.Age    <- df.Menopause.Age[-c(1),]
df.VDYes            <- df.VDYes[-c(1),]

if( age.covariate ){
  df.Race             <- df.Race[-c(2),]
  df.Location         <- df.Location[-c(2),]
  df.BMI              <- df.BMI[-c(2),]
  df.SmokingYes       <- df.SmokingYes[-c(2),]
  df.Drinking         <- df.Drinking[-c(2),]
  df.Menarche         <- df.Menarche[-c(2),]
  df.Been.PregnantYes <- df.Been.PregnantYes[-c(2),]
  df.Times.Pregnant   <- df.Times.Pregnant[-c(2),]
  df.Parity           <- df.Parity[-c(2),]
  df.Age.FB           <- df.Age.FB[-c(2),]
  df.Menopause        <- df.Menopause[-c(2),]
  df.Menopause.Age    <- df.Menopause.Age[-c(2),]
  df.VDYes            <- df.VDYes[-c(2),]
  
  DNAm.Age   <- as.numeric(as.vector(df.cov.K["DNAm.Age",]))
  sample.age <- as.numeric(as.vector(df.cov.K["Age",]))
  m          <- glm(DNAm.Age ~ sample.age)
  df.age     <- data.frame(summary(m)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
  df.age     <- df.age[-c(1,3),]
  
  
  df.sum <- rbind(df.Race, df.BMI, df.Drinking, df.SmokingYes, df.Menarche, 
                  df.Times.Pregnant, df.Parity, df.Age.FB, df.VDYes, df.Menopause, 
                  df.Location, df.age)
  
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)", "Chronological age")
}else{
  df.sum <- rbind(df.Race, df.BMI, df.Drinking, df.SmokingYes, df.Menarche, 
                  df.Times.Pregnant, df.Parity, df.Age.FB, df.VDYes, df.Menopause, 
                  df.Location)
  
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)")
}

colnames(df.sum) <- c("Coef", "Std. Error", "pval")
stargazer(df.sum, summary = FALSE)

##### DNAm Age Residual Univariate #####
df.Race             <- data.frame(summary(lm.Race.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Location         <- data.frame(summary(lm.Location.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.BMI              <- data.frame(summary(lm.BMI.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.SmokingYes       <- data.frame(summary(lm.SmokingYes.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Drinking         <- data.frame(summary(lm.Drinking.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Menarche         <- data.frame(summary(lm.Menarche.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Been.PregnantYes <- data.frame(summary(lm.Been.PregnantYes.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Times.Pregnant   <- data.frame(summary(lm.Times.Pregnant.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Parity           <- data.frame(summary(lm.Parity.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Age.FB           <- data.frame(summary(lm.Age.FB.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Menopause        <- data.frame(summary(lm.Menopause.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.Menopause.Age    <- data.frame(summary(lm.Menopause.Age.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
df.VDYes            <- data.frame(summary(lm.VDYes.DNAm.Age.Residual)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])

df.Race             <- df.Race[-c(1),]
df.Location         <- df.Location[-c(1),]
df.BMI              <- df.BMI[-c(1),]
df.SmokingYes       <- df.SmokingYes[-c(1),]
df.Drinking         <- df.Drinking[-c(1),]
df.Menarche         <- df.Menarche[-c(1),]
df.Been.PregnantYes <- df.Been.PregnantYes[-c(1),]
df.Times.Pregnant   <- df.Times.Pregnant[-c(1),]
df.Parity           <- df.Parity[-c(1),]
df.Age.FB           <- df.Age.FB[-c(1),]
df.Menopause        <- df.Menopause[-c(1),]
df.Menopause.Age    <- df.Menopause.Age[-c(1),]
df.VDYes            <- df.VDYes[-c(1),]

if( age.covariate ){
  df.Race             <- df.Race[-c(2),]
  df.Location         <- df.Location[-c(2),]
  df.BMI              <- df.BMI[-c(2),]
  df.SmokingYes       <- df.SmokingYes[-c(2),]
  df.Drinking         <- df.Drinking[-c(2),]
  df.Menarche         <- df.Menarche[-c(2),]
  df.Been.PregnantYes <- df.Been.PregnantYes[-c(2),]
  df.Times.Pregnant   <- df.Times.Pregnant[-c(2),]
  df.Parity           <- df.Parity[-c(2),]
  df.Age.FB           <- df.Age.FB[-c(2),]
  df.Menopause        <- df.Menopause[-c(2),]
  df.Menopause.Age    <- df.Menopause.Age[-c(2),]
  df.VDYes            <- df.VDYes[-c(2),]
  
  DNAm.Age.Residual <- as.numeric(as.vector(df.cov.K["DNAm.Age.Residual",]))
  sample.age        <- as.numeric(as.vector(df.cov.K["Age",]))
  m                 <- glm(DNAm.Age.Residual ~ sample.age)
  df.age            <- data.frame(summary(m)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
  df.age            <- df.age[-c(1,3),]
  
  
  df.sum <- rbind(df.Race, df.BMI, df.Drinking, df.SmokingYes, df.Menarche, 
                  df.Times.Pregnant, df.Parity, df.Age.FB, df.VDYes, df.Menopause, 
                  df.Location, df.age)
  
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)", "Chronological age")
}else{
  df.sum <- rbind(df.Race, df.BMI, df.Drinking, df.SmokingYes, df.Menarche, 
                  df.Times.Pregnant, df.Parity, df.Age.FB, df.VDYes, df.Menopause, 
                  df.Location)
  
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)")
}
summary(m)
colnames(df.sum) <- c("Coef", "Std. Error", "pval")
stargazer(df.sum, summary = FALSE)


##########################################################################
# Multivariate Analysis Summary
##########################################################################

##### DNAm Age #####
df.sum <- data.frame(summary(lm.DNAm.Age)$coefficients[,c("Estimate","Pr(>|t|)")])
df.sum <- df.sum[-c(1),]
colnames(df.sum) <- c("Coef", "pval")
df.sum
if(age.covariate){
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)", "Chronological age")
}else{
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)")
}

stargazer(df.sum, summary = FALSE)

##### DNAm Age Residual #####
df.sum <- data.frame(summary(lm.DNAm.Age.Residual)$coefficients[,c("Estimate","Pr(>|t|)")])
df.sum <- df.sum[-c(1),]
colnames(df.sum) <- c("Coef", "pval")
if(age.covariate){
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)", "Chronological age")
}else{
  rownames(df.sum) <- c("Race (White vs African American)", "BMI", "Current alcohol consumption (yes vs no)",
                        "Current Tobacco Consumption (yes vs no)", "Age at menarche", "Times pregnant", "Parity", 
                        "Age at first birth", "Multivitamin use (yes vs no)", "Menopause status (pre vs post)", 
                        "Location (urban vs rural)")
}

stargazer(df.sum, summary = FALSE)


###########################
mean(sample.age); sd(sample.age)
100.* length(Race[which(Race == 1)])/length(Race)
mean(BMI); sd(BMI)
100.* length(Drinking[which(Drinking == 1)])/length(Drinking)
100.* length(Smoking[which(Smoking == 1)])/length(Smoking)
mean(Menarche); sd(Menarche)
mean(Times.Pregnant); sd(Times.Pregnant)
mean(Parity); sd(Parity)
mean(Age.FB); sd(Age.FB)
100.* length(VD[which(VD == 1)])/length(VD)
100.* length(MenopausePre) / (length(MenopausePre) + length(MenopausePost) )
100.* length(Location.Urban) / (length(Location.Urban) + length(Location.Rural) )







