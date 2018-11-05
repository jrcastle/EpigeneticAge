setwd('/home/jrca253/EpigeneticAge/data')

# Cut CpG if it is missing in > missingness.fraction samples
# To ease computational load on large-sample datasets, the 
# maximum number of samples missing a particular CpG is 
# hardcoded at 10.
missingness.fraction = 0.01

# K => Normal
# N => Adjacent Normal
# T => Tumor
tissue.type <- "N"
DATADIR     <- '/home/jrca253/DATA/Truseq/'
#METHFILE    <- paste('meth_', tissue.type, '.txt', sep = '')
METHFILE    <- paste('meth_', tissue.type, '_AllCpGs', '.txt', sep = '')
COVFILE     <- paste('cov_',  tissue.type, '.txt', sep = '')


########################################################################################
# LOAD MASTER FILE
########################################################################################
#cov <- read.csv( paste(DATADIR, "For_Nan_DNA_Methylation_data_ids_001-624.csv", sep='') )
cov <- read.csv( paste(DATADIR, "MasterFile.csv", sep = '') )
cov$ID <- as.numeric(rownames(cov))
cov$Batch <- ifelse(cov$ID <= 48, 'HS_032',
                    ifelse(cov$ID <= 96, 'HS_049',
                           ifelse(cov$ID <= 144, 'HS_056',
                                  ifelse(cov$ID <= 192, 'HS_052',
                                         ifelse(cov$ID <= 240, 'HS_058',
                                                ifelse(cov$ID <= 288, 'HS_061',
                                                       ifelse(cov$ID <= 336, 'HS_066',
                                                              ifelse(cov$ID <= 384, 'HS_067',
                                                                     ifelse(cov$ID <= 432, 'HS_070',
                                                                            ifelse(cov$ID <= 480, 'HS_076',
                                                                                   ifelse(cov$ID <= 528, 'HS_087',
                                                                                          ifelse(cov$ID <= 576, 'HS_089',
                                                                                                 ifelse(cov$ID <= 624, 'HS_091', NA
                                                                                                       )
                                                                                                )
                                                                                         )
                                                                                  )
                                                                           )
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )
                   )


ID                <- as.vector(cov$ID)
cov$Race          <- as.character(cov$Race)
cov$Race          <- ifelse(cov$Race == 'AfricanAmerican', 'African American', cov$Race)
cov$dataframe.name <- paste(cov$Batch,cov$ID, sep = '_')


########################################################################################
# SELECT SAMPLES
########################################################################################
cov.selected <- cov[which(substr(cov$Sample.id,1,1)==tissue.type),]
cov.selected[cov.selected==''] <- NA


########################################################################################
# OMIT MISSING SAMPLES
########################################################################################
cov.selected <- cov.selected[which(cov.selected$ID<65 |cov.selected$ID>72),]
ID <- as.vector(cov.selected$ID)

dataframe.name = as.character(cov.selected[which(cov.selected$ID<=96),]$Sample.id)
dataframe.name = do.call(c, list(dataframe.name, as.character(cov.selected[which(cov.selected$ID>96),]$dataframe.name)))


########################################################################################
# COVARIATES
########################################################################################
Age            <- as.numeric(as.character(cov.selected$Age))
Race           <- as.character(cov.selected$Race)
Batch          <- as.character(cov.selected$Batch)
Sampleid       <- as.character(cov.selected$Sample.id)
                         
BMI            <- as.numeric(as.character(cov.selected$BMI))
Location       <- as.character(cov.selected$Location)
Smoking        <- as.character(cov.selected$Current.Smoker)
Cig.Pack.Years <- as.numeric(as.character(cov.selected$Cigarette.Pack.Years))
Drinking       <- as.character(cov.selected$Currently.Drink)
Menarche       <- as.numeric(as.character(cov.selected$Menarche))
Menopause      <- as.character(cov.selected$Menstrual.Status)
Menopause.Age  <- as.numeric(as.character(cov.selected$Age.at.menopause))
Been.Pregnant  <- as.character(cov.selected$Have.You.Been.Pregnant.)
Times.Pregnant <- as.numeric(as.character(cov.selected$How.Many.Times.))
Age.FB         <- as.numeric(as.character(cov.selected$Age.at.First.Birth))
Parity         <- as.numeric(as.character(cov.selected$Number.of.Live.Births))
VD             <- as.character(cov.selected$Vitaminuse)
Cancer.subtype <- as.character(cov.selected$Cancer.subtypes)
Cancer.grade   <- as.character(cov.selected$GradeBin)
Cancer.stage   <- as.character(cov.selected$Stage)

#Age.FB    <- ifelse(cov.selected$Age.at.First.Birth < 25, 1,
#                    ifelse(cov.selected$Age.at.First.Birth < 30, 2,
#                           ifelse(cov.selected$Age.at.First.Birth < 35, 3,4)
#			                    )
#	           )


########################################################################################
# MAKE COV FILE
########################################################################################
options(na.action='na.pass')
df.selected <- data.frame(Age, Race, Batch)
X           <- model.matrix(~Age + Race + Batch, df.selected)

if(tissue.type == "T"){
  df.selected <- data.frame(Age, Race, Batch, Cancer.subtype, Cancer.grade, Cancer.stage)
  X           <- model.matrix(~Age + Race + Batch + Cancer.subtype + Cancer.grade + Cancer.stage, df.selected)
}

if(tissue.type == "K"){
  df.selected <- data.frame(Age, Race, Batch, BMI, Smoking, Cig.Pack.Years, Drinking, Menarche, Menopause, Menopause.Age, Been.Pregnant, Times.Pregnant, Age.FB, Parity, VD, Location)
  X           <- model.matrix(~Age + Race + Batch + BMI + Smoking + Cig.Pack.Years + Drinking + Menarche + Menopause + Menopause.Age + Been.Pregnant + Times.Pregnant + Age.FB + Parity + VD + Location, df.selected)
}

X           <- t(X)
colnames(X) <- dataframe.name
X           <- X[-1,]

write.table(X, file = COVFILE, quote = F, col.names = NA, row.names = T, sep = '\t')


########################################################################################
# MAKE METH FILE, MERGE SAMPLES
########################################################################################

##### PROCESS THE FIRST SAMPLE #####
i=1
print( paste("Processing sample", i, 'of', length(ID), sep = ' ') )
BATCH    <- Batch[i]
RACE     <- Race[i]
DIR      <- paste(DATADIR, BATCH, "/", sep = '')
SAMPLEID <- ID[i]
NAME     <- paste(Batch[i], ID[i], sep = '_')

if( BATCH == 'HS_032'|| BATCH == 'HS_049' ){
  SAMPLEID <- Sampleid[i]
  NAME <- Sampleid[i]
}

FILENAME <- paste(DIR, SAMPLEID, '_inTruseq.csv', sep = '')

B.seq.all           <- read.csv(FILENAME, header = T)
B.seq.all           <- B.seq.all[ which(B.seq.all$N >= 10),]
B.seq.all$position  <- paste(B.seq.all$chr, B.seq.all$pos, sep = ':')
B.seq.all$SAMPLEID  <- B.seq.all$X / B.seq.all$N
colnames(B.seq.all) <- c('chr','pos','N','X','position', NAME)
B.seq.all           <- B.seq.all[ c('position', NAME) ]


##### PROCESS AND MERGE REMAINING SAMPLES #####
for( i in 2:length(ID) ){
  print( paste("Processing sample", i, 'of', length(ID), sep = ' ') )

  BATCH    <- Batch[i]
  DIR      <- paste(DATADIR, BATCH, "/", sep = '')
  SAMPLEID <- ID[i]
  NAME     <- paste(Batch[i], ID[i], sep = '_')

  if(BATCH == 'HS_032'|| BATCH == 'HS_049' ){
    SAMPLEID <- Sampleid[i]
    NAME <- Sampleid[i]
  }
  
  FILENAME <- paste(DIR, SAMPLEID, '_inTruseq.csv', sep = '')
  
  tmp           <- read.csv(FILENAME, header = T)
  tmp           <- tmp[ which(tmp$N >= 10),]
  tmp$position  <- paste(tmp$chr, tmp$pos, sep = ':')
  tmp$SAMPLEID  <- tmp$X / tmp$N
  colnames(tmp) <- c('chr','pos','N','X','position', NAME)
  tmp           <- tmp[ c('position', NAME) ]
  
  B.seq.all <- Reduce(
    function(x, y) merge(x, y, by = c('position'), all=T),
    list(B.seq.all, tmp)
  )
  rm(tmp)
        
}

##### NA CUTS #####
na.cut <- as.integer( missingness.fraction * length(ID) )
if( na.cut > 10 ) na.cut = 10;
if( na.cut == 0 ) na.cut = 1;

#B.seq.all <- B.seq.all[rowSums(is.na(B.seq.all)) <= na.cut,]

write.table(B.seq.all, file = METHFILE, quote = F, append = FALSE, sep = "\t", row.names=FALSE)
rm(B.seq.all)
