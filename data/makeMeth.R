setwd('/home/jrca253/EpigeneticAge/data')

# K => Normal
# N => Adjacent Normal
# T => Tumor
tissue.type <- "N"
DATADIR     <- '/home/jrca253/DATA/Truseq/'
METHFILE    <- 'meth_N_lt10_missing.txt'
COVFILE     <- 'cov_N.txt'


########################################################################################
# LOAD MASTER FILE
########################################################################################
cov <- read.csv( paste(DATADIR, "For_Nan_DNA_Methylation_data_ids_001-624.csv", sep='') )
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
Age       <- as.numeric(as.character(cov.selected$Age))
Race      <- as.character(cov.selected$Race)
Batch     <- as.character(cov.selected$Batch)
Sampleid  <- as.character(cov.selected$Sample.id)
                         
BMI       <- as.numeric(as.character(cov.selected$BMI))
Smoking   <- as.character(cov.selected$Current.Smoker)
Drinking  <- as.character(cov.selected$Currently.Drink)
Menarche  <- as.numeric(as.character(cov.selected$Menarche))
Menopause <- as.character(cov.selected$Menstrual.Status)
Parity    <- as.numeric(as.character(cov.selected$Number.of.Live.Births))
VD        <- as.character(cov.selected$Vitaminuse)
Age.FB    <- ifelse(cov.selected$Age.at.First.Birth < 25, 1,
                    ifelse(cov.selected$Age.at.First.Birth < 30, 2,
                           ifelse(cov.selected$Age.at.First.Birth < 35, 3,4)
			                    )
	           )


########################################################################################
# MAKE COV FILE
########################################################################################
options(na.action='na.pass')
df.selected <- data.frame(Age, Race, Batch)
X           <- model.matrix(~Age + Race + Batch, df.selected)
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
B.seq.all <- B.seq.all[rowSums(is.na(B.seq.all)) < 10,]
#B.seq.all <- B.seq.all[complete.cases(B.seq.all), ]

write.table(B.seq.all, file = METHFILE, quote = F, append = FALSE, sep = "\t", row.names=FALSE)
rm(B.seq.all)
