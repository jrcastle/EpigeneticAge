setwd('/home/jrca253/EpigeneticAge/data')

DATADIR  <- '/home/jrca253/DATA/Truseq/'
METHFILE <- 'meth.txt'
COVFILE  <- 'cov.txt'

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


ID <- as.vector(cov$ID)
cov$Race <- as.character(cov$Race)
cov$Race <- ifelse(cov$Race == 'AfricanAmerican', 'African American', cov$Race)
cov$dataframe_name <- paste(cov$Batch,cov$ID, sep = '_')


########################################################################################
# SELECT ONLY NORMAL SAMPLES
########################################################################################
cov_norm <- cov[which(substr(cov$Sample.id,1,1)=='K'),]
cov_norm[cov_norm==''] <- NA


########################################################################################
# OMIT MISSING SAMPLES
########################################################################################
cov_norm <- cov_norm[which(cov_norm$ID<65 |cov_norm$ID>72),]
ID <- as.vector(cov_norm$ID)

dataframe_name = as.character(cov_norm[which(cov_norm$ID<=96),]$Sample.id)
dataframe_name = do.call(c, list(dataframe_name, as.character(cov_norm[which(cov_norm$ID>96),]$dataframe_name)))


########################################################################################
# COVARIATES
########################################################################################
Age       <- as.numeric(as.character(cov_norm$Age))
Race      <- as.character(cov_norm$Race)
Batch     <- as.character(cov_norm$Batch)
Sampleid  <- as.character(cov_norm$Sample.id)
                         
BMI       <- as.numeric(as.character(cov_norm$BMI))
Smoking   <- as.character(cov_norm$Current.Smoker)
Drinking  <- as.character(cov_norm$Currently.Drink)
Menarche  <- as.numeric(as.character(cov_norm$Menarche))
Menopause <- as.character(cov_norm$Menstrual.Status)
Parity    <- as.numeric(as.character(cov_norm$Number.of.Live.Births))
VD        <- as.character(cov_norm$Vitaminuse)
Age_FB    <- ifelse(cov_norm$Age.at.First.Birth < 25, 1,
                    ifelse(cov_norm$Age.at.First.Birth < 30, 2,
                           ifelse(cov_norm$Age.at.First.Birth < 35, 3,4)
			  )
	           )


########################################################################################
# MAKE COV FILE
########################################################################################
dfCov <- data.frame(Age, Batch)

X           <- model.matrix(~Age + Race + Batch, dfCov)
X           <- t(X)
colnames(X) <- dataframe_name
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

if(BATCH == 'HS_032'|| BATCH == 'HS_049' ){
  SAMPLEID <- Sampleid[i]
}

FILENAME <- paste(DIR, SAMPLEID, '_inTruseq.csv', sep = '')

B_seq_all           <- read.csv(FILENAME, header = T)
B_seq_all           <- B_seq_all[ which(B_seq_all$N >= 10),]
B_seq_all$position  <- paste(B_seq_all$chr, B_seq_all$pos, sep = ':')
B_seq_all$SAMPLEID  <- B_seq_all$X / B_seq_all$N
colnames(B_seq_all) <- c('chr','pos','N','X','position', SAMPLEID)
B_seq_all           <- B_seq_all[ c('position', SAMPLEID) ]

##### PROCESS AND MERGE REMAINING SAMPLES #####
for( i in 2:length(ID) ){
  print( paste("Processing sample", i, 'of', length(ID), sep = ' ') )

  BATCH    <- Batch[i]
  DIR      <- paste(DATADIR, BATCH, "/", sep = '')
  SAMPLEID <- ID[i]

  if(BATCH == 'HS_032'|| BATCH == 'HS_049' ){
    SAMPLEID <- Sampleid[i]
  }
  
  FILENAME <- paste(DIR, SAMPLEID, '_inTruseq.csv', sep = '')
  
  tmp           <- read.csv(FILENAME, header = T)
  tmp           <- tmp[ which(tmp$N >= 10),]
  tmp$position  <- paste(tmp$chr, tmp$pos, sep = ':')
  tmp$SAMPLEID  <- tmp$X / tmp$N
  colnames(tmp) <- c('chr','pos','N','X','position', SAMPLEID)
  tmp           <- tmp[ c('position', SAMPLEID) ]
  
  B_seq_all <- Reduce(
    function(x, y) merge(x, y, by = c('position'), all=T),
    list(B_seq_all, tmp)
  )
  rm(tmp)
        
}

write.table(B_seq_all, file = METHFILE, quote = F, append = FALSE, sep = "\t", row.names=FALSE)
rm(B_seq_all)
