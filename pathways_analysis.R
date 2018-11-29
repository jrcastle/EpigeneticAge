rm(list=ls()); gc();
setwd("/Users/jrca253/Documents/EpigeneticAge/test_code")

seed      <- "123"
model.dir <- paste("cpgs_in_KNT_imputed_seed", seed, "/", sep = '')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("seq2pathway", version = "3.8")
#BiocManager::install("AnnotationHub", version = "3.8")
#BiocManager::install("Homo.sapiens", version = "3.8")
#BiocManager::install("Organism.dplyr", version = "3.8")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", version = "3.8")
#BiocManager::install("biomaRt", version = "3.8")
#BiocManager::install("TxDb.Athaliana.BioMart.plantsmart22", version = "3.8")

library(ggplot2)
library(seq2pathway)
library(AnnotationHub)
library(Homo.sapiens)
library(Organism.dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(TxDb.Athaliana.BioMart.plantsmart22)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(Homo.sapiens)
library(goseq)


wrap.it <- function(x, len){ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}

###########################################################################################
# GENE FUNCTIONS
###########################################################################################
geneRanges <- function(db, column="SYMBOL")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

splitByOverlap <- function(query, subject, column="ENTREZID", ...)
{
  olaps <- findOverlaps(query, subject, maxgap = 0)
  f1 <- factor(subjectHits(olaps),
               levels=seq_len(subjectLength(olaps)))
  splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}

gns = geneRanges(Homo.sapiens,column = 'SYMBOL')

geneRanges <- function(db, column="ENTREZID")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

gns.Ensembl = geneRanges(Homo.sapiens, column = 'ENTREZID')
promoters_txdb <- promoters(gns)
promoters_txdb_ENTREZID <- promoters(gns.Ensembl)

##################################################################################
# LOAD MODEL COEFFICIENTS
##################################################################################
df.model <- read.csv( paste(model.dir, "model_coefficients.csv", sep = '') )
df.model <- df.model[-c(1),]
df.model$ID <- df.model$model.coefficients.name
df.model <- df.model %>%
  separate(model.coefficients.name, c("chr", "start"), ":")
df.model$start <- as.numeric(df.model$start)
df.model$end <- df.model$start #+ 1
df.model$model.coefficients.x <- NULL
df.model <- df.model[, c('ID', 'chr', 'start', 'end')]

##################################################################################
# PATHWAYS
##################################################################################


b <- runseq2gene(as(df.model,'GRanges'),
                 search_radius=150000, promoter_radius=200, promoter_radius2=100,
                 genome=c("hg19"), adjacent=FALSE, SNP=FALSE,
                 PromoterStop=FALSE,NearestTwoDirection=TRUE,UTR3=FALSE)

c <- runseq2pathway(as(df.model,'GRanges'),
                    search_radius=150000, promoter_radius=200, promoter_radius2=100,
                    genome=c('hg19'), adjacent=FALSE, SNP= FALSE,
                    PromoterStop=FALSE, NearestTwoDirection=TRUE,UTR3=FALSE,
                    DataBase=c("GOterm"), FAIMETest=FALSE, FisherTest=TRUE,
                    collapsemethod=c("Average"),
                    alpha=5, logCheck=FALSE, B=100, na.rm=FALSE, min_Intersect_Count=5)

d <- c$gene2pathway_result.FET$GO_BP
write.csv(d, paste(model.dir, 'pathway_BP.csv', sep = ''), row.names = FALSE)

d <- c$gene2pathway_result.FET$GO_CC
write.csv(d, paste(model.dir, 'pathway_CC.csv', sep = ''), row.names = FALSE)

d <- c$gene2pathway_result.FET$GO_MF
write.csv(d, paste(model.dir, 'pathway_MF.csv', sep = ''), row.names = FALSE)

##### PLOT #####
df.top.hits <- read.csv( paste(model.dir, "pathways_merged_tophits.csv", sep = "") )

df.top.hits$mlog10p <- -log10(df.top.hits$Fisher_Pvalue)
df.top.hits <- data.frame(df.top.hits$Description, df.top.hits$mlog10p, df.top.hits$Type)
colnames(df.top.hits) <- c("Description", "mlog10p", "Type")
df.top.hits <- df.top.hits[order(df.top.hits$mlog10p),]
df.top.hits$Description <- wrap.it(df.top.hits$Description, 80)
df.top.hits$Description <- factor(df.top.hits$Description, levels = unique(as.character(df.top.hits$Description)))

p <- ggplot(data = df.top.hits, aes(x = Description, y = mlog10p, fill = Type) ) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  scale_color_manual(
    name = "Type",
    values = c("red", "blue", "green"), 
    labels = c("BP", "CC", "MF")
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(0, 16)
  ) + 
  labs(x = "") + 
  labs(y = "-log10(p)") + 
  theme_bw(base_size = 15) + 
  theme(
    legend.position=c(0.8, 0.2),
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    axis.text.y = element_text(size = 8.5, face = "bold", margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
    #plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points")
  ); p

png( paste(model.dir, "pathways_result_tophits.png", sep = ''), width = 650, height = 500, units = 'px')
p
dev.off()

##################################################################################
# GENE LIST
##################################################################################

df.model$promoter <- paste((splitByOverlap(promoters_txdb, as(df.model,'GRanges'), "SYMBOL")),collapse="/")
df.model$promoter_ENTREZID <- paste((splitByOverlap(promoters_txdb_ENTREZID, as(df.model,'GRanges'), "ENTREZID")),collapse="/")
df.model$gene <- paste((splitByOverlap(gns, as(df.model,'GRanges'), "SYMBOL")),collapse="/")
df.model$gene_ENTREZID <- paste((splitByOverlap(gns.Ensembl, as(df.model,'GRanges'), "ENTREZID")),collapse="/")

t<-nearest(as(df.model,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
df.model$nearest <- nearest$SYMBOL
df.model$nearest <- ifelse(df.model$gene == '' & df.model$promoter == '', df.model$nearest, '')

t<-nearest(as(df.model,'GRanges'),gns.Ensembl,ignore.strand=T)
nearest <- gns.Ensembl[t]
df.model$nearest_ENTREZID <- nearest$ENTREZID
df.model$nearest_ENTREZID <- ifelse(df.model$gene == '' & df.model$promoter == '', df.model$nearest, '')

df.model <- separate_rows(df.model, gene)
df.model <- separate_rows(df.model, promoter)
df.model$genelist <- ifelse(df.model$promoter == '', df.model$gene, df.model$promoter)
df.model$genelist <- ifelse(df.model$gene == '' & df.model$promoter == '', df.model$nearest, df.model$genelist)
genelist <- as.character(df.model[!duplicated(df.model$genelist),]$genelist)
gene_description <- addDescription(genome="hg19",genevector=genelist)
#df.model$description <- gene_description

write.csv(df.model,paste(model.dir, 'genelist.csv', sep = ''), row.names = FALSE)
head(gene_description)
