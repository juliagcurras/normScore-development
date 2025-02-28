
# Julia G Currás - 2025/02/20

# Just another script to test the normScore 



# Set up ####
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")
library(dplyr)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")




# Example data ####

## Load quantification matrix ####
archivo <- "../Data/trialDatashort.xlsx" 
intensityMatrix <- as.data.frame(readxl::read_excel(path = archivo, sheet = 1,
                                                    col_names = TRUE))
intensityMatrix[,-1] <- apply(intensityMatrix[,-1], 2, as.numeric)
rownames(intensityMatrix) <- intensityMatrix$`Sample Name`
intensityMatrix <- intensityMatrix[,-1]
logIntensityMatrix <- as.matrix(log(intensityMatrix, base = 2))
logIntensityMatrix[is.infinite(logIntensityMatrix)] <- NA

## Load design matrix ####
archivoG <- "../Data/LabelsShort.xlsx"
dfGrupos <- as.data.frame(readxl::read_excel(path = archivoG, sheet = 1, 
                                             col_names = T))[-3]
colnames(dfGrupos) <- c("Samples", "Groups")
grupos <- levels(as.factor(dfGrupos$Groups))

## Normalization ####
mydata <- list(log = logIntensityMatrix,
               Mean = meanNorm(rawMatrix = intensityMatrix), 
               Median = medianNorm(rawMatrix = intensityMatrix), 
               GI = GINorm(rawMatrix = intensityMatrix),
               Quantile = quantileNorm(log2Matrix = logIntensityMatrix),
               VSN = VSNNorm(rawMatrix = intensityMatrix),
               CyclicLoess = cyclicLoessNorm(log2Matrix = logIntensityMatrix),
               RLR = RLRNorm(log2Matrix = logIntensityMatrix),
               MAD = MADNormalization(log2Matrix = logIntensityMatrix))

## Assessment: normScore ####
resultado <- normScore(normMatrixList = mydata, 
                       designMatrix = dfGrupos, 
                       dfRaw = intensityMatrix)

resultado$detailRaking
resultado$detailScore
resultado$finalRanking



# Benchmarking data DIA ####
ruta <- "H:/Mi unidad/Doctorado/SEProt_ProteoRed/Data/NewData/DIA/"
designMatrix <- readRDS(file = paste0(ruta, "3_desingMatrix.rds"))
datosNorm <- readRDS(file = paste0(ruta, "3_normData.rds"))
rawMatrix <- readRDS(file = paste0(ruta, "3_rawData.rds"))
refGroup <- "A"
altGroup <- "B"

## One dataset ####
resultado <- normScore(normMatrixList = datosNorm$`88min_200ng`,
                       designMatrix = designMatrix,
                       dfRaw = rawMatrix$`88min_200ng`,
                       refGroup = refGroup, 
                       altGroup = altGroup)

resultado$finalRanking

## All datasets ####
totalResults <- sapply(names(datosNorm), function (i) normScore( 
  designMatrix = designMatrix, 
  normMatrixList = datosNorm[[i]], 
  dfRaw = rawMatrix[[i]],
  refGroup = refGroup, 
  altGroup = altGroup), simplify = F)

scoreAllDatasets <- sapply(totalResults, "[[", 1, simplify = F)
scoreAllDatasets
View(scoreAllDatasets)

totalResults$`44min_200ng`$detailRaking %>% View
totalResults$`44min_200ng`$detailScore %>% View

rankingAllDatasets <- sapply(scoreAllDatasets, function(i)  names(i))
View(rankingAllDatasets)




# Benchmarking data DDA ####
ruta <- "H:/Mi unidad/Doctorado/SEProt_ProteoRed/Data/NewData/DDA/"
designMatrix <- readRDS(file = paste0(ruta, "3_desingMatrix.rds"))
datosNorm <- readRDS(file = paste0(ruta, "3_normData.rds"))
rawMatrix <- readRDS(file = paste0(ruta, "3_rawData.rds"))
refGroup <- "A"
altGroup <- "B"

## One dataset ####
resultado <- normScore(normMatrixList = datosNorm$`88min 200ng`,
                       designMatrix = designMatrix,
                       dfRaw = rawMatrix$`88min 200ng`,
                       refGroup = refGroup, 
                       altGroup = altGroup)
resultado$finalRanking

## All datasets ####
totalResults <- sapply(names(datosNorm), function (i) normScore( 
  designMatrix = designMatrix, 
  normMatrixList = datosNorm[[i]], 
  dfRaw = rawMatrix[[i]],
  refGroup = refGroup, 
  altGroup = altGroup), simplify = F)

scoreAllDatasets <- sapply(totalResults, "[[", 1, simplify = F)
scoreAllDatasets
View(scoreAllDatasets)

rankingAllDataseta <- sapply(scoreAllDatasets, function(i)  names(i))
View(rankingAllDataseta)

totalResults$`44min 200ng`$detailRaking %>% View
totalResults$`44min 200ng`$detailScore %>% View





















