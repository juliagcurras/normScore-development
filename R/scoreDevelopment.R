##############################################################################-

#############        SCORE TO ASSESS NORMALIZATIONS      #####################-

##############################################################################-


# Julia G Currás - 21/01/2025

# Set up ####
setwd("C:/Users/julia/Documents/GitHub/normScore/R")
library(dplyr)
source(file = "supportFunctions.R", encoding = "UTF-8")


#...........................................................................####
# Trial data ####
## Load quantification matrix 
archivo <- "../Data/trialDatashort.xlsx" ####
archivo <- "../Data/SEProt_11min200ng.xlsx" ####
intensityMatrix <- as.data.frame(readxl::read_excel(path = archivo, sheet = 1,
                                                    col_names = TRUE))
intensityMatrix[,-1] <- apply(intensityMatrix[,-1], 2, as.numeric)
rownames(intensityMatrix) <- intensityMatrix$`Sample Name`
intensityMatrix <- intensityMatrix[,-1]
logIntensityMatrix <- as.matrix(log(intensityMatrix, base = 2))
logIntensityMatrix[is.infinite(logIntensityMatrix)] <- NA

## Load design matrix ####
archivoG <- "../Data/LabelsShort.xlsx"
archivoG <- "../Data/SEProt_11min200ng_groups.xlsx"
dfGrupos <- as.data.frame(readxl::read_excel(path = archivoG, sheet = 1, 
                                               col_names = T))[-3]
colnames(dfGrupos) <- c("Samples", "Groups")
grupos <- levels(as.factor(dfGrupos$Groups))

## Normalization ####
mydata <- list(Log = logIntensityMatrix,
               Mean = meanNorm(rawMatrix = intensityMatrix), 
               Median = medianNorm(rawMatrix = intensityMatrix), 
               GI = GINorm(rawMatrix = intensityMatrix),
               Quantile = quantileNorm(log2Matrix = logIntensityMatrix),
               VSN = VSNNorm(rawMatrix = intensityMatrix),
               CyclicLoess = cyclicLoessNorm(log2Matrix = logIntensityMatrix),
               RLR = RLRNorm(log2Matrix = logIntensityMatrix),
               MAD = MADNormalization(log2Matrix = logIntensityMatrix))



#...........................................................................####
# This score will be divided in items - one item by each metric/graphic of assessment #

scoreFinal <- list()

# ITEM 0 - raw intensity ####

total_intensities <- colSums(intensityMatrix, na.rm = T)
cv(total_intensities, proportion = T, na.rm = T)
cv <- sd(total_intensities) / mean(total_intensities)
cv


# ITEM 1 - PVC ####
dfPCV <- data.frame(lapply(mydata, getPCV, grupos = grupos, dfGrupos = dfGrupos))
dfPCV <- as.data.frame(t(dfPCV))
dfPCV$PCV <- apply(dfPCV, 1, mean, na.rm = T)
item1 <- stats::setNames(dfPCV$PCV, rownames(dfPCV))

scoreFinal[["PCV"]] <- item1


# ITEM 2 - Correlation (Spearman) ####
allVectorsCorr <- lapply(mydata, getCorrelationVector,
                         dfGrupos = dfGrupos,
                         metodo = "spearman")

dfCor <- data.frame(sapply(allVectorsCorr, "length<-", 
                           max(lengths(allVectorsCorr))))

item2 <- sapply(colnames(dfCor), function(j) {
  i <- dfCor[,j]
  # faigo  1-correlation porque todas as métricas restantes siguen o patrón a 
  # menor número, mellor é a métrica, para que está tamén sexa así. 
  1-(median(i, na.rm = T)-IQR(i, na.rm = T)/3) 
}, simplify = T, USE.NAMES = T)

scoreFinal[["Correlation"]] <- item2


# ITEM 3 - MAplot regression line 0 ####
refGroup <- (grupos)[length(grupos)]
altGroup <- (grupos)[1]
samplesG1 <- dfGrupos[dfGrupos$Groups == refGroup, "Samples"]
samplesG2 <- dfGrupos[dfGrupos$Groups == altGroup, "Samples"]

item3 <- sapply(mydata, maDiffAreas, samplesG1 = samplesG1, 
                samplesG2 = samplesG2)
scoreFinal[["MAplot"]] <- item3


# ITEM 4 - MeanSD  regression line B = 0 ####
item4 <- sapply(mydata, meanSDdiffArea)
scoreFinal[["MeanSDplot"]] <- item4


# ITEM 5 - RLE:  ####

# adjustment to distribution N(0,1) but SD not equal to 1 
rleKS <- function(dfDatos) {
  medianaProt <- apply(dfDatos, 1, stats::median, na.rm = T)
  
  rleData <- as.data.frame(log(t(t(dfDatos) / medianaProt), base = 2))
  
  ksStatistic <- apply(rleData, 2, function(i) ks.test(x = i, y = "pnorm", mean = 0, sd = 1)$statistic)
  
  return(median(ksStatistic))
}
item5 <- sapply(mydata, rleMAPE)


# MAPE removing the logarithm
mape <- function(actual, predicted, prop = F){
  metric <- mean(abs((actual - predicted)/actual))
  metric <- ifelse(!prop, metric*100, metric)
  return(metric)
}

rleMAPE <- function(dfDatos) {
  # dfDatos <- mydata$Mean
  medianaProt <- apply(dfDatos, 1, stats::median, na.rm = T)
  
  ## log data
  # rleData <- as.data.frame(log(t(t(dfDatos) / medianaProt), base = 2))
  
  ## non log data
  rleData <- as.data.frame(t(t(dfDatos) / medianaProt))
  
  medianVector <- apply(rleData, 2, median, na.rm = T)
  
  return(mape(actual = 1, predicted = medianVector, prop = T))
}

item5 <- sapply(mydata, rleMAPE)
scoreFinal[["RLEplot"]] <- item5


# ITEM 6 - total intensity: MSE median sample (ref = 0) ####
item6 <- sapply(mydata, tiMSE)
scoreFinal[["totalIntensity"]] <- item6


# Join everything ####

## Bind ####
scoreDF <- dplyr::bind_cols(scoreFinal)
scoreDF <- as.data.frame(scoreDF)
rownames(scoreDF) <- names(scoreFinal[[1]])

## Rank ####
rankingDF <- as.data.frame(apply(scoreDF, 2, dplyr::dense_rank, simplify = T))
rownames(rankingDF) <- rownames(scoreDF)
rankingDF$Total <- rowSums(rankingDF)

## Sort ####
rankingDF <- rankingDF %>% dplyr::arrange(Total)
rankingDF





#...........................................................................####
# Function implementation in scoreFunction.R ####
source(file = "scoreFunction.R", encoding = "UTF-8")

## Trying function ####

resultado <- normScore(normMatrixList = mydata, 
                       designMatrix = dfGrupos, 
                       dfRaw = intensityMatrix)
resultado$detailScore
resultado$detailRaking
resultado$finalRanking


























