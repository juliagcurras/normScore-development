##############################################################################-

#############        SCORE TO ASSESS NORMALIZATIONS      #####################-

##############################################################################-


# Julia G Currás - 21/01/2025

# Set up ####
setwd("C:/Users/julia/Documents/GitHub/normScore/R")
library(dplyr)
source(file = "supportFunctions.R", encoding = "UTF-8")



# Trial data ####
## Load quantification matrix 
archivo <- "../Datos/trialDatashort.xlsx" ####
intensityMatrix <- as.data.frame(readxl::read_excel(path = archivo, sheet = 1,
                                                    col_names = TRUE))
intensityMatrix[,-1] <- apply(intensityMatrix[,-1], 2, as.numeric)
rownames(intensityMatrix) <- intensityMatrix$`Sample Name`
intensityMatrix <- intensityMatrix[,-1]
logIntensityMatrix <- as.matrix(log(intensityMatrix, base = 2))
logIntensityMatrix[is.infinite(logIntensityMatrix)] <- NA

## Load design matrix ####
archivoG <- "../Datos/LabelsShort.xlsx"
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



# This score will be divided in items - one item by each metric/graphic of assessment #

scoreFinal <- list()

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


# ITEM 5 - RLE: MSE median sample (ref = 0) ####
item5 <- sapply(mydata, rleMSE)
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





# Función ####

normScore <- function(normMatrixList, designMatrix, 
                      refGroup = NULL, altGroup = NULL){
  # Input: 
  # 1. List of normalized matrix (normMatrixList)
  # 2. Design matrix (designMatrix)
  # 3. Ref group and alternative group (optional). If they are not provided, 
  # the first group will be used as alternative and the last one, as control.
  
  totalGroups <- levels(as.factor(designMatrix$Groups))
  scoreFinal <- list()
  
  # ITEM 1 - PVC ####
  dfPCV <- data.frame(lapply(normMatrixList, getPCV, grupos = totalGroups, 
                             dfGrupos = designMatrix))
  dfPCV <- as.data.frame(t(dfPCV))
  dfPCV$PCV <- apply(dfPCV, 1, mean, na.rm = T)
  item1 <- stats::setNames(dfPCV$PCV, rownames(dfPCV))
  
  scoreFinal[["PCV"]] <- item1
  
  
  # ITEM 2 - Correlation (Spearman) ####
  allVectorsCorr <- lapply(normMatrixList, getCorrelationVector,
                           dfGrupos = designMatrix,
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
  refGroup <- ifelse(is.null(refGroup), totalGroups[length(totalGroups)], refGroup)
  altGroup <- ifelse(is.null(altGroup), totalGroups[1], altGroup)
  samplesG1 <- designMatrix[designMatrix$Groups == refGroup, "Samples"]
  samplesG2 <- designMatrix[designMatrix$Groups == altGroup, "Samples"]
  
  item3 <- sapply(normMatrixList, maDiffAreas, samplesG1 = samplesG1, 
                  samplesG2 = samplesG2)
  scoreFinal[["MAplot"]] <- item3
  
  
  # ITEM 4 - MeanSD  regression line B = 0 ####
  item4 <- sapply(normMatrixList, meanSDdiffArea)
  scoreFinal[["MeanSDplot"]] <- item4
  
  
  # ITEM 5 - RLE: MSE median sample (ref = 0) ####
  item5 <- sapply(normMatrixList, rleMSE)
  scoreFinal[["RLEplot"]] <- item5
  
  
  # ITEM 6 - total intensity: MSE median sample (ref = global median) ####
  item6 <- sapply(normMatrixList, tiMSE)
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
  finalRank <- stats::setNames(rankingDF$Total, rownames(rankingDF))
  
  
  return(list(finalRanking = finalRank, 
              detailRaking = rankingDF, 
              detailScore = scoreDF))
}

resultado <- normScore(normMatrixList = mydata, designMatrix = dfGrupos)
resultado$detailScore
resultado$detailRaking
resultado$finalRanking


# Trying other datasets ####
ruta <- "../../SEProt_ProteoRed/Data/NewData/DIA/"
designMatrix <- readRDS(file = paste0(ruta, "3_desingMatrix.rds"))
datosNorm <- readRDS(file = paste0(ruta, "3_normData.rds"))
refGroup <- "A"
altGroup <- "B"

resultado <- normScore(normMatrixList = datosNorm$`88min_200ng`,
                       designMatrix = designMatrix, 
                       refGroup = refGroup, 
                       altGroup = altGroup)

totalResults <- lapply(datosNorm, normScore, 
       designMatrix = designMatrix, 
       refGroup = refGroup, 
       altGroup = altGroup)

totalResults$`44min_200ng`$detailRaking %>% View
totalResults$`44min_200ng`$detailScore %>% View
scoreAllDatasets <- sapply(totalResults, "[[", 1, simplify = F)
scoreAllDatasets
View(scoreAllDatasets)
rankingAllDataseta <- sapply(scoreAllDatasets, function(i)  names(i))
View(rankingAllDataseta)






















