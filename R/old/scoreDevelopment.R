##############################################################################-

#############        SCORE TO ASSESS NORMALIZATIONS      #####################-

##############################################################################-


# Julia G Currás - 21/01/2025

# Set up ####
setwd("C:/Users/julia/Documents/GitHub/normScore/R")
library(dplyr)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")


#...........................................................................####
# Trial data ####
## Load quantification matrix 
archivo <- "../Data/trialDatashort.xlsx" ####
# archivo <- "../Data/SEProt_11min200ng.xlsx" ####
intensityMatrix <- as.data.frame(readxl::read_excel(path = archivo, sheet = 1,
                                                    col_names = TRUE))
intensityMatrix[,-1] <- apply(intensityMatrix[,-1], 2, as.numeric)
rownames(intensityMatrix) <- intensityMatrix$`Sample Name`
intensityMatrix <- intensityMatrix[,-1]
logIntensityMatrix <- as.matrix(log(intensityMatrix, base = 2))
logIntensityMatrix[is.infinite(logIntensityMatrix)] <- NA

## Load design matrix ####
archivoG <- "../Data/LabelsShort.xlsx"
# archivoG <- "../Data/SEProt_11min200ng_groups.xlsx"
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
item0 <- sd(total_intensities) / mean(total_intensities)
item0


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
item6 <- sapply(mydata, tiMAPE)
scoreFinal[["totalIntensity"]] <- item6


# Join everything ####

## Bind ####
scoreDF <- dplyr::bind_cols(scoreFinal)
scoreDF <- as.data.frame(scoreDF)
rownames(scoreDF) <- names(scoreFinal[[1]])

## Scale ####
# media <- mean(scoreDF$PCV)
# desEst <- sd(scoreDF$PCV)
# (scoreDF$PCV - media)/desEst
# (scoreDF$PCV - min(scoreDF$PCV)) / (max(scoreDF$PCV) - min(scoreDF$PCV))

scoreDF_norm <- apply(scoreDF, 2, function(col) (col - min(col)) / (max(col) - min(col)))

scoreDF_norm <- as.data.frame(scoreDF_norm)


## Rank ####
rownames(scoreDF_norm) <- rownames(scoreDF)
scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), ] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), ]*item0
scores_matrix <- t(scoreDF_norm)
scoreDF_norm$Total <- rowSums(scoreDF_norm)
# scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "Total"] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "Total"]*item0

scoreDF_norm

## Sort ####
scoreDF_norm <- scoreDF_norm %>% dplyr::arrange(Total)
scoreDF_norm


# CI bootstrap ####

library(boot)

# scores_matrix <- t(scoreDF_norm)

# Computing total score for each normalization after resampling proteins (rows)
bootstrap_score_rows <- function(data, indices) {
  resampled_matrix <- data[indices, , drop = FALSE]
  total_scores <- colSums(resampled_matrix)
  return(total_scores)  # One score per normalization
}

# Bootstrap
n_boot <- 1000
boot_results <- boot(data = scores_matrix,              # data
                     statistic = bootstrap_score_rows,  # function for getting the scores by nomralization
                     R = n_boot)                        # number of resamples  

# Output mean scores and confidence intervals
bootstrap_means <- colMeans(boot_results$t)

score_bootstrap <- as.data.frame(t(sapply(1:ncol(scores_matrix), function(i){
  ci <- boot.ci(boot_results, type = "perc", index = i)
  return(c(colnames(scores_matrix)[i], bootstrap_means[i], 
           ci$percent[4], ci$percent[5]))
}, simplify = T)))

colnames(score_bootstrap) <- c("Normalization", "Mean Total Score", "LL95%", "UL95%")





#...........................................................................####
# Function implementation in scoreFunction.R ####
source(file = "scoreFunction.R", encoding = "UTF-8")

## Trying function ####

resultado <- normScore(normMatrixList = mydata, 
                       designMatrix = dfGrupos, 
                       dfRaw = intensityMatrix)
resultado$detailScore
resultado$detailRanking
resultado$finalRanking
resultado$bootstrapScore
resultado$graphic


























