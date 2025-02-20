#### SUPPORT FUNCTIONS ####

# setwd("~/GitHub/Normalization/app")
# library(shiny)
library(dplyr)
library(tidyr)
library(stats)
library(MASS)

options(repos = BiocManager::repositories())
# library(MSnSet.utils)
library(vsn)
library(preprocessCore)


# Normalization methods ####

meanNorm <- function(rawMatrix){
  
  # Do calculations
  colMeans <- colMeans(rawMatrix, na.rm = TRUE)
  avgColMean <- mean(colMeans, na.rm = TRUE)
  
  # Create empty matrix
  normMatrix <- matrix(nrow = nrow(rawMatrix), ncol = ncol(rawMatrix), 
                       byrow = TRUE)
  # Normalization 
  normFunc <- function(colIndex) {
    (rawMatrix[rowIndex, colIndex]/colMeans[colIndex]) * 
      avgColMean
  }
  for (rowIndex in seq_len(nrow(rawMatrix))) {
    normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), 
                                     normFunc, 0)
  }
  # Log transformation
  normLog2Matrix <- log2(normMatrix)
  colnames(normLog2Matrix) <- colnames(rawMatrix)
  
  # Remove is.infinite data
  normLog2Matrix[sapply(normLog2Matrix, is.infinite)] <- NA
  
  # Output
  return(normLog2Matrix)
}
# meanData <- meanNorm(rawMatrix = intensityMatrix)

medianNorm <- function(rawMatrix) {
  
  rawMatrix <- as.matrix(rawMatrix)
  
  # Do calculations
  colMedians <- matrixStats::colMedians(rawMatrix, na.rm = TRUE)
  meanColMedian <- mean(colMedians, na.rm = TRUE)
  
  # Create empty matrix
  normMatrix <- matrix(nrow = nrow(rawMatrix), ncol = ncol(rawMatrix), 
                       byrow = TRUE)
  
  # Normalization 
  normFunc <- function(colIndex) {
    (rawMatrix[rowIndex, colIndex]/colMedians[colIndex]) * 
      meanColMedian
  }
  for (rowIndex in seq_len(nrow(rawMatrix))) {
    normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), 
                                     normFunc, 0)
  }
  # Log transformation
  normLog2Matrix <- log2(normMatrix)
  colnames(normLog2Matrix) <- colnames(rawMatrix)
  
  # Remove is.infinite data
  normLog2Matrix[sapply(normLog2Matrix, is.infinite)] <- NA
  
  # Output
  normLog2Matrix
}

GINorm <- function(rawMatrix){
  
  # Do calculations
  colSums <- colSums(rawMatrix, na.rm = TRUE)
  colSumsMedian <- stats::median(colSums)
  
  # Create empty matrix
  normMatrix <- matrix(nrow = nrow(rawMatrix), ncol = ncol(rawMatrix), 
                       byrow = TRUE)
  
  # Normalization 
  normFunc <- function(colIndex) {
    (rawMatrix[rowIndex, colIndex]/colSums[colIndex]) * 
      colSumsMedian
  }
  for (rowIndex in seq_len(nrow(rawMatrix))) {
    normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), 
                                     normFunc, 0)
  }
  
  # Log transformation
  normLog2Matrix <- log2(normMatrix)
  colnames(normLog2Matrix) <- colnames(rawMatrix)
  
  # Remove is.infinite data
  normLog2Matrix[sapply(normLog2Matrix, is.infinite)] <- NA
  
  # Output
  normLog2Matrix
}

VSNNorm <- function (rawMatrix) {
  rawMatrix <- as.matrix(rawMatrix)
  
  normMatrix <- suppressMessages(vsn::justvsn(rawMatrix))
  colnames(normMatrix) <- colnames(rawMatrix)
  # Remove is.infinite data
  normMatrix[sapply(normMatrix, is.infinite)] <- NA
  normMatrix
}

quantileNorm <- function(log2Matrix){
  log2Matrix <- as.matrix(log2Matrix)
  normMatrix <- preprocessCore::normalize.quantiles(log2Matrix, 
                                                    copy = TRUE)
  colnames(normMatrix) <- colnames(log2Matrix)
  # Remove is.infinite data
  normMatrix[sapply(normMatrix, is.infinite)] <- NA
  normMatrix
}

cyclicLoessNorm <- function (log2Matrix) {
  log2Matrix <- as.matrix(log2Matrix)
  normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method = "fast")
  colnames(normMatrix) <- colnames(log2Matrix)
  # Remove is.infinite data
  normMatrix[sapply(normMatrix, is.infinite)] <- NA
  normMatrix
}

RLRNorm <- function(log2Matrix){
  
  log2Matrix <- as.matrix(log2Matrix)
  
  log2Matrix[is.infinite(log2Matrix)] <- NA
  
  sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm = TRUE)
  
  calculateRLMForCol <- function(colIndex, sampleLog2Median, 
                                 log2Matrix) {
    lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex]) ~ 
                         sampleLog2Median, na.action = stats::na.exclude)
    coeffs <- lrFit$coefficients
    coefIntercept <- coeffs[1]
    coefSlope <- coeffs[2]
    globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept)/coefSlope
    globalFittedRLRCol
  }
  globalFittedRLR <- vapply(seq_len(ncol(log2Matrix)), calculateRLMForCol, 
                            rep(0, nrow(log2Matrix)), 
                            sampleLog2Median = sampleLog2Median, 
                            log2Matrix = log2Matrix)
  colnames(globalFittedRLR) <- colnames(log2Matrix)
  
  # Remove is.infinite data
  globalFittedRLR[sapply(globalFittedRLR, is.infinite)] <- NA
  
  globalFittedRLR
}

MADNormalization <- function(log2Matrix){
  
  log2Matrix <- as.matrix(log2Matrix)
  
  sampleLog2Median <- matrixStats::colMedians(log2Matrix, 
                                              na.rm = TRUE)
  sampleMAD <- matrixStats::colMads(log2Matrix, na.rm = TRUE)
  madMatrix <- t(apply(log2Matrix, 1, function(row) ((row - 
                                                        sampleLog2Median)/sampleMAD)))
  madPlusMedianMatrix <- madMatrix + mean(sampleLog2Median)
  colnames(madPlusMedianMatrix) <- colnames(log2Matrix)
  # Remove is.infinite data
  madPlusMedianMatrix[sapply(madPlusMedianMatrix, is.infinite)] <- NA
  
  madPlusMedianMatrix
}


doNormalization <- function(listaNorm, rawData, logData){
  
  # rawData <- apply(rawData, 2, as.numeric)
  # rawData <- as.matrix(rawData)
  
  listaFinal <- list()
  if ("Mean" %in% listaNorm){
      # mediaData <- meanNorm(rawMatrix = rawData)
    listaFinal <- c(listaFinal, list(Mean = meanNorm(rawMatrix = rawData)))
  }
  if ("Median" %in% listaNorm){
    listaFinal <- c(listaFinal, list(Median = medianNorm(rawMatrix = rawData)))
  }
  if ("TI" %in% listaNorm){
    listaFinal <- c(listaFinal, list(TI = GINorm(rawMatrix = rawData)))
  } 
  if ("VSN" %in% listaNorm){
    listaFinal <- c(listaFinal, list(VSN = VSNNorm(rawMatrix = rawData)))
  }
  if ("Quantile" %in% listaNorm){
    listaFinal <- c(listaFinal, list(Quantile = quantileNorm(log2Matrix = logData)))
  }
  if ("CyclicLoess" %in% listaNorm){
    listaFinal <- c(listaFinal, list(CyclicLoess = cyclicLoessNorm(log2Matrix = logData)))
  }
  if ("RLR" %in% listaNorm){
    listaFinal <- c(listaFinal, list(RLR = RLRNorm(log2Matrix = logData)))
  }
  if ("MAD" %in% listaNorm){
    listaFinal <- c(listaFinal, list(MAD = MADNormalization(log2Matrix = logData)))
  }

  
  # Final output
  if (length(listaFinal) == 0){
    return(NULL)
  } else {
    return(listaFinal)
  }
  
}
# normData <- doNormalization(listaNorm = c("media", "mediana", "ti"),
#                             rawData = intensityMatrix, logData = NULL)
# 





# Support funcions for getting score #### 

diffAreas <- function(intPred, coefPred, minRange, maxRange, intExpected = 0){
  
  # Moving regression lines to reach B=0, a=0
  intPred <- intPred - intExpected
  
  if (coefPred < 0){
    a = 1
    b = -1
  } else if (coefPred > 0){
    a = -1
    b = 1
  }
  
  # Cutpoint regression line with expected line
  cpX <- (-intPred)/coefPred
  cpY <- 0
  
  # Is cutpoint located inside the range?
  if (all(cpX >= minRange, cpX <= maxRange)){
    area1 <- coefPred*a*(((minRange + (intPred/coefPred))^2)/2 - ((cpX + (intPred/coefPred))^2)/2)
    area2 <- coefPred*b*(((cpX + (intPred/coefPred))^2)/2 - ((maxRange + (intPred/coefPred))^2)/2)
    areaMetric <- abs(area1+area2)
  } else if (any(cpX < minRange, cpX > maxRange)){
    areaMetric <- abs(coefPred*(((maxRange + (intPred/coefPred))^2)/2 - ((minRange + (intPred/coefPred))^2)/2))
  }
  
  areaMetric <- areaMetric/(maxRange-minRange)
  return(areaMetric)
  
  # outra forma de calculalo
  # recta <- function(x){coefPred*x+intPred}
  # integrate(recta, minRange, maxRange)
}


mse <- function(actual, predicted){
  mean((actual - predicted)^2)
}

rleMSE <- function(dfDatos) {
  medianaProt <- apply(dfDatos, 1, stats::median, na.rm = T)
  
  rleData <- as.data.frame(log(t(t(dfDatos) / medianaProt), base = 2))
  
  medianVector <- apply(rleData, 2, median, na.rm = T)
  
  return(mse(actual = 0, predicted = medianVector))
}

tiMSE <- function(dfDatos) {
  allMetrics <- apply(dfDatos, 2, stats::quantile, na.rm = T, simplify = F)
  q1 <- sapply(allMetrics, "[[", 2)
  q3 <- sapply(allMetrics, "[[", 4)
  medianVector <- sapply(allMetrics, "[[", 3)
  finalMetric <- 
    mse(actual = median(medianVector), predicted = medianVector)/median(medianVector) + 
    mse(actual = median(q1), predicted = q1)/median(q1) +
    mse(actual = median(q3), predicted = q3)/median(q3)
  
  return(finalMetric)
}


meanSDdiffArea <- function(dfDatos){
  
  meanSamples <- apply(dfDatos, 2, mean, na.rm = T)
  sdSamples <- apply(dfDatos, 2, sd, na.rm = T)
  
  dfAux <- data.frame(Media = meanSamples,
                      DesvEst = sdSamples, 
                      Samples = colnames(dfDatos))
  
  dfAux <- dfAux %>% dplyr::arrange(Media)
  dfAux$Orden <- 1:nrow(dfAux)
  
  res <- lm(DesvEst~Orden, data = dfAux)
  coefPred <-  res$coefficients["Orden"]
  intPred <-  res$coefficients["(Intercept)"]
  
  resultArea <- diffAreas(
    intPred = 0, 
    coefPred = coefPred,
    maxRange = max(dfAux$Orden),
    minRange = min(dfAux$Orden),
    intExpected = 0
  )
  # names(resultArea) <- i
  return(unname(resultArea))
  
}


coefVariation <- function(x, na.rm = TRUE) {
  (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))*100
} 

cvGruposProt <- function(grupo, dfGrupos, dfDatos){
  
  dfGrupos <- as.data.frame(dfGrupos)
  dfDatos <- as.data.frame(dfDatos)
  espGroup <- as.vector(unlist(dfGrupos[dfGrupos$Groups == grupo, "Samples"]))

  df <- dfDatos[, as.character(espGroup)]
  apply(df, 1, coefVariation)
}

getPCV <- function(dfDatos, grupos, dfGrupos){
  # Calculate cv by groups for each protein
  dfGruposProt <- as.data.frame(sapply(grupos, cvGruposProt, 
                                       dfGrupos = dfGrupos, 
                                       dfDatos = dfDatos))
  # Calculate mean and mean IC of CV for each protein
  medias <- apply(dfGruposProt, 2, mean, na.rm = T)
  names(medias) <- colnames(dfGruposProt)
  return(medias)
}


maDiffAreas <- function(data, samplesG1, samplesG2){
  data <- as.data.frame(data)
  df <- data[, c(samplesG1, samplesG2)]
  df <- na.omit(df)
  
  df$logFC  <- apply(df, 1, function(x) {
    mean(x[samplesG2], na.rm =T) - mean(x[samplesG1], na.rm = T)
  }
  )
  
  df$AveExpr  <- apply(df, 1, function(x) {
    (mean(x[samplesG2], na.rm =T) + mean(x[samplesG1], na.rm = T))/2
  }
  )
  df <- df[, c("logFC", "AveExpr")]
  
  res <- lm(logFC~AveExpr, data = df)
  coefPred <-  res$coefficients["AveExpr"]
  intPred <-  res$coefficients["(Intercept)"]
  
  resultArea <- diffAreas(
    intPred = intPred, 
    coefPred = coefPred,
    maxRange = max(df$AveExpr),
    minRange = min(df$AveExpr),
    intExpected = 0
  )
  # names(resultArea) <- i
  return(unname(resultArea))
}



getCorrelationVector <- function(df, dfGrupos, metodo = "pearson"){
  df <- as.data.frame(df)
  allCorrs <- lapply(unique(dfGrupos$Groups), function(i){
    # Select samples
    samplesByGroup <- dfGrupos %>% 
      dplyr::filter(Groups == i) %>%
      dplyr::pull(Samples)
    dfCor <- df %>% 
      dplyr::select(all_of(samplesByGroup))
    # Estimate correlation
    tabCor <- stats::cor(dfCor, use = "complete.obs", method = metodo)
    # Get unique pairs of correlation
    tabCor[lower.tri(tabCor)] <- NA # remove duplicated corr
    diag(tabCor) <- NA # remove variance diagonal
    tabCor <- stats::na.omit(reshape2::melt(tabCor))
  }
  )
  vecFinal <- as.data.frame(dplyr::bind_rows(allCorrs)) %>% dplyr::pull(value)
  return(vecFinal)
}


# SCORE: main function #### 

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
  
  
  # ITEM 6 - total intensity: MSE median sample (ref = 0) ####
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


# resultado <- normScore(normMatrixList = mydata, designMatrix = dfGrupos)
# Example in scoreDevelopment





# Probas MSE ####
# vec1 <- rnorm(n = 100, mean = 1000, sd = 50)
# vec2 <- rnorm(n = 100, mean = 3, sd = 5)
# 
# mse(actual = median(vec1), predicted = vec1)
# (mse(actual = median(vec1), predicted = vec1)/median(vec1))*100
# mse(actual = median(vec2), predicted = vec2)
# (mse(actual = median(vec2), predicted = vec2)/median(vec2))*100
# 
# 
# df <- data.frame(vec1, vec2)
# df$orden <- 1:100
# 
# Biostatech::plotScatter(base = df, varX = "orden",
#                         varY = "vec1", adjustLine = T)
# Biostatech::plotBox(base = df, varResumen = "vec1")
# Biostatech::plotBox(base = df, varResumen = "vec2")
# Biostatech::plotScatter(base = df, varX = "orden",
#                         varY = "vec2", adjustLine = T)ç





# Old code ####


# item3 <- sapply(names(mydata), function(i) {
#   df <- mydata[[i]]
#   df <- as.data.frame(df)
#   df <- df[, c(samplesG1, samplesG2)]
#   
#   df$logFC  <- apply(df, 1, function(x) {
#     mean(x[samplesG2], na.rm =T) - mean(x[samplesG1], na.rm = T)
#   }
#   )
#   
#   df$AveExpr  <- apply(df, 1, function(x) {
#     (mean(x[samplesG2], na.rm =T) + mean(x[samplesG1], na.rm = T))/2
#   }
#   )
#   df <- df[, c("logFC", "AveExpr")]
#   
#   res <- lm(logFC~AveExpr, data = df)
#   coefPred <-  res$coefficients["AveExpr"]
#   intPred <-  res$coefficients["(Intercept)"]
#   
#   resultArea <- diffAreas(
#     intPred = intPred, 
#     coefPred = coefPred,
#     maxRange = max(df$AveExpr),
#     minRange = min(df$AveExpr),
#     intExpected = 0
#   )
#   # names(resultArea) <- i
#   return(unname(resultArea))
#   
# }, simplify = T, USE.NAMES = T
# )

