#### normScore FUNCTIONS ####

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


#...........................................................................####
# Support funcions for getting score #### 

cv <- function(x, proportion = T, na.rm = TRUE) {
  coeficiente <- (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
  if (!proportion){
    coeficiente <- coeficiente*100
  }
  return(coeficiente)
} 

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


cvGruposProt <- function(grupo, dfGrupos, dfDatos){
  
  dfGrupos <- as.data.frame(dfGrupos)
  dfDatos <- as.data.frame(dfDatos)
  espGroup <- as.vector(unlist(dfGrupos[dfGrupos$Groups == grupo, "Samples"]))
  
  df <- dfDatos[, as.character(espGroup)]
  apply(df, 1, cv, proportion = F)
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




#...........................................................................####
# SCORE: main function #### 

normScore <- function(normMatrixList, designMatrix, dfRaw, 
                      refGroup = NULL, altGroup = NULL){
  # Input: 
  # 1. List of normalized matrix (normMatrixList)
  # 2. Design matrix (designMatrix)
  # 3. Raw intensities (for item 0)
  # 4. Ref group and alternative group (optional). If they are not provided, 
  # the first group will be used as alternative and the last one, as control.
  
  totalGroups <- levels(as.factor(designMatrix$Groups))
  scoreFinal <- list()
  
  # ITEM 0 - correction factor ####
  totalIntensities <- colSums(dfRaw, na.rm = T)
  item0 <- cv(totalIntensities, proportion = T, na.rm = T)
  
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
  rankingDF[which(rownames(rankingDF) == "Log"), "Total"] <- rankingDF[which(rownames(rankingDF) == "Log"), "Total"]*item0
  
  ## Sort ####
  rankingDF <- rankingDF %>% dplyr::arrange(Total)
  finalRank <- stats::setNames(rankingDF$Total, rownames(rankingDF))
  
  
  return(list(finalRanking = finalRank, 
              detailRaking = rankingDF, 
              detailScore = scoreDF))
}

