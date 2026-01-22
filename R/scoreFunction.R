#### normScore FUNCTIONS ####

# setwd("~/GitHub/Normalization/app")
# library(shiny)
library(dplyr)
library(tidyr)
library(stats)
library(MASS)

# options(repos = BiocManager::repositories())
# library(MSnSet.utils)
library(vsn)
library(preprocessCore)
library(boot)


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
  
  # # Moving regression lines to reach B=0, a=0
  intPred <- intPred - intExpected

  if (coefPred < 0){
    a = 1
    b = -1
  } else if (coefPred > 0){
    a = -1
    b = 1
  } else if (coefPred == 0){
    return(0)
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
  
  # outra forma de calculalo que non da error cando coefPred == 0 PERO non devolve as áreas relativas e 
  # cando a pendiente é negativa devolve valores negativos...
  # recta <- function(x){coefPred*x+intPred}
  # areaMetric <- (integrate(recta, minRange, maxRange)$value)/abs(maxRange - minRange)
  # return(areaMetric)
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

rleKS <- function(dfDatos) {
  medianaProt <- apply(dfDatos, 1, stats::median, na.rm = T)
  
  rleData <- as.data.frame(log(t(t(dfDatos) / medianaProt), base = 2))
  
  ksStatistic <- apply(rleData, 2, function(i) ks.test(x = i, y = "pnorm", mean = 0, sd = 1)$statistic)
  
  return(median(ksStatistic))
}

# MAPE removing the logarithm
mape <- function(actual, predicted, prop = F){
  metric <- mean(abs((actual - predicted)/actual))
  metric <- ifelse(!prop, metric*100, metric)
  return(metric)
}

rleMAPE <- function(dfDatos) {
  medianaProt <- apply(dfDatos, 1, stats::median, na.rm = T)
  
  ## non log data
  rleData <- as.data.frame(t(t(dfDatos) / medianaProt))
  
  medianVector <- apply(rleData, 2, median, na.rm = T)
  
  return(mape(actual = 1, predicted = medianVector, prop = F))
}


tiMAPE <- function(dfDatos) {
  allMetrics <- apply(dfDatos, 2, stats::quantile, na.rm = T, simplify = F)
  q1 <- sapply(allMetrics, "[[", 2)
  q3 <- sapply(allMetrics, "[[", 4)
  medianVector <- sapply(allMetrics, "[[", 3)
  finalMetric <- 
    # mse(actual = median(medianVector), predicted = medianVector)/median(medianVector) + 
    mape(actual = median(medianVector), predicted = medianVector)/median(medianVector) + 
    # mse(actual = median(q1), predicted = q1)/median(q1) +
    mape(actual = median(q1), predicted = q1)/median(q1) +
    # mse(actual = median(q3), predicted = q3)/median(q3)
    mape(actual = median(q3), predicted = q3)/median(q3)
  
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
  
  # Shape - IQR
  # ordenar AveExp, dividir en 10 partes y estimar IQR en logFC para cada una 
  # de esas partes
  quant10 <- quantile(df$AveExpr, probs = seq(0, 1, 0.1))
  iqrList <- sapply(1:(length(quant10)-1), function(i){
    df %>% 
      dplyr::filter(AveExpr>quant10[i] & AveExpr<=quant10[i+1]) %>%
      dplyr::pull(logFC) %>% 
      stats::IQR()
  })
  # comprobar orden esperada con correlación de spearman
  rho <- cor(iqrList, 10:1, method = "spearman")
  rho <- (1 - rho) / 2 # escalado de 1 a 0, siendo 1 el peor y 0 el mejor
  CF <- 0.1 + (1 - 0.1) * rho # escalado a 0.1-1 para que no tenga tanto efecto cuando la forma se cumple
  
  # Slope - diff areas
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
  
  i3Corrected <- unname(resultArea)*CF
  return(i3Corrected)
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
#.........................................................................####

normScore <- function(
    normMatrixList, 
    designMatrix, #dfRaw, 
    refGroup = NULL, 
    altGroup = NULL, 
    onlyFinalRank = T
  ){
  # Input: 
  # 1. List of normalized matrix (normMatrixList)
  # 2. Design matrix (designMatrix)
  # 3. Raw intensities (for item 0)
  # 4. Ref group and alternative group (optional). If they are not provided, 
  # the first group will be used as alternative and the last one, as control.
  
  totalGroups <- levels(as.factor(designMatrix$Groups))
  scoreFinal <- list()
  
  # ITEM 0 - correction factor ####
  # totalIntensities <- colSums(dfRaw, na.rm = T)
  # item0 <- cv(totalIntensities, proportion = T, na.rm = T)
  
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
  # item2["CyclicLoess"] <- item2["CyclicLoess"]*1.2
  scoreFinal[["Correlation"]] <- item2*0.5
  
  
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
  
  
  # ITEM 5 - RLE: MAPE median sample (ref = 1) ####
  item5 <- sapply(normMatrixList, rleMAPE)
  # item5 <- sapply(normMatrixList, rleMSE)
  scoreFinal[["RLEplot"]] <- item5
  
  
  # ITEM 6 - total intensity: MSE median sample (ref = global median) ####
  item6 <- sapply(normMatrixList, tiMAPE)
  scoreFinal[["totalIntensity"]] <- item6
  
  
  
  # Join everything ####
  ## Bind ####
  scoreDF <- dplyr::bind_cols(scoreFinal)
  scoreDF <- as.data.frame(scoreDF)
  rownames(scoreDF) <- names(scoreFinal[[1]])
  
  # Scale ####
  scoreDF_norm <- apply(scoreDF, 2, function(col) (col - min(col)) / (max(col) - min(col)))
  scoreDF_norm <- as.data.frame(scoreDF_norm)
  
  # Corrections ####
  rownames(scoreDF_norm) <- rownames(scoreDF)
  # # 1) Small variability: no need for normalization
  # scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), ] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), ]*item0
  # # 2) CyclicLoess outstands in correlation: small correction
  # scoreDF_norm[which(rownames(scoreDF_norm) != "CyclicLoess"), 2] <- scoreDF_norm[which(rownames(scoreDF_norm) != "CyclicLoess"), 2]*0.5
  # # 3) MAD outstands in PVC: small correction
  # scoreDF_norm[which(rownames(scoreDF_norm) != "MAD"), 1] <- scoreDF_norm[which(rownames(scoreDF_norm) != "MAD"), 1]*0.8
  # # 4) Quantile outstands in TI graphics: small correction
  # scoreDF_norm[which(rownames(scoreDF_norm) != "Quantile"), 6] <- scoreDF_norm[which(rownames(scoreDF_norm) != "Quantile"), 6]*0.8

  # Rank ####
  scores_matrix <- t(scoreDF_norm)
  scoreDF_norm$Total <- rowSums(scoreDF_norm)
  # scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "Total"] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "Total"]*item0
  
  ## Sort ####
  scoreDF_norm <- scoreDF_norm %>% dplyr::arrange(Total)
  finalRank <- stats::setNames(scoreDF_norm$Total, rownames(scoreDF_norm))
  
  if (onlyFinalRank){
    return(list(finalRanking = finalRank))
  } else {
    # CI bootstrap ####
  
    # Computing total score for each normalization after resampling proteins (rows)
    bootstrap_score_rows <- function(data, indices) {
      resampled_matrix <- data[indices, , drop = FALSE]
      total_scores <- colSums(resampled_matrix)
      return(total_scores)  # One score per normalization
    }
    
    # Bootstrap
    # n_boot <- 1000
    boot_results <- boot(data = scores_matrix,              # data
                         statistic = bootstrap_score_rows,  # function for getting the scores by nomralization
                         R = 1000)                        # number of resamples  
    
    # Output mean scores and confidence intervals
    bootstrap_means <- colMeans(boot_results$t)
    
    score_bootstrap <- as.data.frame(t(sapply(1:ncol(scores_matrix), function(i){
      ci <- boot.ci(boot_results, type = "perc", index = i)
      return(c(colnames(scores_matrix)[i], bootstrap_means[i], 
               ci$percent[4], ci$percent[5]))
    }, simplify = T)))
    
    colnames(score_bootstrap) <- c("Normalization", "Mean Total Score", "LL95%", "UL95%")
    score_bootstrap <- score_bootstrap %>%
      dplyr::arrange(`Mean Total Score`)
    
    
    
    # Gráfico ####
    
    p1 <- Biostatech::plotForest(etiquetas = score_bootstrap$Normalization, 
                           estPunt = score_bootstrap$`Mean Total Score`, 
                            LI = score_bootstrap$`LL95%`, 
                            LS = score_bootstrap$`UL95%`, 
                           tituloX = "normScore with bootstrap interval")$grafico
    
    
    # Return ####
    return(list(finalRanking = finalRank, 
                detailRanking = scoreDF_norm, 
                detailScore = scoreDF, 
                bootstrapScore = score_bootstrap, 
                graphic = p1))
  }
  
}

