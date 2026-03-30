
# Principio de parsimonia GO

# Julia G Currás - 2026/03/30

# Initial things ####
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")
library(dplyr)

# Data from normScore ####
resNS <- readRDS(file = "AssessmentFiles/normScore_dataGS_All_Item0_x3_Item2_01.rds")


# CI bootstrap ####

## Current code to build interval ####

# Computing total score for each normalization after resampling proteins (rows)
bootstrap_score_rows <- function(data, indices) {
  resampled_matrix <- data[indices, , drop = FALSE]
  total_scores <- colSums(resampled_matrix)
  return(total_scores)  # One score per normalization
}

# Bootstrap
scoreDF_norm <- resNS$PXD057069
item0 <- scoreDF_norm["Log", "TotalCorrected"]/scoreDF_norm["Log", "Total"]
scores_matrix <- t(scoreDF_norm[, 1:6])
scores_matrix[, "Log"] <- scores_matrix[, "Log"]*item0

# n_boot <- 1000
boot_results <- boot::boot(data = scores_matrix,              # data
                           statistic = bootstrap_score_rows,  # function for getting the scores by nomralization
                           R = 1000)                        # number of resamples  

# Output mean scores and confidence intervals
bootstrap_means <- colMeans(boot_results$t)

score_bootstrap <- as.data.frame(t(sapply(1:ncol(scores_matrix), function(i){
  ci <- boot::boot.ci(boot_results, type = "perc", index = i)
  return(c(colnames(scores_matrix)[i], bootstrap_means[i], 
           ci$percent[4], ci$percent[5]))
}, simplify = T)))

colnames(score_bootstrap) <- c("Normalization", "Mean Total Score", "LL95%", "UL95%")
score_bootstrap <- score_bootstrap %>%
  dplyr::arrange(`Mean Total Score`)


## New code to compare and decide ####

# ordenNorm <- c("Log", "Median", "Mean", "TI", "Quantile",  "CyclicLoess", "RLR", "VSN")

compareNormPairBootstrap <- function(
    bootMatrix, 
    norm1, 
    norm2, 
    epsilon = 0.5, 
    confLevel = 0.95){
  
  # Nivel de confianza
  alpha <- 1 - confLevel
  
  # Diferencias entre norm1 y norm2 de valores de bootstrap
  diffBoot <- bootMatrix[, norm1] - bootMatrix[, norm2]
  
  # Intervalo de confianza de las diferencias
  ciDiff <- stats::quantile(
    diffBoot,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE
  )
  
  # Determinando si hay diferencias o no
  noDifference <- (ciDiff[1] >= -epsilon) & (ciDiff[2] <= epsilon)
  
  # Return
  out <- list(
    meanDiff = mean(diffBoot, na.rm = TRUE),
    ll = unname(ciDiff[1]),
    ul = unname(ciDiff[2]),
    noDifference = noDifference
  )
  
  return(out)
}


selectBestNorm <- function(
    scoreBootstrap, 
    bootMatrix, 
    ordenNorm = c("Log", "Median", "Mean", "TI", "Quantile",  "CyclicLoess", "RLR", "VSN"), 
    # epsilon, 
    confLevel = 0.95
    ){
    
  # Justa a copy
    scoreDf <- scoreBootstrap
  # Just checking type of data (numeric)
    scoreDf$`Mean Total Score` <- as.numeric(scoreDf$`Mean Total Score`)
    
  # Just estimating a good epsilon - porcentaje del rango
    rangeScores <- max(scoreDf$`Mean Total Score`) - min(scoreDf$`Mean Total Score`)
    epsilon <- 0.05 * rangeScores
    
  # Selecting top norm
    currentNorm <- scoreDf$Normalization[1]
    currentIndex <- match(currentNorm, colnames(bootMatrix))
  
  # Only one normalization assessed
    if (nrow(scoreDf) == 1){
      return(currentNorm)
    }
    
  # Comparing by pairs of normalizations (from best to worst)
    for (i in 2:nrow(scoreDf)){
      
      # Next normalization (previous already detected)
      nextNorm <- scoreDf$Normalization[i]
      nextIndex <- match(nextNorm, colnames(bootMatrix))
      
      # Comparison
      comparison <- compareNormPairBootstrap(
        bootMatrix = bootMatrix,
        norm1 = currentIndex,
        norm2 = nextIndex,
        epsilon = epsilon,
        confLevel = confLevel
      )
      
      # Final decision: best and simplest normalization 
      if (comparison$noDifference){
        
        simplicityRank <- match(c(currentNorm, nextNorm), ordenNorm)
        selectedNorm <- c(currentNorm, nextNorm)[which.min(simplicityRank)]
        
        currentNorm <- selectedNorm
        currentIndex <- match(currentNorm, colnames(bootMatrix))
      }
    }
    
    return(currentNorm)
}

bootMatrix <- boot_results$t
colnames(bootMatrix) <- score_bootstrap$Normalization
bestNorm <- selectBestNorm(
  scoreBootstrap = score_bootstrap,
  bootMatrix = bootMatrix,
  # epsilon = 0.5,
  confLevel = 0.95
)
bestNorm




### Returning more info ####

compareNormPairBootstrap <- function(
    bootMatrix, 
    norm1, 
    norm2, 
    epsilon, 
    confLevel = 0.95){
  
  alpha <- 1 - confLevel
  
  diffBoot <- bootMatrix[, norm1] - bootMatrix[, norm2]
  
  ciDiff <- stats::quantile(
    diffBoot,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE
  )
  
  noDifference <- (ciDiff[1] >= -epsilon) & (ciDiff[2] <= epsilon)
  
  out <- list(
    meanDiff = mean(diffBoot, na.rm = TRUE),
    llDiff = unname(ciDiff[1]),
    ulDiff = unname(ciDiff[2]),
    noDifference = noDifference
  )
  
  return(out)
}


selectBestNormBootstrap <- function(
    scoreBootstrap, 
    bootMatrix, 
    ordenNorm, 
    epsilon = NULL, 
    confLevel = 0.95){
  
  scoreDf <- scoreBootstrap
  
  scoreDf$`Mean Total Score` <- as.numeric(scoreDf$`Mean Total Score`)
  scoreDf$`LL95%` <- as.numeric(scoreDf$`LL95%`)
  scoreDf$`UL95%` <- as.numeric(scoreDf$`UL95%`)
  
  if (is.null(epsilon)){
    # Just estimating a good epsilon - porcentaje del rango
    rangeScores <- max(scoreDf$`Mean Total Score`) - min(scoreDf$`Mean Total Score`)
    epsilon <- 0.25 * rangeScores
  }
  
  if (is.null(colnames(bootMatrix))){
    stop("bootMatrix must have column names matching the normalization names.")
  }
  
  # Asegurar que el orden del bootMatrix coincide con el de scoreDf
  bootMatrix <- bootMatrix[, scoreDf$Normalization, drop = FALSE]
  
  # Ranking inicial por score
  scoreDf$provisionalRank <- seq_len(nrow(scoreDf))
  scoreDf$simplicityRank <- match(scoreDf$Normalization, ordenNorm)
  
  # Inicialización del proceso iterativo
  currentNorm <- scoreDf$Normalization[1]
  currentIndex <- match(currentNorm, colnames(bootMatrix))
  
  comparisonTraceList <- vector("list", max(0, nrow(scoreDf) - 1))
  
  if (nrow(scoreDf) > 1){
    
    for (i in 2:nrow(scoreDf)){
      
      nextNorm <- scoreDf$Normalization[i]
      nextIndex <- match(nextNorm, colnames(bootMatrix))
      
      currentMean <- scoreDf$`Mean Total Score`[match(currentNorm, scoreDf$Normalization)]
      nextMean <- scoreDf$`Mean Total Score`[i]
      
      comparison <- compareNormPairBootstrap(
        bootMatrix = bootMatrix,
        norm1 = currentIndex,
        norm2 = nextIndex,
        epsilon = epsilon,
        confLevel = confLevel
      )
      
      if (comparison$noDifference){
        
        simplicityRank <- match(c(currentNorm, nextNorm), ordenNorm)
        selectedNorm <- c(currentNorm, nextNorm)[which.min(simplicityRank)]
        decisionReason <- "Equivalent within epsilon; selected simpler normalization"
        
      } else {
        selectedNorm <- currentNorm
        decisionReason <- "Different beyond epsilon; retained lower-score normalization"
      }
      
      comparisonTraceList[[i - 1]] <- data.frame(
        step = i - 1,
        currentNorm = currentNorm,
        nextNorm = nextNorm,
        currentMean = currentMean,
        nextMean = nextMean,
        meanDiff = comparison$meanDiff,
        llDiff = comparison$llDiff,
        ulDiff = comparison$ulDiff,
        epsilon = epsilon,
        noDifference = comparison$noDifference,
        selectedNorm = selectedNorm,
        decisionReason = decisionReason,
        stringsAsFactors = FALSE
      )
      
      currentNorm <- selectedNorm
      currentIndex <- match(currentNorm, colnames(bootMatrix))
    }
  }
  
  bestNormalization <- currentNorm
  winnerIndex <- match(bestNormalization, colnames(bootMatrix))
  
  # Resumen frente a la ganadora final
  diffSummaryList <- lapply(seq_len(ncol(bootMatrix)), function(j){
    
    diffBoot <- bootMatrix[, j] - bootMatrix[, winnerIndex]
    
    ciDiff <- stats::quantile(
      diffBoot,
      probs = c((1 - confLevel) / 2, 1 - (1 - confLevel) / 2),
      na.rm = TRUE
    )
    
    data.frame(
      Normalization = colnames(bootMatrix)[j],
      meanDiffVsWinner = mean(diffBoot, na.rm = TRUE),
      llDiffVsWinner = unname(ciDiff[1]),
      ulDiffVsWinner = unname(ciDiff[2]),
      equivalentToWinner = (ciDiff[1] >= -epsilon) & (ciDiff[2] <= epsilon),
      stringsAsFactors = FALSE
    )
  })
  
  diffSummaryDf <- do.call(rbind, diffSummaryList)
  
  normSummary <- scoreDf %>%
    dplyr::left_join(diffSummaryDf, by = "Normalization") %>%
    dplyr::mutate(
      selectedFinal = Normalization == bestNormalization
    )
  
  comparisonTrace <- if (length(comparisonTraceList) > 0){
    do.call(rbind, comparisonTraceList)
  } else {
    data.frame(
      step = numeric(0),
      currentNorm = character(0),
      nextNorm = character(0),
      currentMean = numeric(0),
      nextMean = numeric(0),
      meanDiff = numeric(0),
      llDiff = numeric(0),
      ulDiff = numeric(0),
      epsilon = numeric(0),
      noDifference = logical(0),
      selectedNorm = character(0),
      decisionReason = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  out <- list(
    bestNormalization = bestNormalization,
    normSummary = normSummary,
    comparisonTrace = comparisonTrace
  )
  
  return(out)
}


#### Uso ####

ordenNorm <- c("Log", "Median", "Mean", "TI", "Quantile", "CyclicLoess", "RLR", "VSN")

bootMatrix <- boot_results$t
colnames(bootMatrix) <- score_bootstrap$Normalization
bestNormResults <- selectBestNormBootstrap(
  scoreBootstrap = score_bootstrap,
  bootMatrix = bootMatrix,
  ordenNorm = ordenNorm,
  # epsilon = 0.10,
  confLevel = 0.95
)
normSummary <- bestNormResults$normSummary
comparisonTrace <- bestNormResults$comparisonTrace
bestNormResults$bestNormalization


# Representación 1 - forest plot normSummary
library(ggplot2)
ggplot(normSummary, aes(x = `Mean Total Score`, y = reorder(Normalization, `Mean Total Score`))) +
  geom_errorbarh(aes(xmin = `LL95%`, xmax = `UL95%`, color = selectedFinal), height = 0.2) +
  geom_point(aes(color = selectedFinal, shape = equivalentToWinner), size = 3) +
  labs(
    x = "Mean total score",
    y = "Normalization",
    color = "Selected final",
    shape = "Equivalent to winner"
  ) +
 ggplot2::theme_minimal() + 
 ggplot2::theme(text = element_text(family = "Calibri", size = 14), 
 axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
 axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))


# Representación 2 - forest plot comparisonTrace
ggplot(
  comparisonTrace,
  aes(x = meanDiff, y = reorder(paste0("Step ", step, ": ", currentNorm, " vs ", nextNorm), step))
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(aes(xintercept = epsilon), linetype = "dotted") +
  geom_vline(aes(xintercept = -epsilon), linetype = "dotted") +
  geom_segment(
    aes(
      x = llDiff,
      xend = ulDiff,
      yend = reorder(paste0("Step ", step, ": ", currentNorm, " vs ", nextNorm), step),
      color = noDifference
    ),
    linewidth = 0.8
  ) +
  geom_point(aes(color = noDifference, shape = selectedNorm), size = 3) +
  labs(
    x = "Bootstrap mean difference",
    y = "Comparison",
    color = "Equivalent within epsilon",
    shape = "Selected normalization"
  ) +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = element_text(family = "Calibri"), 
                 axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
                 axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

