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


