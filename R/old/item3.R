##############################################################################- 

# Julia G Curras - 21/01/2026

# Evaluación de items individuais do score
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")

library(dplyr)
library(ggplot2)
library(tictoc)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")
source(file = "simulationFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/others/"



# Data 
sigma_lo_vec <- seq(0, 1.5, 0.1875)
sigma_hi_vec <- rev(seq(0, 1.5, 0.1875))*(1.5)
valores_sample_sd_cap <- rev(seq(0, 1.5, 0.1875))
# valores_sample_sd_cap <- seq(0.5, 2, 0.1875)
sigma_lo_vec <- c(0.4, 0.4, 0.4, 0.4, 0.6, 2.9, 3.2, 0.8, 1.5)
sigma_hi_vec <- c(0.05, 0.05, 0.05, 0.05, 0.9, 3.9, 3.2, 2.25, 1.2)
valores_sample_sd_cap <- c(0.5, 1, 1.25, 2, 1, 1.5, 1.9, 1, 1, 0.35)
values_sample_sd_strength <- c(3, 3, 3, 3, rep(1, 3), 0.8, 0)
values_sample_sd_rho <- c(0, 0, 0, 0, 0.8, 0.8, 0.8, 0.7, 0.8)



sigma_lo_vec <-c(3, 2, 1, 0.5, 0.3, 0.3, 0.3, 0.3, 0.1)
sigma_hi_vec <-c(0.05, 0.1, 0.1, 0.4, 0.8, 1.2, 2, 2.5, 4)
valores_sample_sd_cap <- rep(0.2, 9)
values_sample_sd_strength <- rep(0.8, 9)
values_sample_sd_rho <-  rep(0.35, 9)


resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
  semilla = 9693, 
  n_proteins = 10000,
  n_per_group = 15, 
  prop_de = 0.1,
  sigma_lo = sigma_lo_vec[val],
  sigma_hi = sigma_hi_vec[val], 
  sample_sd_strength = values_sample_sd_strength[val],
  sample_sd_rho = values_sample_sd_rho[val],
  sample_sd_cap = valores_sample_sd_cap[val]
))
getResultsByItem(resByItem, "item3")


# ITEM 3 - MAplot regression line 0 ####
dm <- as.data.frame(resByItem[[1]][["metadata"]])
finalData <- sapply(resByItem, "[[", 1, simplify = F)
totalGroups <- unique(dm$Groups)
refGroup <- "G1"
altGroup <- "G2"
samplesG1 <- dm[dm$Groups == refGroup, "Samples"]
samplesG2 <- dm[dm$Groups == altGroup, "Samples"]

item3 <- sapply(finalData, maDiffAreas, samplesG1 = samplesG1, 
                samplesG2 = samplesG2, simplify = T)
item3

data <- finalData[[1]]

order(item3)



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

  item3 <- unname(resultArea)
  i3Corrected <- item3*CF
  return(i3Corrected)
}






