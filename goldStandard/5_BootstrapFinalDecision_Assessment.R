###############################################################################-

############         Refinement using boostrap decision       #################-

###############################################################################-


# Julia G Curras - 2026/03/31
rm(list=ls())
graphics.off()
outDir <-  "C:/Users/julia/Documents/GitHub/normScore-development/goldStandard/ProcessedDatasets/"
setwd("C:/Users/julia/Documents/GitHub/normScore-development/goldStandard")

#...........................................................................####
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
library(boot)

source(file = "../R/simulationFunction.R", encoding = "UTF-8")
source(file = "../R/scoreFunction.R", encoding = "UTF-8")
setwd("C:/Users/julia/Documents/GitHub/normScore-development/goldStandard")



#...........................................................................####
# Data: Gold Standard ####
dataGS <- readxl::read_excel(path = "_Assessment.xlsx", sheet = "Main", col_names = T)[-1,]
allFiles <- dataGS[, "ID", drop = T]
# allFiles <- allFiles[1:3]


## Apply normScore by dataset ####
normScoreList <- sapply(allFiles, function(i){
  cat("\t* ", i , "\n")
  output <- readRDS(file = paste0(outDir, i, ".rds"))
  return(normScore(
    normMatrixList = output$listaNorm,
    designMatrix = output$dm,
    dfRaw = output$data,
    corrected = F,
    onlyFinalRank = F,
    onlyDetailRanking = F, 
    doBootstrap = T
  ))
},simplify = F, USE.NAMES = T)
saveRDS(object = normScoreList, file = "AssessmentFiles/normScore_dataGS_All_Item0_x3_Item2_01_BOOTSTRAP_epsilon025.rds")


## Retrieve best norm ####
# SAME INFO AS THE ONE INCLUDED AT 5_EvaluatingStrangeDatasets.rmd file
allGS <- readRDS(file = "AssessmentFiles/results_Gold_Standard.rds")
resNSNoBoot <- readRDS(file = "AssessmentFiles/normScore_dataGS_All_Item0_x3_Item2_01.rds")
bestGS <- (unlist(allGS$globalBest[names(normScoreList)]))
bestNS <- sapply(lapply(normScoreList, "[[", 2), "[", 1)
bestNSNoBoost <- sapply(resNSNoBoot, function(i) rownames(i)[1])

# Medindo...
## Boost
df <- as.data.frame.matrix(table(bestGS, bestNS))
df$TI <- 0
df$Mean <- 0
df <- df[, rownames(df)]
df <- as.matrix(df)
vcd::Kappa(df)
table(bestGS == bestNS)
okID <- names((bestGS)[bestGS == bestNS])

## No boost
df <- as.data.frame.matrix(table(bestGS, bestNSNoBoost))
df$TI <- 0
df$Mean <- 0
df <- df[, rownames(df)]
df <- as.matrix(df)
vcd::Kappa(df)
table(bestGS == bestNSNoBoost)
okIDNoBoost <- names((bestGS)[bestGS == bestNSNoBoost])

## Boost vs No boost
table(okIDNoBoost %in% okID) # 38 comuns, 6 diferentes en No boost 
table(okID %in% okIDNoBoost) # 38 comuns, 6 diferentes, 6 novos
df <- as.data.frame.matrix(table(bestNS, bestNSNoBoost))
df$TI <- 0
df$Mean <- 0
df[c("TI", "Mean"), ] <- 0
df <- df[, rownames(df)]
df <- as.matrix(df)
vcd::Kappa(df)
table(bestNS == bestNSNoBoost) # hay 28 que cambian de No boost a boost



## Retrieve bestS norm ####
source(file = "3_AssessmentFunctions.R")
resAllNS <- sapply(normScoreList, "[[", 2, simplify = F, USE.NAMES = T)
summaryGlobal <- resultsGlobal(resGlobalNS = resAllNS, 
                               resGlobalGS_BEST = allGS$globalBest,
                               resGlobalGS_BESTS = allGS$globalBests,
                               resGlobalGS_WORST = allGS$globalWorst)
summaryGlobal$Best$hitTable
table(summaryGlobal$Best$mrrRes)
summaryGlobal$Bests$hitTable
table(summaryGlobal$Worst$coverageRes)
table(summaryGlobal$Worst$mrrRes)






