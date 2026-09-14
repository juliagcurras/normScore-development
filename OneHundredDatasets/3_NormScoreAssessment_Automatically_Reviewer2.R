###############################################################################-

############             Refinement using weights             #################-

###############################################################################-


# Julia G Curras - 2026/05/25
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
dataGS <- readxl::read_excel(path = "_Assessment_R2.xlsx", sheet = "Main", col_names = T)[-1,]
allFiles <- dataGS[, "ID", drop = T]

## Apply normScore by dataset ####
# Done at R1 script
# normScoreList <- sapply(allFiles, function(i){
#   cat("\t* ", i , "\n")
#   output <- readRDS(file = paste0(outDir, i, ".rds"))
#   return(normScore(
#     normMatrixList = output$listaNorm,
#     designMatrix = output$dm,
#     dfRaw = output$data,
#     corrected = F,
#     onlyFinalRank = F,
#     onlyDetailRanking = T
#   )$detailRanking)
# },simplify = F, USE.NAMES = T)
# saveRDS(object = normScoreList, file = "AssessmentFiles/normScore_dataR2_All_Item0_x4_Item2_01.rds")



## Retrieve manual ranking GS ####
### BY ITEM  ####
items <- paste0("Item", 1:6) # Sheet names = items from score
resGS <- purrr::map_dfr(items, function(sh) {
  # load
  df <- readxl::read_excel(path = "_Assessment_R2.xlsx", sheet = sh) %>%
    dplyr::select(-Number, -Final) #%>%
  # slice(1:50) %>%
  # filter(ID != "PXD005025")# Number col
  
  # change format
  df %>%
    tidyr::pivot_longer(
      cols = -ID,
      names_to = "Normalization",
      values_to = "RankGS"
    ) %>%
    dplyr::mutate(
      Item = sh,
      RankGS = as.numeric(RankGS)
    ) %>%
    dplyr::select(ID, Item, Normalization, RankGS) %>%
    dplyr::arrange(ID, Item, RankGS)
})


### GLOBAL  ####
normalizationNames <- c("Log", "Mean", "Median","TI", "VSN",
                        "Quantile", "CyclicLoess", "RLR")
dfGS_Main <- readxl::read_excel(path = "_Assessment_R2.xlsx", sheet = "Main")[-1,] 

# Best UNIQUE normalization
dfBest <- dfGS_Main %>%
  dplyr::select(ID, normScore) %>%
  tibble::deframe()
resGlobalGS_BEST <- sapply(dfBest, strsplit, split = ";", fixed = T,
                           simplify = T, USE.NAMES = T)
resGlobalGS_BEST <- sapply(resGlobalGS_BEST, "[[", 1, simplify = F)
if (!all(sapply(resGlobalGS_BEST, function(i) any(i %in% normalizationNames)))){
  stop("Algún nombre de normalization no coincide!!")
}
bestNormVec <- resGlobalGS_BEST
Biostatech::getCatTable(df = as.data.frame(bestNormVec))$tablaFormato



# Best normalizations
dfBestS <- dfGS_Main %>%
  dplyr::select(ID, normScore) %>%
  tibble::deframe()
resGlobalGS_BESTS <- sapply(dfBest, strsplit, split = ";", fixed = T,
                            simplify = T, USE.NAMES = T)
if (!all(sapply(resGlobalGS_BESTS, function(i) any(i %in% normalizationNames)))){
  stop("Algún nombre de normalization no coincide!!")
}

# Excluded normalizations
dfWorst <- dfGS_Main %>%
  dplyr::select(ID, Excluded) %>%
  filter(!is.na(Excluded)) %>%
  tibble::deframe()
resGlobalGS_WORST <- sapply(dfWorst, strsplit, split = ";", fixed = T,
                            simplify = T, USE.NAMES = T)
if (!all(sapply(resGlobalGS_WORST, function(i) any(i %in% normalizationNames)))){
  stop("Algún nombre de normalization no coincide!!")
}


allResultsGS <- list(
  byItem = resGS, 
  globalBest = resGlobalGS_BEST, 
  globalBests = resGlobalGS_BESTS, 
  globalWorst = resGlobalGS_WORST
)

saveRDS(allResultsGS, file = "AssessmentFiles/results_Gold_Standard_R2.rds")





#...........................................................................####
# Assessment ####
source(file = "3_AssessmentFunctions.R")

# SAME INFO AS THE ONE INCLUDED AT 3_NormScoreAssesment.rmd file

## All ####

### Format data ####
resAllGS <- readRDS(file = "AssessmentFiles/results_Gold_Standard.rds") # Manual ranking
resNS <- readRDS(file = "AssessmentFiles/normScore_dataGS_All_Item0_x4_Item2_01.rds") # NS results non corrected

resAllNS <- formatNormScoreResults(finalList = resNS, item0 = T)

### By item ####
summaryByItem <- resultsByItem(resNormScore = resAllNS$resByItem, 
                                        resGS = resAllGS$byItem)
summaryByItem$tabla %>% View

### Global ####
# Non corrected
summaryGlobal <- resultsGlobal(resGlobalNS = resAllNS$resGlobal, 
                               resGlobalGS_BEST = resAllGS$globalBest,
                               resGlobalGS_BESTS = resAllGS$globalBests,
                               resGlobalGS_WORST = resAllGS$globalWorst)
summaryGlobal$Best$hitTable
summaryGlobal$Bests$hitTable
table(summaryGlobal$Worst$coverageRes)
table(summaryGlobal$Worst$mrrRes)



#