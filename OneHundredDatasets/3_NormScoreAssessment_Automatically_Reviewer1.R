###############################################################################-

##########  Comparing normScore vs Reviewer 1 for 100 dataset   ###############-

###############################################################################-


# Julia G Curras - 2026/03/02
rm(list=ls())
graphics.off()
outDir <-  "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/ProcessedDatasets/"
setwd("C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets")

#...........................................................................####
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
library(boot)

source(file = "../R/simulationFunction.R", encoding = "UTF-8")
source(file = "../R/scoreFunction.R", encoding = "UTF-8")
setwd("C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets")


#...........................................................................####
# One hundred datasets ####
dataGS <- readxl::read_excel(path = "_Assessment_R1.xlsx", sheet = "Main", col_names = T)[-1,]
allFiles <- dataGS[, "ID", drop = T]

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
    onlyDetailRanking = T
  )$detailRanking)
},simplify = F, USE.NAMES = T)
saveRDS(object = normScoreList, file = "AssessmentFiles/3_normScore_results_100datasets.rds")



## Retrieve manual ranking - Reviewer 1  ####
### BY ITEM  ####
items <- paste0("Item", 1:6) # Sheet names = items from score
resGS <- purrr::map_dfr(items, function(sh) {
  # load
  df <- readxl::read_excel(path = "_Assessment_R1.xlsx", sheet = sh) %>%
    dplyr::select(-Number, -Final) #%>%
  
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
dfGS_Main <- readxl::read_excel(path = "_Assessment.xlsx", sheet = "Main")[-1,] 

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
table(bestNormVec)


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

saveRDS(allResultsGS, file = "AssessmentFiles/3_Reviewer1_results_100dataset.rds")



