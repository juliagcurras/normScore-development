###############################################################################-

############             Refinement using weights             #################-

###############################################################################-


# Julia G Curras - 2026/03/02
rm(list=ls())
graphics.off()
outDir <-  "C:/Users/julia/Documents/GitHub/normScore/goldStandard/ProcessedDatasets/"
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")

#...........................................................................####
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
library(boot)

source(file = "../R/simulationFunction.R", encoding = "UTF-8")
source(file = "../R/scoreFunction.R", encoding = "UTF-8")
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")


#...........................................................................####
# Data: Gold Standard ####
dataGS <- readxl::read_excel(path = "_Assessment.xlsx", sheet = "Main", col_names = T)[-1,]
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
saveRDS(object = normScoreList, file = "AssessmentFiles/normScore_dataGS_All_Item0_gamma_09_12.rds")
# saveRDS(object = normScoreList, file = "AssessmentFiles/normScore_dataGS_All_AnyCorrection.rds")
# saveRDS(object = normScoreList, file = "AssessmentFiles/normScore_dataGS_All_suavizadoItem0_2.rds")


## Select best normalization from GS ####
normalizationNames <- c("Log", "Mean", "Median","TI", "VSN",
                        "Quantile", "CyclicLoess", "RLR", "MAD")

# Best UNIQUE normalization
dfBest <- dataGS %>%
  dplyr::select(ID, normScore) %>%
  tibble::deframe()
resGlobalGS_BEST <- sapply(dfBest, strsplit, split = ";", fixed = T,
                           simplify = T, USE.NAMES = T)
bestNormVec <- sapply(resGlobalGS_BEST, "[[", 1, simplify = T)
if (!all(sapply(bestNormVec, function(i) any(i %in% normalizationNames)))){
  stop("Algún nombre de normalization no coincide!!")
}

bestNormVec
Biostatech::getCatTable(df = as.data.frame(bestNormVec))$tablaFormato



#...........................................................................####
# Format data ####

## Just binding normScore results from each datasets ####
metricsLong <- purrr::imap_dfr(normScoreList, function(mat, i) { # pasteing all
  df <- as.data.frame(mat)
  df <- df[, 1:6, drop = FALSE]
  normId <- rownames(mat)
  
  df %>%
    dplyr::mutate(
      datasetId = i,
      normId = normId
    ) %>%
    relocate(datasetId, normId)
})

# New column adding info about best normalization  
metricsLong <- metricsLong %>%
  dplyr::group_by(datasetId) %>%
  dplyr::mutate(
    isWinner = {
      best <- bestNormVec[unique(datasetId)]
      if (is.numeric(best)) normId == as.character(best) else normId == best
    }
  ) %>%
  dplyr::ungroup()
table(metricsLong$isWinner) # Comprobando que hay 7*100 normalizacions non mellores e 100 que son a mellor para cada datasets


## Calculating input data for model ####
itemNames <- setdiff(names(metricsLong), c("datasetId","normId","isWinner"))

# Para cada dataset, restamos as puntuacións da normalización ganadora menos as 
# do resto, por item. Despois, engadimos unha columna onde decimos que eses valores
# teóricamente negativos que saen da resta están asociados a gañar, que sería 1 
# chat di que non fai falta que haxa dúas clases para axustar a loxísica porque 
# esto é un problema de buscar pesos óptimos. Polo tanto, obtemos un dataset final
# para pasarlle ao glm de 7*100=700 filas. Chat tamén di que o tema de que 7 filas 
# pertenzan ao mesmo dataset non importa (certa dependencia de puntuacións?). 
# Por que se fai isto da resta e non se usan os valores tal cual dos scores? Para 
# relativizar respecto ao rango de cada score (puntuación de 0.15 a 0 nun dataset, 
# e un item, que pode ser de 0.7 a 0 para outro).
pairsDf <- metricsLong %>%
  dplyr::group_by(datasetId) %>%
  dplyr::group_modify(~{
    df <- .x
    winnerRow <- df %>% filter(isWinner)
    if (nrow(winnerRow) != 1) stop("En algún dataset no hay exactamente 1 winner.")
    rivals <- df %>% filter(!isWinner)
    # Restas winner - rival para cada ítem
    out <- rivals %>%
      dplyr::mutate(across(all_of(itemNames), ~ winnerRow[[cur_column()]] - .x, .names = "d{.col}")) %>%
      dplyr::transmute(
        datasetId = winnerRow$datasetId,
        rivalNormId = normId,
        across(starts_with("d")),
        y = 1L
      )
    out
  }) %>%
  ungroup()


#...........................................................................####
# Model adjustment ####
fit <- glm(
  y ~ 0 + dItem1 + dItem2 + dItem3 + dItem4 + dItem5 + dItem6,
  data = pairsDf,
  family = binomial()
)
summary(fit)
betas <- coef(fit) # coefficients

# Pesos positivos + suman 1 (interpretables)
weightsRaw <- abs(betas)
itemWeights <- weightsRaw / sum(weightsRaw)

itemWeights

itemWeights <- c(0.3, 0.05, 0.3, 0.15, 0.05, 0.15)
sum(itemWeights)
saveRDS(object = itemWeights, file = "AssessmentFiles/weigths_All.rds")


# Just checking results
normScoreListCorrected <- sapply(allFiles, function(i){
  cat("\t* ", i , "\n")
  output <- readRDS(file = paste0(outDir, i, ".rds"))
  return(normScore(
    normMatrixList = output$listaNorm,
    designMatrix = output$dm,
    dfRaw = output$data,
    corrected = T,
    onlyFinalRank = F,
    onlyDetailRanking = T
  )$detailRanking)
},simplify = F, USE.NAMES = T)
saveRDS(object = normScoreListCorrected, file = "AssessmentFiles/normScore_CORRECTED_dataGS_All_ManualWeights.rds")
# normScoreListCorrected <- readRDS(file = "AssessmentFiles/normScore_CORRECTED_dataGS_All_ModItem2_ManualWeights.rds")
# saveRDS(object = normScoreListCorrected, file = "AssessmentFiles/normScore_CORRECTED_dataGS_All.rds")
# normScoreListCorrected <- readRDS(file = "AssessmentFiles/normScore_CORRECTED_dataGS_All.rds")


## Repeated CV ####
fitModel <- function(pairsDF){
  fit <- glm(
    y ~ 0 + dItem1 + dItem2 + dItem3 + dItem4 + dItem5 + dItem6,
    data = pairsDF,
    family = binomial()
  )
  summary(fit)
  betas <- coef(fit) # coefficients
  
  # Pesos positivos + suman 1 (interpretables)
  weightsRaw <- abs(betas)
  itemWeights <- weightsRaw / sum(weightsRaw)
  
  return(itemWeights)
}


### Train - test split ####
set.seed(9396)
idsTrain <- sample(x = dataGS$ID, size = 70, replace = F)
idsTest <- dataGS$ID[!(dataGS$ID %in% idsTrain)]
pairsDfTrain <- pairsDf %>% dplyr::filter(datasetId %in% idsTrain)
unique(pairsDfTrain$datasetId)

### Folds ####
set.seed(9396)
folds <- caret::createMultiFolds(unique(pairsDfTrain$datasetId), k = 3, times = 20)

## RUN! ####
cvWeigths <- sapply(folds, function(x){
  idsFold <- unique(pairsDfTrain$datasetId)[x] # selecting datasets ID
  training_fold <- pairsDfTrain %>% dplyr::filter(datasetId %in% idsFold) # Selecting rows from included datasets
  pesos <- fitModel(pairsDF = training_fold) # adjusting model and extracting weights
  return(pesos)
}, simplify = T, USE.NAMES = T)

### Show results ####
cvWeigths <- as.data.frame(t(cvWeigths))
cvWeigths[,1:ncol(cvWeigths)] <- apply(cvWeigths[,1:ncol(cvWeigths)], 2, as.numeric)

# Summary #
res <- apply(cvWeigths, 2, Biostatech::getMeanIC)
res <- sapply(res, "[[", 1, simplify = T)
rownames(res) <- c("Mean", "SD", "N", "Quantile t", "IL", "SL")
res

# Graphical representation #
Biostatech::plotForest(etiquetas = colnames(res),
                       estPunt = res["Mean",], 
                       LI = res["IL", ], 
                       LS = res["SL", ])$grafico

Biostatech::plotBoxMultivar(base = cvWeigths, varResumen = colnames(cvWeigths))$grafico





#...........................................................................####
# Assessment ####
source(file = "3_AssessmentFunctions.R")

## All ####

### Format data ####
resAllGS <- readRDS(file = "AssessmentFiles/results_Gold_Standard.rds") # Manual ranking
resNS <- readRDS(file = "AssessmentFiles/normScore_dataGS_All_Item0_suavizado_alpha09.rds") # NS results non corrected
resNSCorrected <- readRDS(file = "AssessmentFiles/normScore_CORRECTED_dataGS_All_MOdItem2_ManualWeights.rds") # NS results corrected

resAllNS_Corrected <- formatNormScoreResults(finalList = resNSCorrected, item0 = T)
resAllNS <- formatNormScoreResults(finalList = resNS, item0 = T)

### By item ####
# Same between corrected and non corrected
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

# Corrected
summaryGlobalCorrected <- resultsGlobal(resGlobalNS = resAllNS_Corrected$resGlobal, 
                                        resGlobalGS_BEST = resAllGS$globalBest,
                                        resGlobalGS_BESTS = resAllGS$globalBests,
                                        resGlobalGS_WORST = resAllGS$globalWorst)
summaryGlobalCorrected$Best$hitTable
summaryGlobalCorrected$Bests$hitTable
table(summaryGlobalCorrected$Worst$coverageRes)
table(summaryGlobalCorrected$Worst$mrrRes)


summaryGlobal$Best$hitTable
summaryGlobalCorrected$Best$hitTable

## Training ####


## Test ####


## Non-weitghted score ####



#


































