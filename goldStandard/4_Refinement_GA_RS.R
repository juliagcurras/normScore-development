###############################################################################-

############             Refinement using weights             #################-

###############################################################################-


# Julia G Curras - 2026/03/11
rm(list=ls())
graphics.off()
outDir <-  "C:/Users/julia/Documents/GitHub/normScore/goldStandard/ProcessedDatasets/"
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")

source(file = "../R/scoreFunction.R", encoding = "UTF-8")

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
    onlyDetailRanking = F, 
    detailRankItem0 = T
  ))
},simplify = F, USE.NAMES = T)
saveRDS(object = normScoreList, file = "AssessmentFiles/Refinement_normScore_dataGS_All.rds")
# normScoreList <- readRDS(file = "AssessmentFiles/normScore_dataGS_All.rds")


## Select best normalization from GS ####
normalizationNames <- c("Log", "Mean", "TI", "Median", "Quantile", "CyclicLoess", 
                        "RLR", "VSN")

# Best UNIQUE normalization
dfBest <- dataGS %>%
  dplyr::select(ID, normScore) %>%
  tibble::deframe()
resGlobalGS_BEST <- sapply(dfBest, strsplit, split = ";", fixed = T,
                           simplify = T, USE.NAMES = T)
bestNormVec <- sapply(resGlobalGS_BEST, "[[", 1, simplify = F)
bestNormVec <- lapply(bestNormVec, function(x) list(solution = x))
resGlobalGS_BEST <- lapply(resGlobalGS_BEST, function(x) list(solutions = x))
if (!all(sapply(bestNormVec, function(i) any(i %in% normalizationNames)))){
  stop("Algún nombre de normalization no coincide!!")
}

listaFinal <- Map(c, normScoreList, bestNormVec)
listaFinal <- Map(c, listaFinal, resGlobalGS_BEST)
saveRDS(object = listaFinal, file = "AssessmentFiles/Refinement_inputData.rds")


#...........................................................................####
# Auxiliar functions ####
maxFitness <- function(listDatasets, fitness = 1){
  if (fitness == 1){
    maxFitness = length(listDatasets)
  } else if (fitness == 2){
    maxFitness = length(listDatasets)*(1.001)
  } else if (fitness == 3){
    maxFitness = length(listDatasets)*(1.1)
  }
  return(maxFitness)
} 

normScoreAssembly <- function(matriz, item0, itemWeights){
  scoreDF_norm <- sweep(matriz, 2, itemWeights, `*`)

  # Rank #
  scores_matrix <- t(scoreDF_norm)
  scoreDF_norm$Total <- rowSums(scoreDF_norm)
  scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "Total"] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "Total"]*item0
  
  # Sort #
  ordenFilas <- c("Log", "Median", "Mean", "TI", "Quantile", "CyclicLoess", 
                  "RLR", "VSN")
  scoreDF_norm <- scoreDF_norm %>%
    dplyr::arrange(Total, match(rownames(.), ordenFilas)) # en caso de empates, se ordena por sencillez
  
  # Return better
  return(rownames(scoreDF_norm))

}

# normScoreDecision <- function(normScoreMatriz){
#   # Sort #
#   ordenFilas <- c("Log", "Median", "Mean", "TI", "Quantile", "CyclicLoess", 
#                   "RLR", "VSN")
#   scoreDF_norm <- normScoreMatriz %>%
#     dplyr::arrange(Total, match(rownames(.), ordenFilas)) # en cso de empates, se ordena por sencillez
#   
#   return(rownames(scoreDF_norm)[1])
# }


softmaxWeights  <- function(weights, fitness = 1) {
  # O GA() xera números reales que nos adaptamos ás restriccions dos pesos: sumar 1 e ser positivas
  weights <- weights - max(weights) # Esto é solo para que non de +Inf con valores moi altos
  weights <- exp(weights)
  weights <- weights / sum(weights)
  return(weights)
}


# assessNormScore <- function(matriz, item0, itemWeights, solution){
assessNormScore <- function(datasetList, itemWeights, fitness = 1){
  # Get info
  matriz <- datasetList$matriz
  item0 <- datasetList$item0
  solution <- datasetList$solution
  solutions <- datasetList$solutions
  
  # Calcular normScore con pesos e obter mellor normalización
  rankedNorms <- normScoreAssembly(matriz = matriz, item0 = item0, 
                                   itemWeights = itemWeights)
  finalNorm <- rankedNorms[1]
  
  
  if (fitness == 1){
    #----------------------- Fitness 1
    ## C1: match único
    C1 <- as.integer(finalNorm == solution)
    ## C2: if C1 == 0, is norm in the list of other better norms? (desempate)
    C2 <- as.integer((finalNorm %in% solutions) & (finalNorm != solution))
    ## Combinamos
    fitness <- C1 + 0.001*C2 # Si C1 = 0 (non hai coincidencia), pero C2 == 1, mellor peso
  } else if (fitness == 2){
    #----------------------- Fitness 2
    C2 <- as.integer(finalNorm %in% solutions)
    C1 <- as.integer(finalNorm == solution)
    fitness <- C2 + 0.001 * C1
  } else if (fitness == 3){
    ## C1: a mellor normalización segundo normScore está na lista de mellores do gold standard
    C1 <- as.integer(finalNorm %in% solutions)
    ## C2: a primeira opción do gold standard aparece no top 3 do ranking segundo normScore
    rankPos <- match(solution, rankedNorms)
    if (is.na(rankPos) || rankPos > 3) {
      C2 <- 0
    } else {
      C2 <- (4 - rankPos) / 3
    }
    ## Combinación
    fitness <- C1 + 0.1 * C2
  }
  
  return(fitness)
}


generateRandomWeights <- function(nWeights = 6) {
  x <- runif(nWeights)
  x / sum(x)
}

# Main function GA ####

getFitness <- function(
    itemWeights, 
    fitness = 3
){
  # transform weights
  weights <- softmaxWeights(itemWeights)
  sum(weights)
  # estimate fitness 
  allFitness <- sapply(datasetsList, assessNormScore, itemWeights = weights,
                       fitness = fitness, simplify = T, USE.NAMES = F)
  sumFitness <- sum(allFitness)
  
  # Return fitness
  return(sumFitness)
}


# Main function RS ####

randomSearchWeights <- function(nIter, nWeights = 6, seed = 123, 
                                nCores = 4, fitness = 1) {
  
  set.seed(seed) # general seed
  seeds <- sample.int(.Machine$integer.max, nIter) # Creating random seeds for each iteration
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(future::multisession, workers = nCores)
  
  cat("Starting iterations...\n")
  results <- future.apply::future_lapply( # lapply in parallel -> decreasing execution time
    X = seq_len(nIter), # for each iterarion
    FUN = function(i) {
      set.seed(seeds[i]) # Set a seed
      weights <- generateRandomWeights(nWeights = nWeights) # generate weights
      fitness <- getFitness(weights, fitness = fitness) # estimate fitness in 100 datasets
      list(weights = weights, fitness = fitness) # return
    },
    future.seed = TRUE 
  )
  
  cat("Selecting best iteration...\n")
  # generate a vector just for fitness results
  fitnessVec <- vapply(results, function(x) x$fitness, numeric(1))
  bestIdx <- which.max(fitnessVec) # get best weight
  
  cat("End\n")
  list(
    bestWeights = results[[bestIdx]]$weights,
    bestFitness = results[[bestIdx]]$fitness,
    fitnessHistory = fitnessVec
  )
}


#...........................................................................####
# Genetic algorithm ejecution ####
# Data required #
allDatasetsList <- readRDS(file = "AssessmentFiles/Refinement_inputData.rds") # IMPORTANTE QUE SE CHAME datasetsList!!!


# Applying train-test #
  # Generation data partition 70 - 30
set.seed(1165)
n <- names(allDatasetsList)
idTrain <- sample(x = n, size = 70, replace = F)
idTest <- n[!(n %in% idTrain)]
any(idTest %in% idTrain) # just a checkpoint
length(idTest) + length(idTrain) == 100 # just another checkpoint

  # Selecting datasets for each part
datasetsListTrain <- allDatasetsList[idTrain]
datasetsListTest <- allDatasetsList[idTest]


# Execute for train
# Execute with different fitness values (changing default value for fitness 
# parameter in getFitness function)
library(GA)
datasetsList <- datasetsListTrain # IMPORTANTE QUE O OBXECTO COS DATOS SE CHAME datasetsList!!!
gaModel <- ga(
  type = "real-valued",
  fitness = getFitness,
  lower = rep(-5, 6), # limitación espacio de numeros reales a generar (mínimo por peso)
  upper = rep(5, 6), # limitación espacio de numeros reales a generar (máximos por peso)
  popSize = 200,
  maxiter = 500,
  run = 50,
  pcrossover = 0.8,
  pmutation = 0.1,
  elitism = 10,
  monitor = TRUE,
  seed = 9396,
  parallel = T 
)
tipoFitness <- 3 # cambiar en getFitness manualmente
# gaModel <- readRDS(file = "AssessmentFiles/Refinement_results_GA_FITNESS3_TRAIN.rds") # 500 iteracións

# Assess for test
summary(gaModel)
softmaxWeights(gaModel@solution[1,])
plot(gaModel)
fitnessTRAIN <- gaModel@fitnessValue/maxFitness(listDatasets = datasetsListTrain, 
                                                fitness = tipoFitness) # fitness relativo

# Assess for train
datasetsList <- datasetsListTest
pesos <- softmaxWeights(gaModel@solution[1,])
res <- getFitness(itemWeights = pesos, fitness = tipoFitness)
fitnessTEST <-res/maxFitness(listDatasets = datasetsList, fitness = tipoFitness)

# Saving
resGA <- list(
  model = gaModel, 
  fitnessTrain = fitnessTRAIN,
  fitnessTest = fitnessTEST
)
saveRDS(resGA, file = paste0("AssessmentFiles/Refinement_results_GA_FITNESS", 
                             tipoFitness, "_TRAIN.rds")) # 500 iteracións


# Assess all
datasetsList <- allDatasetsList
pesos <- softmaxWeights(gaModel@solution[1,])
res <- getFitness(itemWeights = pesos, fitness = tipoFitness)
res/maxFitness(listDatasets = datasetsList, fitness = tipoFitness)


#...........................................................................####
# Random search ####
tictoc::tic()
randomResult <- randomSearchWeights(
  nIter = 50000, # 1000 iterations = 49.5s // 50000 iteracións 1609 s
  nWeights = 6,
  seed = 1409,
  nCores = 8, 
  fitness = 2 # 1
)
tictoc::toc()

randomResult$bestWeights
randomResult$bestFitness

saveRDS(randomResult, file = "AssessmentFiles/Refinement_results_RS_FITNESS2.rds")  # 50000 iteracións
#...........................................................................####


