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
saveRDS(object = normScoreList, file = "AssessmentFiles/GA_normScore_dataGS_All.rds")
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
saveRDS(object = listaFinal, file = "AssessmentFiles/GA_inputData.rds")


#...........................................................................####
# Auxiliar functions ####

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
    dplyr::arrange(Total, match(rownames(.), ordenFilas)) # en cso de empates, se ordena por sencillez
  
  # Return better
  return(rownames(scoreDF_norm)[1])

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


softmaxWeights  <- function(weights) {
  # O GA() xera números reales que nos adaptamos ás restriccions dos pesos: sumar 1 e ser positivas
  weights <- weights - max(weights) # Esto é solo para que non de +Inf con valores moi altos
  weights <- exp(weights)
  weights <- weights / sum(weights)
  return(weights)
}


# assessNormScore <- function(matriz, item0, itemWeights, solution){
assessNormScore <- function(datasetList, itemWeights){
  # Get info
  matriz <- datasetList$matriz
  item0 <- datasetList$item0
  solution <- datasetList$solution
  solutions <- datasetList$solutions
  
  # Calcular normScore con pesos e obter mellor normalización
  finalNorm <- normScoreAssembly(matriz = matriz, item0 = item0, 
                                   itemWeights = itemWeights)
  
  # Avaliar
  ## C1: match único
  C1 <- as.integer(finalNorm == solution)
  
  ## C2: if C1 == 0, is norm in the list of other better norms? (desempate)
  C2 <- as.integer((finalNorm %in% solutions) & (finalNorm != solution))
  
  ## Combinamos
  fitness <- C1 + 0.001*C2 # Si C1 = 0 (non hai coincidencia), pero C2 == 1, mellor peso
  
  return(fitness)
  
}


generateRandomWeights <- function(nWeights = 6) {
  x <- runif(nWeights)
  x / sum(x)
}

# Main function GA ####

getFitness <- function(
    itemWeights
){
  # transform weights
  weights <- softmaxWeights(itemWeights)
  sum(weights)
  # estimate fitness 
  allFitness <- sapply(datasetsList, assessNormScore, itemWeights = weights, 
                       simplify = T, USE.NAMES = F)
  sumFitness <- sum(allFitness)
  
  # Return fitness
  return(sumFitness)
}


# Main function RS ####

randomSearchWeights <- function(nIter, nWeights = 6, seed = 123, nCores = 4) {
  
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
      fitness <- getFitness(weights) # estimate fitness in 100 datasets
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
datasetsList <- readRDS(file = "AssessmentFiles/GA_inputData.rds") # IMPORTANTE QUE SE CHAME datasetsList!!!
# Applying train-test here must be necessary

# Execute
library(GA)
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

summary(gaModel)
softmaxWeights(gaModel@solution[1,])
plot(gaModel)
saveRDS(randomResult, file = "AssessmentFiles/GA_results_GA.rds") # 500 iteracións

#...........................................................................####
# Random search ####
tictoc::tic()
randomResult <- randomSearchWeights(
  nIter = 50000, # 1000 iterations = 49.5s // 50000 iteracións 1609 s
  nWeights = 6,
  seed = 1409,
  nCores = 8
)
tictoc::toc()

randomResult$bestWeights
randomResult$bestFitness

saveRDS(randomResult, file = "AssessmentFiles/GA_results_RS.rds")  # 50000 iteracións



