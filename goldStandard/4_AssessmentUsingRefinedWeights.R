###############################################################################-

############             Assessment using weights             #################-

###############################################################################-


# Julia G Curras - 2026/03/24
rm(list=ls())
graphics.off()
outDir <-  "C:/Users/julia/Documents/GitHub/normScore/goldStandard/AssessmentFiles/"
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")


#..........................................................................####
# Functions ####
source(file = "../R/scoreFunction.R", encoding = "UTF-8")
source(file = "3_AssessmentFunctions.R", encoding = "UTF-8")
softmaxWeights  <- function(weights, fitness = 1) {
  # O GA() xera números reales que nos adaptamos ás restriccions dos pesos: sumar 1 e ser positivas
  weights <- weights - max(weights) # Esto é solo para que non de +Inf con valores moi altos
  weights <- exp(weights)
  weights <- weights / sum(weights)
  return(weights)
}

normScoreAssembly <- function(datasetList, itemWeights){
  matriz <- datasetList$matriz
  item0 <- datasetList$item0
  solution <- datasetList$solution
  solutions <- datasetList$solutions
  
  # Add weights #
  scoreDF_norm <- sweep(matriz, 2, itemWeights, `*`)
  
  # Rank #
  scores_matrix <- t(scoreDF_norm)
  scoreDF_norm$Total <- rowSums(scoreDF_norm)
  scoreDF_norm$TotalCorrected <- scoreDF_norm$Total
  scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "TotalCorrected"] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "TotalCorrected"]*item0
  
  # Sort #
  ordenFilas <- c("Log", "Median", "Mean", "TI", "Quantile", "CyclicLoess", 
                  "RLR", "VSN")
  scoreDF_norm <- scoreDF_norm %>%
    dplyr::arrange(Total, match(rownames(.), ordenFilas)) # en caso de empates, se ordena por sencillez
  
  # Return better
  return(scoreDF_norm)
}



#..........................................................................####
# Datasets ####
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



#..........................................................................####
# Fitness 1 ####
fitness <- 1 
resGA <- readRDS(file = paste0(outDir, "Refinement_results_GA_FITNESS", 
                               fitness, "_TRAIN.rds"))
pesos <- softmaxWeights(resGA$model@solution[1,])
resNS1 <- lapply(datasetsListTrain, normScoreAssembly, itemWeights = pesos)













