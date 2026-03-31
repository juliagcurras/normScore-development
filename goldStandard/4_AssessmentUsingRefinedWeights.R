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
# Este input crease en 4_Refinement_GA_RS
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
resNS1_train <- lapply(datasetsListTrain, normScoreAssembly, itemWeights = pesos)
resNS1_test <- lapply(datasetsListTest, normScoreAssembly, itemWeights = pesos)
resNS1_all <- lapply(allDatasetsList, normScoreAssembly, itemWeights = pesos)

out1 <- list(
  train  = resNS1_train,
  test  = resNS1_test, 
  all = resNS1_all
)

saveRDS(out1, file = "AssessmentFiles/Refinement_AssessingWeigths_fitness1.RDS")


#..........................................................................####
# Fitness 2 ####
fitness <- 2 
resGA <- readRDS(file = paste0(outDir, "Refinement_results_GA_FITNESS", 
                               fitness, "_TRAIN.rds"))
pesos <- softmaxWeights(resGA$model@solution[1,])
resNS1_train <- lapply(datasetsListTrain, normScoreAssembly, itemWeights = pesos)
resNS1_test <- lapply(datasetsListTest, normScoreAssembly, itemWeights = pesos)
resNS1_all <- lapply(allDatasetsList, normScoreAssembly, itemWeights = pesos)

out1 <- list(
  train  = resNS1_train,
  test  = resNS1_test, 
  all = resNS1_all
)

saveRDS(out1, file = "AssessmentFiles/Refinement_AssessingWeigths_fitness2.RDS")



#..........................................................................####
# Fitness 3 ####
fitness <- 3 
resGA <- readRDS(file = paste0(outDir, "Refinement_results_GA_FITNESS", 
                               fitness, "_TRAIN.rds"))
pesos <- softmaxWeights(resGA$model@solution[1,])
resNS1_train <- lapply(datasetsListTrain, normScoreAssembly, itemWeights = pesos)
resNS1_test <- lapply(datasetsListTest, normScoreAssembly, itemWeights = pesos)
resNS1_all <- lapply(allDatasetsList, normScoreAssembly, itemWeights = pesos)

out1 <- list(
  train  = resNS1_train,
  test  = resNS1_test, 
  all = resNS1_all
)

saveRDS(out1, file = "AssessmentFiles/Refinement_AssessingWeigths_fitness3.RDS")

#..........................................................................####
# Fitness 3 ####
fitness <- 4
resGA <- readRDS(file = paste0(outDir, "Refinement_results_GA_FITNESS", 
                               fitness, "_ALL.rds"))
pesos <- softmaxWeights(resGA@solution[1,])
resNS1_all <- lapply(allDatasetsList, normScoreAssembly, itemWeights = pesos)

out1 <- list(
  all = resNS1_all
)

saveRDS(out1, file = "AssessmentFiles/Refinement_AssessingWeigths_fitness4.RDS")












