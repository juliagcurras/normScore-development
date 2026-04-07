
#.........
# Loading 

source(file = "simulationFunction.R")
source(file = "normScore.R")


#...................
# Simulate data ####

data <- simulate_proteomics_clean()
normList <- Biomics::doNormalization(rawData = data$rawData, logData = data$logData)


# Assess data
normScore(
  normMatrixList = normList, 
  designMatrix = data$metadata,
  dfRaw = data$rawData, 
  onlyFinalRank = T)

res <- normScore(
  normMatrixList = normList, 
  designMatrix = data$metadata,
  dfRaw = data$rawData, 
  onlyFinalRank = F)

res



#............................................
# Simulate data with higher variability ####

data <- simulate_proteomics_clean(sample_shift_cap = 0.55, sample_shift_sd = 0.5)
normList <- Biomics::doNormalization(rawData = data$rawData, logData = data$logData)

# Assess data
normScore(
  normMatrixList = normList, 
  designMatrix = data$metadata,
  dfRaw = data$rawData, 
  onlyFinalRank = T)

res <- normScore(
  normMatrixList = normList, 
  designMatrix = data$metadata,
  dfRaw = data$rawData, 
  onlyFinalRank = F)

res
