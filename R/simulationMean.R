
# Julia G Currás - 2025/10/06
# Computational benchmarking datasets to test the score

rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")

# Functions
# install.packages("future.apply")
install.packages("progressr")
library(factoextra)
library(dplyr)
library(future.apply)
library(progressr)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")

# Global variables
pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/"
norm <- "mean"



# MEAN ####
## Function ####
simulate_mean_simple <- function(n_proteins = 1000, n_samples = 20, k_groups = 2) {
  group <- rep(1:k_groups, length.out = n_samples)
  base <- matrix(rnorm(n_proteins * n_samples, mean = 20, sd = 2),
                 nrow = n_proteins)
  offsets <- rnorm(n_samples, 0, 3)
  mat <- as.data.frame(sweep(base, 2, offsets, "+"))
  colnames(mat) <- paste0("Sample_", 1:n_samples)
  rownames(mat) <- paste0("Protein_", 1:n_proteins)
  # mat <- mat %>% dplyr::select(Proteins, everything())
  meta <- data.frame(
    Samples = colnames(mat),
    Groups = factor(group)
  )
  matRaw <- as.data.frame(2^mat)
  
  return(list(
    rawData = matRaw,
    logData = mat,
    metadata = meta
  ))
}

simular_proteomica <- function(n_proteinas = 500, 
                               n_muestras = 10, 
                               rango_valores = c(15, 25), 
                               sd_base = 1, 
                               sd_desplazamiento = 0.5, 
                               semilla = 42) {
  set.seed(semilla)
  
  # 1. Generamos matriz base de proteínas (normal)
  base <- matrix(rnorm(n_proteinas * n_muestras, 
                       mean = mean(rango_valores), 
                       sd = sd_base), 
                 nrow = n_proteinas, ncol = n_muestras)
  
  # 2. Introducimos efectos por muestra (desplazamientos)
  desplazamientos <- rnorm(n_muestras, mean = 0, sd = sd_desplazamiento)
  data <- sweep(base, 2, desplazamientos, FUN = "+")
  
  # 3. Escalamos para mantener el rango aproximado
  min_val <- min(data)
  max_val <- max(data)
  data <- (data - min_val) / (max_val - min_val) * (rango_valores[2] - rango_valores[1]) + rango_valores[1]
  
  # 4. Convertimos en data frame y agregamos nombres
  proteinas <- paste0("Protein_", 1:n_proteinas)
  muestras <- paste0("Sample_", 1:n_muestras)
  df <- as.data.frame(data)
  colnames(df) <- muestras
  rownames(df) <- proteinas
  matRaw <- 2^df
  
  meta <- data.frame(
    Samples = colnames(df),
    Groups = factor(group)
  )
  
  # Retornamos la matriz
  return(list(
    rawData = matRaw,
    logData = df,
    metadata = meta
  ))
}


## Individual trial ####
new <- simulate_mean_simple()
new <- simular_proteomica()
newTab <- new$rawData
newTab$Proteins <- rownames(newTab)
newTab <- newTab %>% dplyr::select(Proteins, everything())
newDesign <- new$metadata
writexl::write_xlsx(newTab, path = paste0(pathToData, norm, "/meanProba.xlsx"), col_names = T)
writexl::write_xlsx(newDesign, path = paste0(pathToData, norm, "/meanProbaDesign.xlsx"), col_names = T)

write.table(x = newTab, file = paste0(pathToData, norm, "/meanProba.txt"),
            sep = "\t", row.names = F, col.names = T)
write.table(x = newDesign, file = paste0(pathToData, norm, "/meanProbaDesign.txt"),
            sep = "\t", row.names = F, col.names = T)



## Multiple simulations ####
# Settings
nProts <- c(1000, 5000, 10000)
sSamples <- c(20, 60, 90)
k_groups <- c(2, 4, 6)
bio_var <- c(0.5, 1)
tech_var <- c(1, 2, 3)
noise_sd <- c(0.05, 1, 2)

# All combinations
grid <- expand.grid(
  n_proteins = nProts,
  n_samples = sSamples,
  k_groups = k_groups,
  bio_var = bio_var,
  tech_var = tech_var,
  noise_sd = noise_sd
)
nrow(grid)  # 486 combinacions
saveRDS(grid, file = paste0(pathToData, norm, "/optionsSim.rds"))

# Execution of simulations
set.seed(9396)
handlers(global = TRUE)
handlers("txtprogressbar")  
plan(multisession, workers = parallel::detectCores() - 16)
ordenNorm <- c("CyclicLoess", "GI", "Log","MAD", "Mean", "Median", "Quantile", "RLR", "VSN")

with_progress({
  p <- progressor(steps = 100)
  datasets <- future_lapply(seq_len(nrow(grid)), function(i) {
  # datasets <- future_lapply(1:50, function(i) {
  # datasets <- sapply(1:5, function(i) {
    # Simulate data
    params <- grid[i, ]
    sim <- simulate_mean(
      n_proteins = params$n_proteins,
      n_samples = params$n_samples,
      k_groups = params$k_groups,
      bio_var = params$bio_var,
      tech_var = params$tech_var,
      noise_sd = params$noise_sd,
      seed = 1000 + i
    )
    
    # Retrieving data
    dfGrupos <- sim$metadata
    logIntensityMatrix <- sim$logData
    intensityMatrix <- sim$rawData
    grupos <- levels(as.factor(dfGrupos$Groups))
    
    # Normalization...
    mydata <- list(Log = logIntensityMatrix,
                   Mean = meanNorm(rawMatrix = intensityMatrix), 
                   Median = medianNorm(rawMatrix = intensityMatrix), 
                   GI = GINorm(rawMatrix = intensityMatrix),
                   Quantile = quantileNorm(log2Matrix = logIntensityMatrix),
                   VSN = VSNNorm(rawMatrix = intensityMatrix),
                   CyclicLoess = cyclicLoessNorm(log2Matrix = logIntensityMatrix),
                   RLR = RLRNorm(log2Matrix = logIntensityMatrix),
                   MAD = MADNormalization(log2Matrix = logIntensityMatrix))
    
    # Assessing normalization
    finalRank <- normScore( 
      designMatrix = dfGrupos, 
      normMatrixList = mydata, 
      dfRaw = intensityMatrix)$finalRanking
    
    # finalRank <- finalRank[ordenNorm]
    p()  # avanza la barra
    return(finalRank)
  })
})

# saveRDS(datasets, file = paste0(pathToData, norm, "/results_sim486.RDS"))

# Results into df
dfRes <- sapply(datasets, function(i) i[ordenNorm])
# dfRes <- as.data.frame(dplyr::bind_cols(datasets))
colnames(dfRes) <- paste0("Sim ", 1:ncol(dfRes))
rownames(dfRes) <- ordenNorm

# Top 1 norm: mean is the second one
topNorm <- sapply(datasets, function(i) names(sort(i))[1])
table(topNorm)
# Top 3 norm: mean is the second one
top3Norm <- sapply(datasets, function(i) "Mean" %in% names(sort(i))[1:3])
table(top3Norm)

# Top 5 norm: mean is the second one
top5Norm <- sapply(datasets, function(i) "Mean" %in% names(sort(i))[1:5])
table(top5Norm)

# 






















