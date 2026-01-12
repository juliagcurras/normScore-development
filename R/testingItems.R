
##############################################################################- 

#######                 TESTING ITEMS                          ########## 

##############################################################################- 

# Julia G Curras - 12/01/2026

# Evaluación de items individuais do score
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")

library(MASS)
library(future)
library(future.apply)
library(progressr)
library(dplyr)
library(ggplot2)
library(tictoc)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")
source(file = "simulationFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/others/"



#.............................................................................
# Items 0, 1, 5, 6 ####
#.............................................................................

## Parámetros ####

### Parámetros específicos ####
goldStandard <- 1:9
names(goldStandard) <- paste0("Simulation_", 1:9)

valores_sample_shift_sd <- list(
  seq(0, 0.05, 0.00625),
  seq(0, 0.1, 0.0125),
  seq(0, 0.5, 0.0625), 
  seq(0, 0.5, 0.0625), 
  seq(0.5, 1.5, 0.125),
  seq(0.5, 1.5, 0.125),
  seq(1.5, 3.5, 0.25),
  seq(1.5, 3.5, 0.25))
names(valores_sample_shift_sd) <- paste0("Sim", 1:length(valores_sample_shift_sd))

offset <- c(0.05, 0.05, rep(c(0.05, -0.05), (length(valores_sample_shift_sd)/2)-1))
# offset <- rep(c(0.05, -0.05), (length(valores_sample_shift_sd)/2))
valores_sample_shift_cap <- sapply(1:length(valores_sample_shift_sd), function(x){
  valores_sample_shift_sd[[x]] + offset[x]
}, simplify = F, USE.NAMES = T)
names(valores_sample_shift_cap) <- paste0("Sim", 1:length(valores_sample_shift_cap))


### Parámetros globales ####
especificos <- 0:length(valores_sample_shift_sd)
n_proteins <- c(100, 1000, 5000, 10000)
n_per_group <- c(5, 10, 20, 50)

### All combinations of settings ####
grid <- expand.grid(
  n_proteins = n_proteins,
  n_per_group = n_per_group,
  especificos = especificos
)
dim(grid)
grid$semilla <- 1:nrow(grid)



## Ejecución ####
plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))

allSim <- future_sapply(1:nrow(grid), function(i) {
# allSim <- sapply(1:nrow(grid), function(i){
  params <- grid[i, ]
  
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    tictoc::tic()
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
    tictoc::toc()
  } else { # ejecuciones en sí, de comprobación
    valores <- length(valores_sample_shift_sd[[params$especifico]])
    tictoc::tic()
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_shift_sd = valores_sample_shift_sd[[params$especifico]][[val]],
      sample_shift_cap = valores_sample_shift_cap[[params$especifico]][[val]])
      )
    tictoc::toc()
  }
  names(resByItem) <- names(goldStandard)
  
  # checking item 0
  item0List <- lapply(resByItem, function(res){
    dfRaw <- res[["rawData"]]
    totalIntensities <- colSums(dfRaw, na.rm = T)
    item0 <- cv(totalIntensities, proportion = T, na.rm = T)
  })
  
  item0List <- sort(unlist(item0List)) # Ascending
  item0 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item0List)))
  names(item0) <- names(item0List)
  corRes <- cor(x = goldStandard, y = item0, method = "kendall")
  
  return(corRes)
}, future.seed=TRUE)

allSim



# ITEM 0 - correction factor ####

item0List <- lapply(resByItem, function(res){
  dfRaw <- res[["rawData"]]
  totalIntensities <- colSums(dfRaw, na.rm = T)
  item0 <- cv(totalIntensities, proportion = T, na.rm = T)
})

item0List <- sort(unlist(item0List)) # Ascending
item0 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item0List)))
names(item0) <- names(item0List)
cor(x = goldStandard, y = item0, method = "kendall")


