
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
# Prametros globales ####
#.............................................................................

## Variables ####
goldStandard <- 1:9
names(goldStandard) <- paste0("Simulation_", 1:9)
especificos <- 0:8
n_proteins <- 1000 #c(100, 1000, 5000, 10000)
n_per_group <- 15 #c(5, 10, 20, 50)

## All combinations of settings ####
grid <- expand.grid(
  n_proteins = n_proteins,
  n_per_group = n_per_group,
  especificos = especificos
)
dim(grid)
grid$semilla <- 1:nrow(grid)


#.............................................................................
# Items 0, 1, 5, 6 ####
#.............................................................................

## Parámetros específicos ####
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


## Ejecución ####
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))

# allSim <- future_sapply(1:nrow(grid), function(i) {
allSim <- sapply(1:nrow(grid), function(i){
    tictoc::tic()
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else { # ejecuciones en sí, de comprobación
    valores <- length(valores_sample_shift_sd[[params$especifico]])
    # tictoc::tic()
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_shift_sd = valores_sample_shift_sd[[params$especifico]][[val]],
      sample_shift_cap = valores_sample_shift_cap[[params$especifico]][[val]])
      )
    # tictoc::toc()
  }
  names(resByItem) <- names(goldStandard)
  
  #--- Checking item 0 ---
  item0List <- lapply(resByItem, function(res){
    dfRaw <- res[["rawData"]]
    totalIntensities <- colSums(dfRaw, na.rm = T)
    item0 <- cv(totalIntensities, proportion = T, na.rm = T)
  })
  
  item0List <- sort(unlist(item0List)) # Ascending
  item0 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item0List)))
  names(item0) <- names(item0List)
  corResItem0 <- cor(x = goldStandard, y = item0, method = "kendall")
  
  
  #--- Cheking item 1 ---
  dm <- as.data.frame(resByItem[[1]][["metadata"]])
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  dfPCV <- data.frame(lapply(lista, getPCV, grupos = unique(dm$Groups), 
                             dfGrupos = dm))
  dfPCV <- as.data.frame(t(dfPCV))
  dfPCV$PCV <- apply(dfPCV, 1, mean, na.rm = T)
  item1All <- stats::setNames(dfPCV$PCV, rownames(dfPCV))
  item1 <- sort(item1All) # Ascending
  item1 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item1)))
  names(item1) <- names(item1All)
  corResItem1 <- cor(x = goldStandard, y = item1, method = "kendall")
  
  
  #--- Cheking item 5 ---
  item5All <- sapply(lista, rleMAPE)
  item5 <- sort(item5All) # Ascending
  item5 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item5)))
  names(item5) <- names(item5All)
  corResItem5 <- cor(x = goldStandard, y = item5, method = "kendall")

  
    #--- Cheking item 6 ---
  item6All <- sapply(lista, tiMAPE)
  item6 <- sort(item6All) # Ascending
  item6 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item6)))
  names(item6) <- names(item6All)
  corResItem6 <- cor(x = goldStandard, y = item6, method = "kendall")
  
  
  #--- Return ----
  kendallAll <- c(corResItem0, corResItem1, corResItem5, corResItem6)
  names(kendallAll) <- c("Item0", "Item1", "Item5", "Item6")
  tictoc::toc()
  return(kendallAll)
}, simplify = T, USE.NAMES = T)
# }, future.seed=TRUE)

allSim # traspoñer e adxuntar á grid
# grid <- cbind(grid, as.data.frame(t(allSim)))






#.............................................................................
# Item 2 ####
#.............................................................................

## Parámetros específicos ####
valores_values_rho_between <- list(
  seq(0, 0.05, 0.00625),
  seq(0, 1.25, 0.15625), 
  seq(0, 1, 0.125), 
  seq(0, 0.75, 0.09375), 
  seq(0, 0.5, 0.0625),
  seq(0, 0.5, 0.0625),
  seq(0, 0.25, 0.03125), 
  seq(0, 0.05, 0.00625)) 
names(valores_values_rho_between) <- paste0("Sim", 1:length(valores_values_rho_between))

valores_values_rho_within <- list(
  seq(3.5, 1.5, -0.25),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 0.75, 0.09375)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)))
names(valores_values_rho_within) <- paste0("Sim", 1:length(valores_values_rho_within))


## Ejecución ####
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))

# allSim <- future_sapply(1:nrow(grid), function(i) {
allSim <- sapply(1:nrow(grid), function(i){
# allSim <- sapply(1:4, function(i){
  tictoc::tic()
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else { # ejecuciones en sí, de comprobación
    valores <- length(valores_values_rho_between[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      rho_between = valores_values_rho_between[[params$especifico]][[val]],
      rho_within = valores_values_rho_within[[params$especifico]][[val]])
      # rho_within = -valores_values_rho_between[[params$especifico]][[val]])
    )
  }
  names(resByItem) <- names(goldStandard)
  print(getResultsByItem(resByItem, "item2"))
  
  
  #--- Cheking item 6 ---
  dm <- as.data.frame(resByItem[[1]][["metadata"]])
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  allVectorsCorr <- lapply(lista, getCorrelationVector,
                           dfGrupos = dm,
                           metodo = "spearman")
  
  dfCor <- data.frame(sapply(allVectorsCorr, "length<-", 
                             max(lengths(allVectorsCorr))))
  
  item2 <- sapply(colnames(dfCor), function(j) {
    i <- dfCor[,j]
    # faigo  1-correlation porque todas as métricas restantes siguen o patrón a 
    # menor número, mellor é a métrica, para que está tamén sexa así. 
    1-(median(i, na.rm = T)-IQR(i, na.rm = T)/3) 
  }, simplify = T, USE.NAMES = T)
  # item2["CyclicLoess"] <- item2["CyclicLoess"]*1.2
  item2All <- item2
  item2 <- sort(item2All) # Ascending
  item2 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item2)))
  names(item2) <- names(item2All)
  corResItem2 <- cor(x = goldStandard, y = item2, method = "kendall")
  
  
  #--- Return ----
  names(corResItem2) <- c("Item2")
  tictoc::toc()
  return(corResItem2)
}, simplify = T, USE.NAMES = T)
# }, future.seed=TRUE)

allSim #



#.............................................................................
# Item 4: meanSD ####
#.............................................................................

## Parámetros específicos ####
valores_sample_sd_cap <- list(
  seq(0, 0.5, 0.0625),
  seq(0, 1, 0.125),
  seq(0, 1.5, 0.1875),
  seq(0.5, 2, 0.1875),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9)
  # seq(0.5, 1.5, 0.125), # problemas ->
  # seq(1, 2, 0.125), 
  # seq(1, 3, 0.25), 
  # seq(1.5, 3, 0.1875)
  )
names(valores_sample_sd_cap) <- paste0("Sim", 1:length(valores_sample_sd_cap))

values_sample_sd_strength <- list( # rep(c(1, 2), each = 4) 
  rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  seq(0, 0.25, 0.03125)*(-1),
  seq(0, 0.5, 0.0625)*(-1),
  seq(0, 1, 0.125)*(-1),
  seq(0, 1.5, 0.1875)*(-1)
  ) # 2
values_sample_sd_rho <- rep(c(1.5, 0.9), each = 4) #0

## Ejecución ####
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))

# allSim <- future_sapply(1:nrow(grid), function(i) {
allSim <- sapply(1:nrow(grid), function(i){
# allSim <- sapply(1, function(i){
  tictoc::tic()
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else { # ejecuciones en sí, de comprobación
    valores <- length(valores_sample_sd_cap[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_shift_sd = 0.5,
      sample_sd_cap = valores_sample_sd_cap[[params$especifico]][[val]],
      sample_sd_strength = values_sample_sd_strength[[params$especifico]][[val]],
      # sample_sd_strength = values_sample_sd_strength[params$especifico],
      sample_sd_rho = values_sample_sd_rho[params$especifico])
    )
  }
  names(resByItem) <- names(goldStandard)
  cat("\n\tSimulación ", params$especifico, "\n")
  print(getResultsByItem(resByItem, "item4"))
  
  
  #--- Cheking item 4 ---
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  item4All <- sapply(lista, meanSDdiffArea)
  item4 <- sort(item4All) # Ascending
  item4 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item4)))
  names(item4) <- names(item4All)
  corResItem4 <- cor(x = goldStandard, y = item4, method = "kendall")
  
  #--- Return ----
  tictoc::toc()
  return(c(item4All, " - ", corResItem4))
}, simplify = F, USE.NAMES = T)
# }, future.seed=TRUE)

allSim #



#.............................................................................
# Item 3: MA plot ####
#.............................................................................

# Cambio forma de nubes / cambio pendiente.

## Parámetros específicos ####
### Opcion 1 ####
valores_sigma_lo <- list(
  seq(0, 1, 0.125), # simulacion non corresponde con gold standard
  seq(0, 0.5, 0.0625), 
  seq(0, 1.5, 0.1875),
  seq(0.5, 2, 0.1875), # simulacion non corresponde con gold standard
  seq(0.5, 1.5, 0.125),
  seq(1, 2, 0.125), 
  seq(1, 3, 0.25), # simulacion non corresponde con gold standard
  seq(1.5, 3, 0.1875)
)
names(valores_sigma_lo) <- paste0("Sim", 1:length(valores_sigma_lo))
values_sigma_hi <- lapply(rev(valores_sigma_lo), `*`, 1.5)

### Opcion 2 ####
opcion2 <- T
valores_sample_sd_cap <- list(
  seq(0, 0.25, 0.03125),
  seq(0, 0.5, 0.0625),
  seq(0, 1, 0.125),
  seq(0, 1.5, 0.1875),
  seq(0.5, 2, 0.1875),
  seq(0.5, 1.5, 0.125),
  seq(1, 2, 0.125), 
  seq(1, 2.5, 0.1875) 
)
names(valores_sample_sd_cap) <- paste0("Sim", 1:length(valores_sample_sd_cap))

## Ejecución ####
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))

# allSim <- future_sapply(1:nrow(grid), function(i) {
allSim <- sapply(1:nrow(grid), function(i){
# allSim <- sapply(1:2, function(i){
  tictoc::tic()
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else if (!opcion2){ # ejecuciones en sí, de comprobación
    valores <- length(valores_sigma_lo[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      prop_de = 0.1,
      sigma_lo = valores_sigma_lo[[params$especifico]][[val]],
      sigma_hi = values_sigma_hi[[params$especifico]][[val]])
    )
  } else {
    valores <- length(valores_sigma_lo[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_sd_strength = 3,
      sample_sd_rho = 0,
      sample_sd_cap = valores_sample_sd_cap[[params$especifico]][[val]])
    )
  }
  names(resByItem) <- names(goldStandard)
  cat("\n\tSimulación ", params$especifico, "\n")
  print(getResultsByItem(resByItem, "item3"))
  
  #--- Cheking item 3 ---
  dm <- as.data.frame(resByItem[[1]][["metadata"]])
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  refGroup <- "G1"
  altGroup <- "G2"
  samplesG1 <- dm[dm$Groups == refGroup, "Samples"]
  samplesG2 <- dm[dm$Groups == altGroup, "Samples"]
  
  item3All <- sapply(lista, maDiffAreas, samplesG1 = samplesG1, 
                  samplesG2 = samplesG2)
  item3 <- sort(item3All) # Ascending
  item3 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item3)))
  names(item3) <- names(item3All)
  corResItem3 <- cor(x = goldStandard, y = item3, method = "kendall")
  
  #--- Return ----
  tictoc::toc()
  return(corResItem3)
}, simplify = T, USE.NAMES = T)
# }, future.seed=TRUE)

allSim # algunhas simulacións non corresponden coa orde de gold standar (cambiar)





#.............................................................................
# All score ####
#.............................................................................
sev <- 0.1 # 0-
results <- simulate_proteomics_clean(
  sample_shift_sd = 0.12 * sev,
  sample_sd_cap = 1*sev,
  sample_sd_strength = 1.2 * sev,
  sample_sd_rho = 0.3 * sev,
  rho_within = 0.85 - 0.25 * sev,
  rho_between = 0.55 - 0.25 * sev,
  loading_sd = 0.25 + 0.20 * sev,
  sigma_hi = 0.05 - 0.02 * sev,
  sigma_lo = 0.40 + 0.40 * sev,
  gamma_sigma = 2.5 + 1.3 * sev,
  target_missing = 0.001 + 0.019 * sev,
  k_mnar = 1.2 + 0.4 * sev,
  missing_by_sample_sd = 0.05 + 0.15 * sev
)
getResults(results)



## Ejecución ####
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))

# allSim <- future_sapply(1:nrow(grid), function(i) {
allSim <- sapply(1:nrow(grid), function(i){
  # allSim <- sapply(1:2, function(i){
  tictoc::tic()
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else if (!opcion2){ # ejecuciones en sí, de comprobación
    valores <- length(valores_sigma_lo[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      prop_de = 0.1,
      sigma_lo = valores_sigma_lo[[params$especifico]][[val]],
      sigma_hi = values_sigma_hi[[params$especifico]][[val]])
    )
  } else {
    valores <- length(valores_sigma_lo[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_sd_strength = 3,
      sample_sd_rho = 0,
      sample_sd_cap = valores_sample_sd_cap[[params$especifico]][[val]])
    )
  }
  names(resByItem) <- names(goldStandard)
  cat("\n\tSimulación ", params$especifico, "\n")
  print(getResultsByItem(resByItem, "item3"))
  
  #--- Cheking item 3 ---
  dm <- as.data.frame(resByItem[[1]][["metadata"]])
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  refGroup <- "G1"
  altGroup <- "G2"
  samplesG1 <- dm[dm$Groups == refGroup, "Samples"]
  samplesG2 <- dm[dm$Groups == altGroup, "Samples"]
  
  item3All <- sapply(lista, maDiffAreas, samplesG1 = samplesG1, 
                     samplesG2 = samplesG2)
  item3 <- sort(item3All) # Ascending
  item3 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item3)))
  names(item3) <- names(item3All)
  corResItem3 <- cor(x = goldStandard, y = item3, method = "kendall")
  
  #--- Return ----
  tictoc::toc()
  return(corResItem3)
}, simplify = T, USE.NAMES = T)
# }, future.seed=TRUE)




