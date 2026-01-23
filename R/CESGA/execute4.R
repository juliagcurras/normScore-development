
# Executions 
# Julia G Curras - 22/01/2026

library(MASS)
library(future.apply)
library(future)
library(dplyr)
source(file = "scoreFunction.R", encoding = "UTF-8")
source(file = "simulationFunction.R", encoding = "UTF-8")
source(file = "inputs.R")
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))
dirOut <- "/mnt/lustre/scratch/nlsas/home/ulc/es/jgc"
plan(multisession)

#.............................................................................
# Item 4  ####
#.............................................................................

allSim <- future_sapply(1:nrow(grid), function(i) {
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else { # ejecuciones en sí, de comprobación
    valores <- length(item4_sample_sd_cap[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_shift_sd = 0.5,
      sample_sd_cap = item4_sample_sd_cap[[params$especifico]][[val]],
      sample_sd_strength = item4_sample_sd_strength[[params$especifico]][[val]],
      sample_sd_rho = item4_sample_sd_rho[params$especifico])
    )
  }
  names(resByItem) <- names(goldStandard)
  # cat("\n\tSimulación ", params$especifico, "\n")
  # print(getResultsByItem(resByItem, "item4"))
  
  
  #--- Cheking item 4 ---
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  item4All <- sapply(lista, meanSDdiffArea)
  item4 <- sort(item4All) # Ascending
  item4 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item4)))
  names(item4) <- names(item4All)
  corResItem4 <- cor(x = goldStandard, y = item4, method = "kendall")
  
  #--- Return ----
  # return(c(item4All, " - ", corResItem4))
  return(corResItem4)
}, future.seed=TRUE)

grid$item4 <- allSim
saveRDS(object = grid, file = "LUSTRE/normScore/Simulations/item4_result.rds")






