
# Executions 
# Julia G Curras - 22/01/2026

library(MASS)
library(future.apply)
library(future)
library(dplyr)
source(file = "scoreFunction.R", encoding = "UTF-8")
source(file = "simulationFunction.R", encoding = "UTF-8")
source(file = "inputs.R")

plan(multisession)


#.............................................................................
# Item 3  ####
#.............................................................................

tictoc::tic()
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
    valores <- length(item3_sigma_lo[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      prop_de = 0.1,
      sigma_lo = item3_sigma_lo[[params$especifico]][[val]],
      sigma_hi = item3_sigma_hi[[params$especifico]][[val]], 
      sample_sd_strength = item3_sample_sd_strength[[params$especifico]][[val]],
      sample_sd_rho = item3_sample_sd_rho[[params$especifico]][[val]],
      sample_sd_cap = item3_sample_sd_cap[[params$especifico]][[val]])
    )
  } 
  names(resByItem) <- names(goldStandard)
  # cat("\n\tSimulación ", params$especifico, "\n")
  # print(getResultsByItem(resByItem, "item3"))
  
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
  cat("Correlation: ", corResItem3, "\n")
  
  #--- Return ----
  return(corResItem3)
}, future.seed=TRUE)
tictoc::toc()

grid$item3 <- allSim
saveRDS(object = grid, file = "LUSTRE/normScore/Simulations/item3_result.rds")
# saveRDS(object = grid, file = "../Simulations/grid_items3.rds")

