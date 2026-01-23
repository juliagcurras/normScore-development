
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
plan(multisession)


#.............................................................................
# Items 0, 1, 5, 6  ####
#.............................................................................

allSim <- future.apply::future_sapply(1:nrow(grid), function(i) {
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else { # ejecuciones en sí, de comprobación
    valores <- length(item0156_sample_shift_sd[[params$especifico]])
    # tictoc::tic()
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      sample_shift_sd = item0156_sample_shift_sd[[params$especifico]][[val]],
      sample_shift_cap = item0156_sample_shift_cap[[params$especifico]][[val]])
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
  return(kendallAll)
}, future.seed=TRUE)

# Traspoñer e adxuntar á grid
grid <- cbind(grid, as.data.frame(t(allSim)))
saveRDS(object = grid, file = "LUSTRE/normScore/Simulations/item0156_result.rds")

