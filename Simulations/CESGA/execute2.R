
# Executions 
# Julia G Curras - 22/01/2026

library(MASS)
library(future.apply)
library(future)
library(dplyr)
source(file = "/home/ulc/es/jgc/Simulacions/scoreFunction.R", encoding = "UTF-8")
source(file = "/home/ulc/es/jgc/Simulacions/simulationFunction.R", encoding = "UTF-8")
source(file = "/home/ulc/es/jgc/Simulacions/inputs.R", encoding = "UTF-8")
# plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))
dirOut <- "/mnt/lustre/scratch/nlsas/home/ulc/es/jgc/Simulations/"
plan(multisession)

#.............................................................................
# Item 2  ####
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
    valores <- length(item2_values_rho_between[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
      semilla = params$semilla, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group,
      rho_between = item2_values_rho_between[[params$especifico]][[val]],
      rho_within = item2_values_rho_within[[params$especifico]][[val]])
    )
  }
  names(resByItem) <- names(goldStandard)
  
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
  item2All <- item2
  item2 <- sort(item2All) # Ascending
  item2 <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(item2)))
  names(item2) <- names(item2All)
  corResItem2 <- cor(x = goldStandard, y = item2, method = "kendall")
  
  
  #--- Return ----
  names(corResItem2) <- c("Item2")
  return(corResItem2)
}, future.seed=TRUE)

grid$item2 <- allSim
saveRDS(object = grid, file = paste0(dirOut, "item2_result_n16.rds"))


