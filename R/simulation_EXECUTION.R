
# Executions 
# Julia G Curras - 22/01/2026
library(MASS)
library(future)
library(future.apply)
library(progressr)
library(dplyr)
library(tictoc)
source(file = "scoreFunction.R", encoding = "UTF-8")
source(file = "simulation_INPUTS.R")
source(file = "simulationFunction.R", encoding = "UTF-8")
plan(multicore, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 8)))



#.............................................................................
# Items 0, 1, 5, 6  ####
#.............................................................................

allSim <- future_sapply(1:nrow(grid), function(i) {
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
  tictoc::toc()
  return(kendallAll)
}, future.seed=TRUE)

# Traspoñer e adxuntar á grid
grid <- cbind(grid, as.data.frame(t(allSim)))
saveRDS(object = grid, file = "../Simulations/grid_items0156.rds")



#.............................................................................
# Item 2  ####
#.............................................................................


allSim <- future_sapply(1:nrow(grid), function(i) {
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
  tictoc::toc()
  return(corResItem2)
}, future.seed=TRUE)

grid$item2 <- allSim
# saveRDS(object = grid, file = "../Simulations/grid_items2.rds")




#.............................................................................
# Item 3  ####
#.............................................................................

allSim <- future_sapply(1:nrow(grid), function(i) {
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
  tictoc::toc()
  return(corResItem3)
}, future.seed=TRUE)

grid$item3 <- allSim
# saveRDS(object = grid, file = "../Simulations/grid_items3.rds")






#.............................................................................
# Item 4  ####
#.............................................................................

allSim <- future_sapply(1:nrow(grid), function(i) {
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
  tictoc::toc()
  # return(c(item4All, " - ", corResItem4))
  return(corResItem4)
}, future.seed=TRUE)

grid$item4 <- allSim
saveRDS(object = grid, file = "../Simulations/grid_ALLITEMS.rds")







