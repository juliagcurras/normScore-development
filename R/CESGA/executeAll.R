
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
# All Score  ####
#.............................................................................


allSim <- sapply(1:nrow(grid), function(i){
  
  #--- Simulations ---
  params <- grid[i, ]
  if (params$especificos == 0){ # control negativo, no hay diferencias entre muestras => kendall = 0
    resByItem <- lapply(1:9, function(val) simulate_proteomics_clean(
      semilla = params$semilla+val, 
      n_proteins = params$n_proteins,
      n_per_group = params$n_per_group)
    )
  } else { # ejecuciones en sí, de comprobación
    valores <- length(sevList[[params$especifico]])
    resByItem <- lapply(1:valores, function(val) {
      sev <- sevList[[params$especifico]][[val]]
      simulate_proteomics_clean(
        semilla = params$semilla+val,
        n_proteins = params$n_proteins,
        n_per_group = params$n_per_group,
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
    }
    )
  }
  names(resByItem) <- names(goldStandard)
  
  #--- Cheking score ---
  dm <- as.data.frame(resByItem[[1]][["metadata"]])
  lista <- sapply(resByItem, "[[", 1, simplify = F, USE.NAMES = T)
  refGroup <- "G1"
  altGroup <- "G2"
  
  finalRank <- normScore(
    lista, 
    dm,  
    refGroup = refGroup, 
    altGroup = altGroup 
    # onlyFinalRank = T
  )$finalRanking
  
  outRank <- as.numeric(gsub(pattern = "Simulation_", replacement = "", x = names(finalRank)))
  names(outRank) <- names(finalRank)
  corFinal <- cor(x = goldStandard, y = outRank, method = "kendall")
  
  #--- Return ----
  return(corFinal)
}, simplify = T, USE.NAMES = T)

grid$normScore <- allSim
saveRDS(object = grid, file = paste0(dirOut, "normScore_result.rds"))