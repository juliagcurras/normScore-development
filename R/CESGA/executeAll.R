
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




saveRDS(object = grid, file = "LUSTRE/normScore/Simulations/all_result.rds")