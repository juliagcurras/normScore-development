
##############################################################################- 

#######               SIMULATIONS: PREPARING INPUTS                 ########## 

##############################################################################- 

# Julia G Curras - 22/01/2026

# Preparar parámetros simulacións
library(dplyr)
setwd("C:/Users/julia/Documents/GitHub/normScore/R")


#.............................................................................
# Prametros globales ####
#.............................................................................

## Variables ####
goldStandard <- 1:9
names(goldStandard) <- paste0("Simulation_", 1:9)
especificos <- 0:8
n_proteins <- c(100, 1000, 5000, 10000)
n_per_group <- c(5, 10, 20, 50)
semilla <- c(2012, 265, 301, 9396, 107, 11792, 2511, 242, 1165, 301) # X10
             # 149, 182, 702, 261, 233, 211, 814, 662, 126, 888)

## All combinations of settings ####
grid <- expand.grid(
  n_proteins = n_proteins,
  n_per_group = n_per_group,
  especificos = especificos,
  semilla = semilla
)
dim(grid)

# ordenar
vars_grupo <- setdiff(names(grid), "especificos")  
grid <- grid %>% dplyr::arrange(across(all_of(vars_grupo)), especificos)
grid <- grid[1:36,]

# saveRDS(object = grid, file = "../Simulations/grid.rds")




#.............................................................................
# Parámetros específicos ####
#.............................................................................

## Items 0, 1, 5, 6 ####

### sample_shift_sd ####
item0156_sample_shift_sd <- list(
  seq(0, 0.05, 0.00625),
  seq(0, 0.1, 0.0125),
  seq(0, 0.5, 0.0625), 
  seq(0, 0.5, 0.0625), 
  seq(0.5, 1.5, 0.125),
  seq(0.5, 1.5, 0.125),
  seq(1.5, 3.5, 0.25),
  seq(1.5, 3.5, 0.25))
names(item0156_sample_shift_sd) <- paste0("Sim", 1:length(item0156_sample_shift_sd))

### sample_shift_cap ####
offset <- c(0.05, 0.05, rep(c(0.05, -0.05), (length(item0156_sample_shift_sd)/2)-1))
item0156_sample_shift_cap <- sapply(1:length(item0156_sample_shift_sd), function(x){
  item0156_sample_shift_sd[[x]] + offset[x]
}, simplify = F, USE.NAMES = T)
names(item0156_sample_shift_cap) <- paste0("Sim", 1:length(item0156_sample_shift_cap))




## Item 2####

### rho_between ####
item2_values_rho_between <- list(
  seq(0, 0.05, 0.00625),
  seq(0, 1.25, 0.15625), 
  seq(0, 1, 0.125), 
  seq(0, 0.75, 0.09375), 
  seq(0, 0.5, 0.0625),
  seq(0, 0.5, 0.0625),
  seq(0, 0.25, 0.03125), 
  seq(0, 0.05, 0.00625)) 
names(item2_values_rho_between) <- paste0("Sim", 1:length(item2_values_rho_between))

### rho_within ####
item2_values_rho_within <- list(
  seq(3.5, 1.5, -0.25),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 0.75, 0.09375)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)))
names(item2_values_rho_within) <- paste0("Sim", 1:length(item2_values_rho_within))






## Item 3 ####

### sigma_lo ####
item3_sigma_lo <- list(
  rep(0.4, 9), # Control (-1)
  rep(0.4, 9), 
  rep(0.4, 9)*1.5, 
  rep(0.4, 9), 
  rep(0.4, 9), 
  rev(seq(0, 1.5, 0.1875)), 
  rev(seq(0.4, 2, 0.2)), 
  c(3, 2, 1, 0.5, 0.3, 0.3, 0.3, 0.3, 0.1) # solo forma
)

### sigma_hi ####
item3_sigma_hi <- list(
  rep(0.05, 9), # Control
  rep(0.05, 9), 
  rep(0.05, 9)*1.5, 
  rep(0.05, 9), 
  rep(0.05, 9), 
  seq(0, 1.5, 0.1875)*(1.5),
  seq(0.4, 2, 0.2)*(1.5),
  c(0.05, 0.1, 0.2, 0.4, 1.2, 1.2, 1.2, 1.9, 2.5)
)

### sample_sd_cap ####
item3_sample_sd_cap <- list(
  seq(2.5, 1, -0.1875), # Control - orden inverso, ten que dar 0 
  seq(1, 2.5, 0.1875), # Control - orden inverso, ten que dar 0 
  seq(1, 2.5, 0.1875),  
  seq(0.5, 2, 0.1875),
  seq(0, 1, 0.125),
  seq(0, 1.5, 0.1875), 
  seq(0, 1.5, 0.1875), 
  rep(0.2, 9)
  # 8 
)

### sample_sd_strength ####
item3_sample_sd_strength = list(
  rep(3, 9), 
  rep(3, 9), 
  rep(3, 9), 
  rep(3, 9), 
  rep(3, 9), 
  rep(3, 9), 
  rep(2.5, 9), 
  rep(0.8, 9)
)

### sample_sd_rho ####
item3_sample_sd_rho = list(
  rep(0, 9), 
  rep(0, 9), 
  rep(0, 9), 
  rep(0, 9), 
  rep(0, 9), 
  rep(0, 9), 
  rep(0.2, 9), 
  rep(0.35, 9)
)







## Item 4 ####

### sample_sd_cap ####
item4_sample_sd_cap <- list(
  seq(0, 0.5, 0.0625),
  seq(0, 1, 0.125),
  seq(0, 1.5, 0.1875),
  seq(0.5, 2, 0.1875),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9)
)
names(item4_sample_sd_cap) <- paste0("Sim", 1:length(item4_sample_sd_cap))

### sample_sd_cstrength ####
item4_sample_sd_strength <- list( # rep(c(1, 2), each = 4) 
  rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  seq(0, 0.25, 0.03125)*(-1),
  seq(0, 0.5, 0.0625)*(-1),
  seq(0, 1, 0.125)*(-1),
  seq(0, 1.5, 0.1875)*(-1)
) # 2
names(item4_sample_sd_strength) <- paste0("Sim", 1:length(item4_sample_sd_strength))


### sample_sd_rho ####
item4_sample_sd_rho <- rep(c(1.5, 0.9), each = 4) #0











