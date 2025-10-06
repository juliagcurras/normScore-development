
# Julia G Currás - 2025/10/06
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")
# Computational benchmarking datasets to test the score
library(factoextra)
library(dplyr)
pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/"



# 1. LOG ####
simulate_log2 <- function(n_proteins = 1000, n_samples = 40) {
  mu <- runif(n_proteins, 1e3, 1e6)
  data <- sapply(1:n_samples, function(i) rnorm(n_proteins, mu, sd = mu * 0.2))
  return(data)
}
abc <- simulate_log2()


# 2. MEAN ####
simulate_mean <- function(n_proteins = 1000, n_samples = 40) {
  base <- matrix(rnorm(n_proteins * n_samples, mean = 20, sd = 2),
                 nrow = n_proteins)
  offsets <- rnorm(n_samples, 0, 2)
  base <- sweep(base, 2, offsets, "+")
  base <- 2^(base)
  base <- as.data.frame(base)
  colnames(base) <- paste0("Sample ", 1:ncol(base))
  base$Names <- paste0(1:nrow(base), LETTERS)
  base <- base %>% dplyr::select(Names, everything())
  write.table(x = base, file = "meanProba.txt", sep = "\t", row.names = F, col.names = T)
  return(base)
}


new <- simulate_mean()
pcares <- prcomp(t(new))
fviz_pca_ind(pcares, geom.ind = "point", 
             col.ind = "#FC4E07", 
             axes = c(1, 2), 
             pointsize = 1.5) 
newLog <- meanNorm(new)
pcares <- prcomp(t(newLog))
fviz_pca_ind(pcares, geom.ind = "point", 
             col.ind = "#FC4E07", 
             axes = c(1, 2), 
             pointsize = 1.5) 


simulate_mean <- function(
    n_proteins = 1000,
    n_samples = 12,
    k_groups = 3,             # número de grupos biológicos
    bio_var = 1,              # magnitud de variabilidad biológica
    tech_var = 3,             # magnitud del sesgo técnico
    noise_sd = 1,             # ruido residual
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Asignamos las muestras a los grupos biológicos
  group <- rep(1:k_groups, length.out = n_samples)
  
  # Nivel medio global
  mu <- 20
  
  # Efecto biológico (cada grupo tiene una desviación distinta)
  bio_effects <- rnorm(k_groups, 0, bio_var)
  
  # Efecto técnico aditivo por muestra
  tech_effects <- rnorm(n_samples, 0, tech_var)
  
  # Generamos la matriz
  mat <- matrix(NA, nrow = n_proteins, ncol = n_samples)
  for (j in 1:n_samples) {
    mat[, j] <- rnorm(n_proteins, mean = mu + bio_effects[group[j]] + tech_effects[j], 
                      sd = noise_sd)
  }
  mat <- as.data.frame(2^mat)
  
  # Añadimos nombres y metadatos
  colnames(mat) <- paste0("Sample_", 1:n_samples)
  rownames(mat) <- paste0("Protein_", 1:n_proteins)
  # mat <- mat %>% dplyr::select(Proteins, everything())
  meta <- data.frame(
    # sample = colnames(mat)[-1],
    sample = colnames(mat),
    group = factor(group),
    tech_offset = tech_effects
  )
  
  return(list(
    data = mat,
    metadata = meta,
    params = list(
      n_proteins = n_proteins,
      n_samples = n_samples,
      k_groups = k_groups,
      bio_var = bio_var,
      tech_var = tech_var,
      noise_sd = noise_sd
    )
  ))
}

nProts <- c(1000, 5000, 10000)
sSamples <- c(20, 60, 90)
k_groups <- c(2, 4, 6)
bio_var <- c(0.5, 1)
tech_var <- c(1, 2, 3)
noise_sd <- c(0.05, 1, 2)

# Creamos todas las combinaciones
grid <- expand.grid(
  n_proteins = nProts,
  n_samples = sSamples,
  k_groups = k_groups,
  bio_var = bio_var,
  tech_var = tech_var,
  noise_sd = noise_sd
)

nrow(grid)  # 27 combinaciones


# Lista para guardar datasets (o se puede guardar en disco directamente)
datasets <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  params <- grid[i, ]
  
  sim <- simulate_mean(
    n_proteins = params$n_proteins,
    n_samples = params$n_samples,
    k_groups = params$k_groups,
    bio_var = params$bio_var,
    tech_var = params$tech_var,
    noise_sd = params$noise_sd,
    seed = 1000 + i
  )
  
  # Guardamos en la lista (opcional: se puede guardar a disco en lugar de memoria)
  datasets[[i]] <- sim
}

saveRDS(datasets, file = paste0(pathToData, "/sim486.rds"))

# for (i in seq_len(nrow(grid))) {
#   filename <- paste0("sim_dataset_", i,
#                      "_nProts", grid$n_proteins[i],
#                      "_nSamples", grid$n_samples[i],
#                      "_k", grid$k_groups[i],
#                      "_b", grid$bio_var[i],
#                      "_t", grid$tech_var[i],
#                      "_noise", grid$noise_sd[i], ".rds")
#   
#   saveRDS(datasets[[i]], file = filename)
# }





# new <- simulate_mean()
# newTab <- new$data
# newDesign <- new$metadata[,-3]
# write.table(x = newTab, file = paste0(pathToData, "mean/meanProba.txt"), 
#             sep = "\t", row.names = F, col.names = T)
# write.table(x = newDesign, file = paste0(pathToData, "mean/meanProbaDesign.txt"), 
#             sep = "\t", row.names = F, col.names = T)
#



#






















# 3. MEDIAN ####
simulate_median <- function(n_proteins = 1000, n_samples = 10) {
  base <- matrix(rnorm(n_proteins * n_samples, mean = 10, sd = 2),
                 nrow = n_proteins)
  offsets <- rnorm(n_samples, 0, 3)
  base <- sweep(base, 2, offsets, "+")
  # añadimos outliers
  for (i in 1:n_samples) {
    idx <- sample(1:n_proteins, 10)
    base[idx, i] <- base[idx, i] + rnorm(10, 0, 15)
  }
  return(base)
}


# 4. TI ####
simulate_total_intensity <- function(n_proteins = 1000, n_samples = 10) {
  base <- matrix(rnorm(n_proteins * n_samples, mean = 100, sd = 10),
                 nrow = n_proteins)
  scale_factors <- runif(n_samples, 0.5, 1.5)
  base <- sweep(base, 2, scale_factors, "*")
  return(base)
}


# 5. VSN ####
simulate_vsn <- function(n_proteins = 1000, n_samples = 10) {
  mu <- runif(n_proteins, 50, 1e4)
  data <- sapply(1:n_samples, function(i) rnorm(n_proteins, mu, sd = sqrt(mu)))
  return(data)
}

