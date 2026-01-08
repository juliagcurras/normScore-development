
##############################################################################- 

#######                         SIMULACIÓNS                         ########## 

##############################################################################- 

# Julia G Curras - 08/01/2026

# Simulación de datos de proteómica (log2-scale)
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")

library(MASS)
library(future.apply)
library(progressr)
library(dplyr)
library(tictoc)
# source(file = "supportFunctions.R", encoding = "UTF-8")
# source(file = "scoreFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/others/"



#.........................................................................####
# Auxiliar Functions ####
#.........................................................................####


# Helper: SD en log2 que produce CV objetivo en factores multiplicativos (aprox.)
# Si factor = 2^offset, offset~N(0,sd), entonces CV(factor)=sqrt(exp((ln2^2)*sd^2)-1)
sd_from_target_cv <- function(cv_target) {
  sqrt(log(cv_target^2 + 1)) / log(2)
}

apply_missing_mnar_by_sample <- function(M, target_missing, k, sample_sd) {
  rng <- range(M, na.rm = TRUE)
  
  # x0 global para aproximar la tasa media objetivo
  f <- function(x0) mean(plogis(k * (x0 - M))) - target_missing
  x0_global <- tryCatch(
    uniroot(f, lower = rng[1] - 5, upper = rng[2] + 5)$root,
    error = function(e) median(M, na.rm = TRUE)
  )
  x0_j <- x0_global + rnorm(ncol(M), 0, sample_sd)
  
  X0 <- matrix(rep(x0_j, each = nrow(M)), nrow = nrow(M))
  prob <- plogis(k * (X0 - M))
  mask <- matrix(runif(length(prob)) < prob, nrow = nrow(M))
  M2 <- M
  M2[mask] <- NA_real_
  list(M = M2, x0 = x0_j, prob = prob, mask = mask)
}

#.........................................................................####
# Function  simulation ####
#.........................................................................####


simulate_proteomics_v2 <- function(
    # ---- Tamaños
  n_proteins  = 1000,
  n_per_group = 20,
  k_groups    = 2,
  
  # ---- Mezcla de medias proteicas (log2)
  p_high       = 0,
  mu_low_mean  = 18,
  mu_low_sd    = 1.2,
  mu_high_mean = 20,
  mu_high_sd   = 0.2,
  
  # ---- Correlación entre muestras (bloques)
  rho_samples  = 0.6,   # se usa como "centro"
  rho_boost_intra = 0.12,
  rho_drop_inter  = 0.12,
  
  # ---- Heterocedasticidad por proteína (para forma del MAplot)
  hetero_resid = TRUE,
  sigma_resid  = 0.35,   # escala base del ruido en log2
  sigma_min_mult = 0.5, # multiplicador en alta expresión
  sigma_max_mult = 1.5, # multiplicador en baja expresión
  
  # ---- DE
  prop_de      = 0.05,
  logFC_mean   = 0,
  logFC_sd     = 2,
  hetero_logFC = TRUE,
  logFC_min_mult = 0.50, # DE más pequeñas a alta expresión
  logFC_max_mult = 1.5, # DE más grandes a baja expresión
  
  # ---- “RLE-safe”
  enforce_rle_center = TRUE,
  cap_rle = 1.5,         # asegura |RLE| <= 1.5 (en log2)
  
  # ---- Missing data (por defecto: NO)
  add_missing        = TRUE,
  target_missing     = 0.05,
  k_mnar             = 1.2,
  missing_by_sample_sd = 0.2,
  
  semilla = NULL
){
  
  semilla <- ifelse(is.null(semilla), 9396, semilla + 9396)
  set.seed(semilla)
  
  # ----------------- 1) Diseño -----------------
  groups <- paste0("G", rep(1:k_groups, each = n_per_group))
  m <- length(groups)
  
  # ----------------- 2) Medias por proteína (log2) -----------------
  mu <- rnorm(n_proteins, mean = mu_low_mean, sd = mu_low_sd)
  high_idx <- if (p_high > 0) sample.int(n_proteins, size = round(n_proteins * p_high)) else integer(0)
  if (length(high_idx) > 0) {
    mu[high_idx] <- rnorm(length(high_idx), mean = mu_high_mean, sd = mu_high_sd)
  }
  mu <- pmin(pmax(mu, 15), 25)
  
  # Escalado 0..1 donde 1 = baja expresión, 0 = alta (para cuña MA)
  mu_min <- 15; mu_max <- 25
  w_low <- (mu_max - mu) / (mu_max - mu_min)  # 1 en 15, 0 en 25
  w_low <- pmin(pmax(w_low, 0), 1)
  
  # ----------------- 3) Correlación entre muestras (bloques) -----------------
  rho_within  <- pmin(0.99, rho_samples + rho_boost_intra)
  rho_between <- pmax(-0.20, rho_samples - rho_drop_inter)
  
  R <- matrix(rho_between, nrow = m, ncol = m)
  diag(R) <- 1
  for (g in unique(groups)) {
    idx <- which(groups == g)
    R[idx, idx] <- rho_within
    diag(R[idx, idx]) <- 1
  }
  
  # Asegurar PD (jitter si hiciera falta)
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) {
    R <- R + diag(abs(min(ev)) + 1e-6, m)
  }
  
  # ----------------- 4) Ruido por proteína: mayor si baja expresión -----------------
  if (hetero_resid) {
    sigma_i <- sigma_resid * (sigma_min_mult + w_low * (sigma_max_mult - sigma_min_mult))
  } else {
    sigma_i <- rep(sigma_resid, n_proteins)
  }
  
  # Matriz de desviaciones correladas por muestra (una fila por proteína)
  E <- MASS::mvrnorm(n = n_proteins, mu = rep(0, m), Sigma = R)
  E <- E * sigma_i  # escala fila a fila (por proteína)
  
  # ----------------- 5) Efecto de grupo (DE) con magnitud dependiente de abundancia -----------------
  n_de <- round(n_proteins * prop_de)
  de_idx <- if (n_de > 0) sample.int(n_proteins, size = n_de) else integer(0)
  
  group_effect <- rep(0, n_proteins)
  if (n_de > 0) {
    if (hetero_logFC) {
      fc_sd_i <- logFC_sd * (logFC_min_mult + w_low[de_idx] * (logFC_max_mult - logFC_min_mult))
    } else {
      fc_sd_i <- rep(logFC_sd, n_de)
    }
    de_logFC <- rnorm(n_de, mean = logFC_mean, sd = fc_sd_i) * sample(c(-1, 1), n_de, TRUE)
    group_effect[de_idx] <- de_logFC
  }
  
  # ----------------- 6) Construcción de la matriz log2 -----------------
  mat <- sweep(E, 1, mu, "+")
  is_G1 <- as.numeric(groups == "G1")
  mat <- mat + outer(group_effect, is_G1)
  
  rownames(mat) <- paste0("P", sprintf("%05d", 1:n_proteins))
  colnames(mat) <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  
  # ----------------- 7) Enforzar propiedades RLE -----------------
  # RLE en log2: rle_ij = x_ij - median_i(x_i.)
  enforce_rle <- function(X) {
    prot_med <- apply(X, 1, median)
    rle <- sweep(X, 1, prot_med, "-")
    samp_shift <- apply(rle, 2, median)  # queremos que sea 0
    sweep(X, 2, samp_shift, "-")
  }
  
  if (isTRUE(enforce_rle_center)) {
    mat <- enforce_rle(mat)
  }
  
  if (!is.null(cap_rle) && is.finite(cap_rle)) {
    prot_med <- apply(mat, 1, median)
    upper <- prot_med + cap_rle
    lower <- prot_med - cap_rle
    mat <- pmin(pmax(mat, lower), upper)
    
    # re-centrar tras el cap (por si el clipping movió medianas)
    if (isTRUE(enforce_rle_center)) {
      mat <- enforce_rle(mat)
    }
  }
  
  V1 <- mat
  
  # ----------------- 8) Missingness (si se pide) -----------------
  miss_info <- NULL
  if (isTRUE(add_missing)) {
    miss_info <- apply_missing_mnar_by_sample(
      V1, target_missing = target_missing, k = k_mnar, sample_sd = missing_by_sample_sd
    )
    V1 <- miss_info$M
  }
  
  rawV1 <- 2^V1
  
  metadata <- data.frame(
    Samples = colnames(V1),
    Groups  = factor(groups, levels = paste0("G", seq_len(k_groups)))
  )
  
  list(
    originalMat = mat,     # V0 log2 (sin NA)
    logData     = V1,      # V1 log2 (con NA si add_missing=TRUE)
    rawData     = rawV1,
    metadata    = metadata,
    de_info     = data.frame(
      ProteinID = rownames(mat)[de_idx],
      logFC_expected = group_effect[de_idx],
      stringsAsFactors = FALSE
    ),
    sim_info = list(
      R = R,
      rho_within = rho_within,
      rho_between = rho_between,
      sigma_i = sigma_i
    )
  )
}



#.........................................................................####
# Execution simulations ####
#.........................................................................####

## Settings 
n_proteins <- c(1000, 5000, 10000)
n_per_group <- 20
k_groups <- 2
sigma_resid <- c(0.2, 0.6) # variación residual
seed <- sample(1:10000, 5)

## All combinations of settings
grid <- expand.grid(
  n_proteins = n_proteins,
  n_per_group = n_per_group,
  k_groups = k_groups,
  sigma_resid = sigma_resid, 
  seed = 9396 #seed
)
# grid <- do.call(rbind, replicate(5, grid, simplify = FALSE)) # repetir 3 veces
nrow(grid)  

# Simulate data - no effect
# results <- sapply(1:nrow(grid), function(i) {
#   params <- grid[i, ]
#   sim <- simulate_proteomics_scoreaware(
#     n_proteins = params$n_proteins,
#     n_per_group =  params$n_per_group,
#     k_groups = params$k_groups,
#     sigma_resid = params$sigma_resid,
#     semilla = i, 
#     add_missing = T
#   )
#   return(sim)
# }, 
# simplify = F
# )

results <- list(simulate_proteomics_v2())

i <- 1
datos <- as.data.frame(results[[i]]["rawData"])
colnames(datos) <- gsub(colnames(datos), pattern = "rawData.", replacement = "")
datos$ProteinID <- rownames(datos)
datos <- datos %>% dplyr::select(ProteinID, everything())
writexl::write_xlsx(x = datos, path = "../Simulations/trial_new_error_.xlsx")
dm <- as.data.frame(results[[i]]["metadata"])
colnames(dm) <- c("Samples", "Groups")
writexl::write_xlsx(x = dm, path = "../Simulations/trial_new_error_DESIGN_.xlsx")




