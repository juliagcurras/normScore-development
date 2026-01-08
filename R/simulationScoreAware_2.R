
##############################################################################- 

#######                         SIMULACIONS                         ########## 

##############################################################################- 

# Julia G Curras - 08/01/2026

# Simulacion de datos de proteomica (log2-scale)
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")

library(MASS)
library(future.apply)
library(progressr)
library(dplyr)
library(tictoc)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/others/"
norm <- "log"


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

simulate_proteomics_scoreaware <- function(
    # ---- Tamanos 
    n_proteins  = 1000,
    n_per_group = 20,
    k_groups    = 2,
    
    # ---- Mezcla de medias proteicas (log2)
    p_high       = 0.02,
    mu_low_mean  = 18,
    mu_low_sd    = 1.2,
    mu_high_mean = 23,
    mu_high_sd   = 0.6,
    # ---- Estructura de correlacion entre muestras
    rho_samples  = 0.6,    # correlacion base
    sigma_sample = 1.0,    # escala del efecto por muestra
    loading_sd   = 0.6,    # heterogeneidad entre proteinas frente al efecto de muestra
    sigma_resid  = 0.6,    # ruido residual prot~muestra (log2)
    # ---- DE
    prop_de      = 0.05,
    logFC_mean   = 2.0,
    logFC_sd     = 0.6,
    
    # ---- switches de errores (V1)
    add_additive_shift = FALSE,   # "loading" por muestra (global shift en log2)
    ti_cv_target       = 0.10,    # objetivo aproximado de CV en total intensities (raw)
    shift_sd_log2      = NULL,    # si NULL se calcula desde ti_cv_target
    
    add_scale_variance = FALSE,   # sd por muestra ligado al loading
    extra_noise_base   = 0.20,    # sd base del ruido tecnico extra (log2)
    hetero_strength    = 0.9,     # cuanto crece sd con el loading (MeanSD item)
    
    add_intensity_bias = FALSE,   # sesgo dependiente de intensidad (no-lineal) + group-specific
    bias_group_shift   = 0.18,    # diferencia media (en pendiente) entre grupos -> MAplot item
    bias_intercept_sd  = 0.08,
    bias_slope_sd      = 0.08,
    bias_quad_sd       = 0.03,    # componente no-lineal -> reduce Spearman
    bias_share_shape   = FALSE,   # TRUE: misma forma dentro de cada grupo (menos dano a corr)
    
    add_shape_mixture  = FALSE,   # outliers/glitches
    outlier_prop_prot  = 0.01,    # % proteinas afectadas
    outlier_prop_samp  = 0.25,    # % muestras donde se manifiesta el glitch
    outlier_sd_log2    = 1.3,     # magnitud del outlier (log2)
    
    add_missing        = TRUE,    # MNAR + componente por muestra
    target_missing     = 0.15,
    k_mnar             = 1.2,
    missing_by_sample_sd = 0.35,
    
    seed = NULL
){
  
  if (!is.null(seed)) set.seed(seed)
  
  semilla <- ifelse(is.null(semilla), 9396, semilla+9396)
  set.seed(semilla)
  
  # ----------------- 1) Inicio simulacion base -----------------
  groups <- paste0("G", rep(1:k_groups, length.out = n_per_group * k_groups))
  m <- length(groups)
  
  # medias de proteina (log2)
  mu <- rnorm(n_proteins, mean = mu_low_mean, sd = mu_low_sd)
  high_idx <- if (p_high > 0) sample(1:n_proteins, size = round(n_proteins * p_high)) else integer(0)
  if (length(high_idx) > 0) {
    mu[high_idx] <- rnorm(length(high_idx), mean = mu_high_mean, sd = mu_high_sd)
  }
  mu <- pmin(pmax(mu, 15), 25)
  
  # correlacion entre muestras + cargas por proteina
  Sigma <- matrix(rho_samples, nrow = m, ncol = m); diag(Sigma) <- 1
  sample_effects <- as.numeric(MASS::mvrnorm(n = 1, mu = rep(0, m), Sigma = Sigma)) * sigma_sample
  loadings <- rnorm(n_proteins, mean = 0, sd = loading_sd)
  
  # DE
  candidate_idx <- setdiff(1:n_proteins, high_idx)
  n_de <- round(n_proteins * prop_de)
  de_idx <- if (n_de > 0) sample(candidate_idx, size = n_de) else integer(0)
  de_logFC <- if (n_de > 0) rnorm(n_de, mean = logFC_mean, sd = logFC_sd) * sample(c(-1,1), n_de, TRUE) else numeric(0)
  group_effect <- rep(0, n_proteins); if (n_de > 0) group_effect[de_idx] <- de_logFC
  
  # matriz base (log2)
  mat <- matrix(NA_real_, nrow = n_proteins, ncol = m)
  for (i in 1:n_proteins) {
    base     <- mu[i]
    samp_eff <- loadings[i] * sample_effects
    grp_eff  <- ifelse(groups == "G1", group_effect[i], 0)
    resid    <- rnorm(m, mean = 0, sd = sigma_resid)
    mat[i, ] <- base + samp_eff + grp_eff + resid
  }
  rownames(mat) <- paste0("P", sprintf("%05d", 1:n_proteins))
  colnames(mat) <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  V1 <- mat
  
  
  #................................
  # 2) Anade errores tecnicos ----
  
  ## --- Error A: additive shift por muestra (loading/injection) ----
  offsets <- rep(0, m)
  if (add_additive_shift) {
    sd_off <- if (is.null(shift_sd_log2)) sd_from_target_cv(ti_cv_target) else shift_sd_log2
    offsets <- rnorm(m, 0, sd_off)
    V1 <- sweep(V1, 2, offsets, "+")
  }
  
  ## --- Error B: heteroscedasticidad por muestra ligada al loading (MeanSD + PCV) ----
  if (add_scale_variance) {
    # sd_j crece (o decrece) con offsets: esto crea tendencia SD~mean entre muestras
    offs_c <- offsets - mean(offsets)
    sd_j <- extra_noise_base * exp(hetero_strength * offs_c)
    noise <- matrix(rnorm(n_proteins * m, 0, rep(sd_j, each = n_proteins)), nrow = n_proteins)
    V1 <- V1 + noise
  }
  
  ## --- Error C: sesgo dependiente de intensidad (no-lineal) + group-specific (MAplot + corr) ----
  if (add_intensity_bias) {
    A <- rowMeans(V0, na.rm = TRUE)                  # gabundancia verdaderah
    A_c <- as.numeric(A - mean(A))                   # centrado
    Q <- A_c^2 - mean(A_c^2)                         # componente no-lineal centrada
    
    # si share_shape=TRUE, compartimos coeficientes por grupo (menos dano a corr intra-grupo)
    if (bias_share_shape) {
      coefs_by_group <- lapply(unique(groups), function(g) {
        mu_slope <- if (g == paste0("G", k_groups)) bias_group_shift else 0
        list(
          a0 = rnorm(1, 0, bias_intercept_sd),
          a1 = rnorm(1, mu_slope, bias_slope_sd),
          a2 = rnorm(1, 0, bias_quad_sd)
        )
      })
      names(coefs_by_group) <- unique(groups)
      
      for (j in seq_len(m)) {
        cc <- coefs_by_group[[groups[j]]]
        V1[, j] <- V1[, j] + cc$a0 + cc$a1 * A_c + cc$a2 * Q
      }
    } else {
      for (j in seq_len(m)) {
        mu_slope <- if (groups[j] == paste0("G", k_groups)) bias_group_shift else 0
        a0 <- rnorm(1, 0, bias_intercept_sd)
        a1 <- rnorm(1, mu_slope, bias_slope_sd)
        a2 <- rnorm(1, 0, bias_quad_sd)
        V1[, j] <- V1[, j] + a0 + a1 * A_c + a2 * Q
      }
    }
  }
  
  ## --- Error D: outliers/glitches (mixture) ----
  if (add_shape_mixture) {
    n_out <- max(1, round(n_proteins * outlier_prop_prot))
    prot_idx <- sample(seq_len(n_proteins), n_out)
    samp_idx <- which(runif(m) < outlier_prop_samp)
    if (length(samp_idx) > 0) {
      V1[prot_idx, samp_idx] <- V1[prot_idx, samp_idx] +
        matrix(rnorm(n_out * length(samp_idx), 0, outlier_sd_log2),
               nrow = n_out, ncol = length(samp_idx))
    }
  }
  
  #................................
  # 3) Missingness MNAR (+ componente por muestra) ----
  miss_info <- NULL
  if (add_missing) {
    miss_info <- apply_missing_mnar_by_sample(
      V1, target_missing = target_missing, k = k_mnar, sample_sd = missing_by_sample_sd
    )
    V1 <- miss_info$M
  }
  
  #................................
  # 4) Salida: log2 + raw (para item0) ----
  rawV1 <- 2^V1
  
  metadata <- data.frame(
    Samples = colnames(V1),
    Groups  = factor(groups, levels = paste0("G", seq_len(k_groups)))
  )
  
  list(
    originalMat = mat,          # V0 log2 (sin sesgo tecnico)
    logData     = V1,          # V1 log2 (con sesgo tecnico + NA)
    rawData     = rawV1,       # intensidades lineales (para item0 en normScore)
    metadata    = metadata,
    de_info     = data.frame(
      ProteinID = rownames(mat)[de_idx],
      logFC_expected = group_effect[de_idx],
      stringsAsFactors = FALSE
    ),
    tech_info   = list(
      offsets_log2 = offsets,
      miss = miss_info
    )
  )
}




#.........................................................................####
# Function  assessment normalization ####
#.........................................................................####

# Evalua una normalizacion comparando la matriz inicial (rawData) con la final (normData)
# - Ambas pueden contener NA.
# - Por defecto trabaja en log2 con pseudoconteo.
# Devuelve metricas por muestra y un resumen.

evaluate_normalization <- function(rawData_init,
                                   normData_final
) {
  
  # --- 0) Alineacion por IDs (filas=proteinas, columnas=muestras) ---
  stopifnot(is.matrix(rawData_init) || is.data.frame(rawData_init))
  stopifnot(is.matrix(normData_final) || is.data.frame(normData_final))
  X0  <- as.matrix(rawData_init)
  Xh  <- as.matrix(normData_final)
  
  # Alinear por nombres si existen
  common_rows <- intersect(rownames(X0), rownames(Xh))
  common_cols <- intersect(colnames(X0), colnames(Xh))
  if (length(common_rows) == 0 || length(common_cols) == 0) {
    stop("No hay interseccion de filas o columnas entre rawData_init y normData_final.")
  }
  X0w <- X0[common_rows, common_cols, drop = FALSE]
  Xhw <- Xh[common_rows, common_cols, drop = FALSE]
  
  
  # --- 1) Metricas por muestra: MAE sin ajuste afin + Spearman ---
  m <- ncol(X0w)
  per_sample <- vector("list", m)
  # Xhat_aligned <- Xhw  # guardaremos la version alineada (ajustada) por muestra
  
  for (j in seq_len(m)) {
    x0  <- X0w[, j]
    xh  <- Xhw[, j]
    ok  <- is.finite(x0) & is.finite(xh)
    
    if (!any(ok)) {
      per_sample[[j]] <- data.frame(
        sample     = colnames(X0w)[j],
        MAE_affine = NA_real_,
        a_hat      = NA_real_,
        b_hat      = NA_real_,
        Spearman   = NA_real_,
        n          = 0
      )
      next
    }
    
    x0o <- x0[ok]
    xho <- xh[ok]
    
    # AJUSTE MAE AFIN #
    # Ajuste afin por minimos cuadrados: x0 ~ a + b * xhat (efectos de escala)
    # v <- stats::var(xho)
    # b <- if (is.finite(v) && v > 0) stats::cov(x0o, xho) / v else 1
    # a <- mean(x0o) - b * mean(xho)
    # 
    # xh_adj <- a + b * xh
    # Xhat_aligned[, j] <- xh_adj
    # 
    # # Calculo de metricas
    # mae <- mean(abs(x0o - (a + b * xho)))
    # sp  <- suppressWarnings(stats::cor(x0o, xho, method = "spearman"))
    # 
    # per_sample[[j]] <- data.frame(
    #   sample     = colnames(X0w)[j],
    #   MAE_affine = mae,
    #   a_hat      = a,
    #   b_hat      = b,
    #   Spearman   = sp,
    #   n          = sum(ok)
    # )
    
    
    # MAE NORMAL #
    per_sample[[j]] <- data.frame(
      sample   = colnames(X0w)[j],
      MAE      = mean(abs(x0o - xho)),
      Spearman = suppressWarnings(stats::cor(x0o, xho, method = "spearman")),
      n        = sum(ok)
    )
  }
  
  per_sample_df <- do.call(rbind, per_sample)
  
  # --- 2) Resumen global ---
  summary_df <- data.frame(
    MAE_mean = mean(per_sample_df$MAE, na.rm = TRUE),
    Spearman_mean   = mean(per_sample_df$Spearman,   na.rm = TRUE),
    samples_used    = sum(per_sample_df$n > 0),
    entries_used    = sum(per_sample_df$n)
  )
  
  # Salida:
  return(list(
    summary      = summary_df         # resumen global
    # per_sample   = per_sample_df      # metricas por muestra
    # X0_work      = X0w,                # matriz "verdad" en escala de trabajo
    # Xhat_aligned = Xhat_aligned        # normData ajustada (afin) en escala de trabajo
  ))
}



#.........................................................................####
# Execution simulations ####
#.........................................................................####

## Settings 
n_proteins <- c(1000, 5000, 10000)
n_per_group <- 20
k_groups <- 2
sigma_resid <- c(0.2, 0.6) # variacion residual
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
results <- sapply(1:nrow(grid), function(i) {
  params <- grid[i, ]
  sim <- simulate_proteomics(
    n_proteins = params$n_proteins,
    n_per_group =  params$n_per_group,
    k_groups = params$k_groups,
    sigma_resid = params$sigma_resid,
    semilla = i, 
    add_additive_shift = F
  )
  return(sim)
}, 
simplify = F
)

