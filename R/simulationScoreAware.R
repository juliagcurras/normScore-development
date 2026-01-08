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
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/others/"
norm <- "log"


# Adaptado por chatGPT á función orixinal pero non me convence a adaptación dos erros. 
# Non me ten sentido mesturalos. 

#.........................................................................####
# SCORE-AWARE ERROR MODELS (drop-in replacement) ####
#.........................................................................####

# Mantiene la interfaz de simulate_proteomics() y su salida


# --- 1) additive shift (loading / injection factor) ---
add_shift_additive <- function(M, offsets = NULL, sd_shift = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(offsets)) offsets <- rnorm(ncol(M), mean = 0, sd = sd_shift)
  out <- sweep(M, 2, offsets, FUN = "+")
  attr(out, "offsets") <- offsets
  out
}

# --- 2) scale/variance: SCORE-AWARE heteroscedasticity by sample ---
# Reinterpreta 'sd_logscale' como magnitud base; usa offsets si existen para inducir tendencia
add_scale_variance <- function(M, scales = NULL, sd_logscale = 0.3, center = c("row","global"), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  center <- match.arg(center)
  
  # offsets si vienen de add_shift_additive
  offsets <- attr(M, "offsets")
  if (is.null(offsets)) offsets <- rep(0, ncol(M))
  offsets_c <- offsets - mean(offsets)
  
  # sd_j: heteroscedasticidad por muestra ligada a offsets (dispara MeanSD item)
  # sd_logscale controla la "fuerza" global del efecto
  # center se mantiene solo por compatibilidad (no lo usamos para multiplicar desviaciones, sino para "modo")
  base_sd <- pmax(1e-6, sd_logscale)
  hetero_strength <- 0.9  # fijo (puedes exponerlo si quieres más adelante)
  sd_j <- base_sd * exp(hetero_strength * offsets_c)
  
  noise <- matrix(rnorm(nrow(M) * ncol(M), mean = 0, sd = rep(sd_j, each = nrow(M))),
                  nrow = nrow(M), ncol = ncol(M))
  out <- M + noise
  
  attr(out, "sd_by_sample") <- sd_j
  out
}

# --- 3) intensity bias: SCORE-AWARE + group-specific to trigger MA-slope ---
# Mantiene firma (n_knots, bias_sd, share_shape)
add_intensity_bias_spline <- function(M, n_knots = 5, bias_sd = 0.5, share_shape = FALSE,
                                      seed = NULL, groups = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Intensidad "verdadera" proxy
  A <- rowMeans(M, na.rm = TRUE)
  A_c <- as.numeric(A - mean(A, na.rm = TRUE))
  Q <- A_c^2 - mean(A_c^2, na.rm = TRUE)  # componente no lineal centrada
  
  # Group-specific: si no pasan groups, cae a comportamiento no específico
  if (is.null(groups)) groups <- rep("G1", ncol(M))
  groups <- as.character(groups)
  levs <- unique(groups)
  ref_g <- levs[length(levs)]  # coincide con refGroup por defecto en normScore
  alt_g <- levs[1]
  
  # Magnitudes: bias_sd controla todo
  # Hacemos que el grupo "alt" tenga pendiente distinta a "ref" -> MA slope != 0
  group_shift <- 0.35 * bias_sd
  intercept_sd <- 0.15 * bias_sd
  slope_sd     <- 0.20 * bias_sd
  quad_sd      <- 0.08 * bias_sd
  
  out <- M
  
  if (share_shape) {
    # misma forma dentro de grupo -> menos daño a corr intra-grupo
    coef_by_group <- lapply(levs, function(g) {
      mu_slope <- if (g == alt_g) group_shift else 0
      list(
        a0 = rnorm(1, 0, intercept_sd),
        a1 = rnorm(1, mu_slope, slope_sd),
        a2 = rnorm(1, 0, quad_sd)
      )
    })
    names(coef_by_group) <- levs
    
    for (j in seq_len(ncol(M))) {
      cc <- coef_by_group[[groups[j]]]
      out[, j] <- out[, j] + cc$a0 + cc$a1 * A_c + cc$a2 * Q
    }
  } else {
    # coeficientes por muestra -> baja más la correlación
    for (j in seq_len(ncol(M))) {
      mu_slope <- if (groups[j] == alt_g) group_shift else 0
      a0 <- rnorm(1, 0, intercept_sd)
      a1 <- rnorm(1, mu_slope, slope_sd)
      a2 <- rnorm(1, 0, quad_sd)
      out[, j] <- out[, j] + a0 + a1 * A_c + a2 * Q
    }
  }
  
  out
}

# Opción linear: también group-specific (mantiene firma original)
add_intensity_bias_linear <- function(M, sd_slope = 0.1, sd_intercept = 0.2, seed = NULL, groups = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  A <- rowMeans(M, na.rm = TRUE)
  A_c <- A - mean(A, na.rm = TRUE)
  
  if (is.null(groups)) groups <- rep("G1", ncol(M))
  groups <- as.character(groups)
  levs <- unique(groups)
  ref_g <- levs[length(levs)]
  alt_g <- levs[1]
  
  # grupo alt con shift en slope
  group_shift <- 0.8 * sd_slope
  
  a1 <- rnorm(ncol(M), 0, sd_slope) + ifelse(groups == alt_g, group_shift, 0)
  a0 <- rnorm(ncol(M), 0, sd_intercept)
  
  out <- M
  for (j in seq_len(ncol(M))) {
    out[, j] <- out[, j] + (a0[j] + a1[j] * A_c)
  }
  out
}

# --- 4) shape mixture: SCORE-AWARE outliers localized (protein x subset samples) ---
add_shape_mixture_fn <- function(M, prop_mix = 0.15, delta_mean = 1.5, sd_mix = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- M
  
  n_out <- max(1, round(nrow(M) * prop_mix))
  prot_idx <- sample(seq_len(nrow(M)), n_out)
  
  # subset de muestras afectadas (fijo 25% para generar degradación de corr + RLE)
  samp_idx <- which(runif(ncol(M)) < 0.25)
  if (length(samp_idx) == 0) samp_idx <- sample(seq_len(ncol(M)), 1)
  
  # outliers: mezcla con desplazamiento y dispersión
  glitch <- matrix(rnorm(n_out * length(samp_idx), mean = delta_mean, sd = sd_mix),
                   nrow = n_out, ncol = length(samp_idx))
  out[prot_idx, samp_idx] <- out[prot_idx, samp_idx] + glitch
  out
}

# --- 5) SAS transform (lo mantenemos tal cual, ya altera shape/colas) ---
add_shape_sas_fn <- function(M, skew = 0.0, tail = 1.0) {
  mu <- mean(M, na.rm = TRUE)
  sdv <- stats::sd(as.vector(M), na.rm = TRUE)
  Z  <- (M - mu) / sdv
  Y  <- sinh((asinh(Z) + skew) * tail)
  Y  <- Y * sdv + mu
  Y
}

# --- 6) mean-variance lognorm (ya es score-relevante: MeanSD) ---
add_mean_variance_lognorm_fn <- function(M, a0 = 0.25, a1 = -0.02, min_sd = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  A <- rowMeans(M, na.rm = TRUE)
  sd_row <- pmax(a0 + a1 * A, min_sd)
  noise <- matrix(rnorm(length(M), 0, rep(sd_row, times = ncol(M))), nrow = nrow(M))
  M + noise
}

# --- 7) mean-variance Poisson-Gamma (ya mete discreción + heteroscedasticidad) ---
add_mean_variance_pg_fn <- function(M, disp = 0.2, pseudocount = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- 2^M
  shape <- 1 / disp
  rate  <- shape / X
  lambda <- matrix(stats::rgamma(length(X), shape = shape, rate = rate), nrow = nrow(M))
  Y <- matrix(stats::rpois(length(lambda), lambda = as.vector(lambda)), nrow = nrow(M))
  logY <- log2(Y + pseudocount)
  colnames(logY) <- colnames(X)
  rownames(logY) <- rownames(X)
  logY
}

# --- 8) MNAR missingness: SCORE-AWARE with per-sample threshold jitter ---
apply_missing_mnar <- function(mat_log, target_missing = 0.15, k = 1.2) {
  # calibración global
  f_mean_diff <- function(x0, M) mean(plogis(k * (x0 - M))) - target_missing
  rng <- range(mat_log, na.rm = TRUE)
  lo <- rng[1] - 5; hi <- rng[2] + 5
  
  x0_star <- tryCatch(
    uniroot(function(z) f_mean_diff(z, mat_log), lower = lo, upper = hi)$root,
    error = function(e) median(mat_log, na.rm = TRUE)
  )
  
  # jitter por muestra (muy importante para TI/RLE/corr)
  sample_sd <- 0.35
  x0_j <- x0_star + rnorm(ncol(mat_log), 0, sample_sd)
  X0 <- matrix(rep(x0_j, each = nrow(mat_log)), nrow = nrow(mat_log))
  
  prob_miss <- plogis(k * (X0 - mat_log))
  mask_na   <- matrix(runif(length(prob_miss)) < prob_miss, nrow = nrow(mat_log))
  mat_out   <- mat_log
  mat_out[mask_na] <- NA_real_
  
  list(mat_log_na = mat_out, prob_miss = prob_miss, mask_na = mask_na, x0 = x0_star, x0_by_sample = x0_j)
}

# =============================================================================
# simulate_proteomics(): misma firma y salida, inyección de errores score-aware
# =============================================================================
simulate_proteomics <-
  function(
    # ---- Tamaños ----
    n_proteins   = 1000,
    n_per_group  = 10,
    k_groups     = 2,
    # ---- Mezcla de medias proteicas (log2) ----
    p_high       = 0.02,
    mu_low_mean  = 18,
    mu_low_sd    = 1.2,
    mu_high_mean = 23,
    mu_high_sd   = 0.6,
    # ---- Estructura de correlación entre muestras ----
    rho_samples  = 0.6,
    sigma_sample = 1.0,
    loading_sd   = 0.6,
    sigma_resid  = 0.6,
    # ---- DE ----
    prop_de      = 0.05,
    logFC_mean   = 2.0,
    logFC_sd     = 0.6,
    # ---- Faltantes MNAR ----
    add_missing    = TRUE,
    target_missing = 0.20,
    k         = 1.2,
    # ---- Switches de errores (se aplican en este orden) ----
    add_additive_shift     = FALSE,
    additive_sd_shift      = 0.5,
    add_scale_variance     = FALSE,
    scale_sd_logscale      = 0.3,
    scale_center           = c("row","global"),
    add_intensity_bias     = FALSE,
    bias_mode              = c("spline","linear"),
    bias_n_knots           = 5,
    bias_sd                = 0.5,
    bias_share_shape       = FALSE,
    bias_sd_slope          = 0.1,
    bias_sd_intercept      = 0.2,
    add_shape_mixture      = FALSE,
    shape_prop_mix         = 0.15,
    shape_delta_mean       = 1.5,
    shape_sd_mix           = 0.8,
    add_shape_sas          = FALSE,
    sas_skew               = 0.0,
    sas_tail               = 1.0,
    add_meanvar_lognorm    = FALSE,
    mv_a0                  = 0.25,
    mv_a1                  = -0.02,
    mv_min_sd              = 0.05,
    add_meanvar_pg         = FALSE,
    pg_disp                = 0.2,
    pg_pseudocount         = 1,
    # ---- Semilla ----
    semilla        = NULL
  ){
    
    semilla <- ifelse(is.null(semilla), 9396, semilla + 9396)
    set.seed(semilla)
    
    # ----------------- inicio simulación base (V0) -----------------
    groups <- paste0("G", rep(1:k_groups, length.out = n_per_group * k_groups))
    m <- length(groups)
    
    # medias proteicas (log2)
    mu <- rnorm(n_proteins, mean = mu_low_mean, sd = mu_low_sd)
    high_idx <- if (p_high > 0) sample(1:n_proteins, size = round(n_proteins * p_high)) else integer(0)
    if (length(high_idx) > 0) {
      mu[high_idx] <- rnorm(length(high_idx), mean = mu_high_mean, sd = mu_high_sd)
    }
    mu <- pmin(pmax(mu, 15), 25)
    
    # correlación entre muestras + cargas por proteína
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
    
    originalMat <- mat  # V0 (sin errores técnicos)
    
    # ----------------- inyección de errores (V1 score-aware) -----------------
    # 1) Aditivo (loading)
    if (add_additive_shift) {
      mat <- add_shift_additive(mat, sd_shift = additive_sd_shift, seed = semilla)
    }
    
    # 2) Varianza por muestra ligada a loading (MeanSD/PCV/corr)
    if (add_scale_variance) {
      mat <- add_scale_variance(mat, sd_logscale = scale_sd_logscale,
                                center = match.arg(scale_center), seed = semilla)
    }
    
    # 3) Sesgo dependiente de intensidad + componente group-specific (MA slope + corr)
    if (add_intensity_bias) {
      mode <- match.arg(bias_mode)
      if (mode == "spline") {
        mat <- add_intensity_bias_spline(
          mat, n_knots = bias_n_knots, bias_sd = bias_sd,
          share_shape = bias_share_shape, seed = semilla, groups = groups
        )
      } else {
        mat <- add_intensity_bias_linear(
          mat, sd_slope = bias_sd_slope, sd_intercept = bias_sd_intercept,
          seed = semilla, groups = groups
        )
      }
    }
    
    # 4) Cambio de forma: outliers localizados (mixture)
    if (add_shape_mixture) {
      mat <- add_shape_mixture_fn(mat, prop_mix = shape_prop_mix,
                                  delta_mean = shape_delta_mean, sd_mix = shape_sd_mix, seed = semilla)
    }
    if (add_shape_sas) {
      mat <- add_shape_sas_fn(mat, skew = sas_skew, tail = sas_tail)
    }
    
    # 5) Relación media–varianza
    if (add_meanvar_lognorm) {
      mat <- add_mean_variance_lognorm_fn(mat, a0 = mv_a0, a1 = mv_a1, min_sd = mv_min_sd, seed = semilla)
    }
    if (add_meanvar_pg) {
      mat <- add_mean_variance_pg_fn(mat, disp = pg_disp, pseudocount = pg_pseudocount, seed = semilla)
    }
    
    # ----------------- faltantes  MNAR (con jitter por muestra) -----------------
    if (add_missing) {
      miss <- apply_missing_mnar(mat, target_missing = target_missing, k = k)
      mat <- miss$mat_log_na
    }
    
    # ----------------- salida -----------------
    meta <- data.frame(
      Samples = colnames(mat),
      Groups  = factor(groups, levels = paste0("G", 1:k_groups))
    )
    matRaw <- 2^mat
    
    de_info <- data.frame(
      ProteinID      = rownames(mat)[de_idx],
      logFC_expected = group_effect[de_idx],
      ref_group      = "G1",
      other_group    = if (k_groups == 2) "G2" else "others",
      stringsAsFactors = FALSE
    )
    
    return(list(
      originalMat = originalMat,
      rawData  = matRaw,
      logData  = as.data.frame(mat),
      metadata = meta,
      de_info  = de_info
    ))
  }
