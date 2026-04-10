

# Simulación de datos de proteómica (log2-scale)
library(dplyr)
# library(MASS)


#.............................................................................
# Main function ####
#.............................................................................

simulate_proteomics_clean <- function(
    n_proteins = 10000,
    n_per_group = 20,
    
    # Medias por proteína (log2)
    mu_mean = 18,
    mu_sd   = 1.2,
    mu_clip = c(15, 25),
    
    # Correlación (intra > inter) vía factor correlacionado entre muestras (ITEM2)
    rho_within  = 0.85,
    rho_between = 0.55,
    loading_sd  = 0.25,
    
    # “Cuña” MA: varianza residual depende de abundancia (ITEM3)
    sigma_hi = 0.05,
    sigma_lo = 0.4,
    gamma_sigma = 2.5,
    
    # DE (simétrica entre grupos)
    prop_de = 0.35,
    logFC_sd = 1,
    logFC_mean = 0,
    hetero_logFC = TRUE,
    fc_hi = 0.55,
    fc_lo = 2.5,
    gamma_fc = 7,
    
    # Item 0/1/2/3/5/6 (shift global por muestra, afecta sumas)
    sample_shift_sd = 0,
    sample_shift_cap = 0.2,
    
    # Item 4 (dependencia media–SD por muestra)
    sample_sd_strength = 0,   # 0 = independencia
    sample_sd_rho      = 0.8, # 0..1
    sample_sd_cap      = 0.35,# cap en log-multiplicador
    
    # Missing
    add_missing = TRUE,
    target_missing = 0.001,
    k_mnar = 1.2,
    missing_by_sample_sd = 0.05,
    
    semilla = 9396
){
  set.seed(semilla)
  
  # --------- grupos y nombres ----------
  groups <- rep(c("G1","G2"), each = n_per_group)
  m <- length(groups)
  
  # --------- medias por proteína ----------
  mu <- rnorm(n_proteins, mu_mean, mu_sd)
  mu <- pmin(pmax(mu, mu_clip[1]), mu_clip[2])
  w_low <- (mu_clip[2] - mu) / (mu_clip[2] - mu_clip[1])
  w_low <- pmin(pmax(w_low, 0), 1)
  
  # --------- DE simétrica (+/− logFC/2) ----------
  n_de <- round(n_proteins * prop_de)
  de_idx <- if (n_de > 0) sample.int(n_proteins, n_de) else integer(0)
  logFC <- rep(0, n_proteins)
  if (n_de > 0) {
    sd_i <- rep(logFC_sd, n_de)
    if (hetero_logFC) {
      mult <- fc_hi + (w_low[de_idx]^gamma_fc) * (fc_lo - fc_hi)
      sd_i <- logFC_sd * mult
    }
    logFC[de_idx] <- rnorm(n_de, logFC_mean, sd_i) * sample(c(-1,1), n_de, TRUE)
  }
  gvec <- ifelse(groups == "G1", +0.5, -0.5)
  DE_mat <- outer(logFC, gvec)
  
  # --------- factor correlacionado para correlación entre muestras ----------
  # (change to modify item 2)
  R <- matrix(rho_between, m, m); diag(R) <- 1
  idx1 <- which(groups=="G1"); idx2 <- which(groups=="G2")
  R[idx1, idx1] <- rho_within; diag(R[idx1, idx1]) <- 1
  R[idx2, idx2] <- rho_within; diag(R[idx2, idx2]) <- 1
  
  ev <- eigen(R, symmetric=TRUE, only.values=TRUE)$values
  if (min(ev) <= 1e-8) R <- R + diag(abs(min(ev)) + 1e-6, m)
  
  s <- as.numeric(MASS::mvrnorm(1, mu = rep(0, m), Sigma = R))
  s[idx1] <- s[idx1] - mean(s[idx1])
  s[idx2] <- s[idx2] - mean(s[idx2])
  
  load <- rnorm(n_proteins, 0, loading_sd)
  FACT_mat <- outer(load, s)
  
  # --------- ruido residual heterocedástico (cuña MA) ----------
  sigma_i <- sigma_hi + (w_low^gamma_sigma) * (sigma_lo - sigma_hi)
  EPS <- matrix(rnorm(n_proteins*m), nrow=n_proteins, ncol=m) * sigma_i
  
  # --------- Item0/1/5/6: shift global por muestra ----------
  b_shift <- rep(0, m)
  if (sample_shift_sd > 0) {
    b_shift <- rnorm(m, mean = 0, sd = sample_shift_sd)
    b_shift <- b_shift - mean(b_shift)
    if (!is.null(sample_shift_cap) && is.finite(sample_shift_cap)) {
      b_shift <- pmin(pmax(b_shift, -sample_shift_cap), sample_shift_cap)
      b_shift <- b_shift - mean(b_shift)
    }
  }
  
  # --------- matriz base (log2) ----------
  X <- matrix(mu, nrow=n_proteins, ncol=m) + FACT_mat + DE_mat + EPS
  if (any(b_shift != 0)) X <- sweep(X, 2, b_shift, "+")
  
  # --------- Item4 (EFECTIVO): SD por muestra dependiente de su media ----------
  if (sample_sd_strength != 0) {
    m_mean <- colMeans(X, na.rm = TRUE)
    z_mean <- as.numeric(scale(rank(m_mean, ties.method = "average")))
    if (anyNA(z_mean)) z_mean <- rep(0, m)
    u <- as.numeric(scale(rnorm(m)))
    if (anyNA(u)) u <- rep(0, m)
    rho <- max(0, min(1, sample_sd_rho))
    z_sd <- rho * z_mean + sqrt(1 - rho^2) * u
    log_mult <- sample_sd_strength * z_sd
    if (!is.null(sample_sd_cap) && is.finite(sample_sd_cap)) {
      log_mult <- pmin(pmax(log_mult, -sample_sd_cap), sample_sd_cap)
    }
    sd_scale <- exp(log_mult)
    X_center <- sweep(X, 2, m_mean, "-")
    X <- sweep(X_center, 2, sd_scale, "*")
    X <- sweep(X, 2, m_mean, "+")
  }
  
  # --------- Aesthetics ----------
  rownames(X) <- paste0("P", sprintf("%05d", 1:n_proteins))
  colnames(X) <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  
  # --------- missing MNAR opcional ----------
  miss_info <- NULL
  if (add_missing) {
    b_miss <- rnorm(m, 0, missing_by_sample_sd)
    b_miss <- b_miss - mean(b_miss)
    f <- function(a){
      p <- plogis(a - k_mnar * X + matrix(b_miss, nrow=n_proteins, ncol=m, byrow=TRUE))
      mean(p, na.rm=TRUE) - target_missing
    }
    a_hat <- uniroot(f, interval=c(-50, 50))$root
    P <- plogis(a_hat - k_mnar * X + matrix(b_miss, nrow=n_proteins, ncol=m, byrow=TRUE))
    M <- matrix(runif(n_proteins*m), nrow=n_proteins, ncol=m) < P
    X[M] <- NA_real_
    miss_info <- list(P = P, mask = M, b_miss = b_miss)
  }
  
  list(
    logData  = X,
    rawData  = 2^X,
    metadata = data.frame(Samples=colnames(X),
                          Groups=factor(groups, levels=c("G1","G2")))
  )
}
