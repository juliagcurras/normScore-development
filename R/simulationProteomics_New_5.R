
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
library(ggplot2)
library(tictoc)
# source(file = "supportFunctions.R", encoding = "UTF-8")
# source(file = "scoreFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/others/"


#.............................................................................
# Auxiliar functions ####
#.............................................................................

getResults <- function(results){
  # Extract data
  dfRaw <- as.data.frame(results[["rawData"]])
  datos <- as.data.frame(results[["logData"]])
  dm <- as.data.frame(results[["metadata"]])
  grupos <- unique(dm$Groups)
  # colnames(dfRaw) <- gsub(colnames(dfRaw), pattern = "rawData.", replacement = "")
  # colnames(datos) <- gsub(colnames(datos), pattern = "logData.", replacement = "")
  # colnames(dm) <- c("Samples", "Groups")
  
  # Plot graphics
  p0 <- Biomics::plotBarTI(data = dfRaw, interact = F)$grafico
  p1 <- Biomics::plotBoxMulti(base = datos, varResumen = colnames(datos),
                               interact = F, tituloX = "TI distribution")$grafico
  p2 <- Biomics::plotRLE(df = datos, normalizacion = "log", interact = F)$grafico
  p3 <- Biomics::plotMeanSD(df = datos, interact = F)$grafico
  p4 <- Biomics::plotMA(df = datos, dfGrupos = dm, gControl = "G1",
                        gCase = "G2", showR2 = F, interact = F)$grafico
  
  # Get metrics and plot graphics
      # PVC #
  dfPCV <- Biomics::getPCVSimple(dfDatos = datos, grupos = grupos, dfGrupos = dm)
  dfPCV <-  as.data.frame(dfPCV)
  if (length(grupos) < 5) {
    p5 <- Biostatech::plotForest(
      etiquetas = rep(colnames(dfPCV), length(grupos)), 
      estPunt = as.vector(t(as.matrix(dfPCV[seq(1, nrow(dfPCV), 3), ]))), 
      LI = as.vector(t(as.matrix(dfPCV[seq(2, nrow(dfPCV), 3), ]))), 
      LS = as.vector(t(as.matrix(dfPCV[seq(3, nrow(dfPCV), 3), ]))), 
      grupos = rep(grupos, each = ncol(dfPCV)), 
      vertical = T, 
      tituloX = "Mean of the pooled variation coefficient - PVC (%)", 
      referenceLine = F, interact = F)$grafico
  }
  else {
    p5 <- Biomics::plotBoxMulti(base = dfPCV[seq(1, nrow(dfPCV), 3), ], 
                                  varResumen = colnames(dfPCV), 
                                  tituloX = "Normalizations", tituloY = "PVC (%)", 
                                  interact = F)$grafico
  }
      # Correlations #
  correlations <- data.frame(Biomics::getPooledCor(df = datos, dfGrupos = dm, metodo = "spearman"))
  colnames(correlations) <- "cor"
  p6 <- Biostatech::plotBox(base = correlations, tituloX = "Correlation (Spearman)",
                            varResumen = "cor", interact = F)$grafico
  
  # Patch together
  ggpubr::ggarrange(p0, p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 4)
}


getResultsByItem <- function(lista, item){
  # just a common object
  dm <- as.data.frame(lista[[1]][["metadata"]])
  
  # Starting with data extraction and metrics/graphs estimation... 
  if (item == "item0"){ # only item 0
    a <<- 0
    finalPlots <- lapply(lista, function(res){
      a <<- a+1
      nome <- paste0("SimDataset_", a)
      dfRaw <- as.data.frame(res[["rawData"]])
      Biomics::plotBarTI(data = dfRaw, interact = F)$grafico + ggplot2::xlab(nome)
    })
    output <- ggpubr::ggarrange(plotlist = finalPlots, 
                                ncol = 2, nrow = ceiling(length(lista)/2))
  } else { # other items
    # First extrating individual data and generating individual plots
    a <<- 0
    finalData <- lapply(lista, function(res){
      a <<- a+1
      nome <- paste0("SimDataset_", a)
      datos <- as.data.frame(res[["logData"]])
      dm <- as.data.frame(res[["metadata"]])
      if (item == "item6"){
        Biomics::plotBoxMulti(base = datos, varResumen = colnames(datos),
                              interact = F, 
                              tituloX = nome
        )$grafico
      } else if (item == "item5"){
        Biomics::plotRLE(df = datos, normalizacion = nome, interact = F, 
                         tituloX = nome)$grafico
      } else if (item == "item4"){
        Biomics::plotMeanSD(df = datos, interact = F, tituloX = nome)$grafico
      } else if (item == "item3"){
        Biomics::plotMA(df = datos, dfGrupos = dm, gControl = "G1", titulo = nome,
                        gCase = "G2", showR2 = F, interact = F)$grafico
      } else {
        as.data.frame(datos)
      }
    })
    
    names(finalData) <- paste0("SimDataset_",1:length(finalData))
    # Mix plots into a single image or...
    if (item %in% c("item3", "item4", "item5", "item6")){ # plot graphs together: MAplot, RLEplot, meanSDplot, TIboxplot
      output <- ggpubr::ggarrange(plotlist = finalData, 
                                  ncol = 2, 
                                  nrow = ceiling(length(lista)/2))
    } else if (item == "item1"){ # ...or generate other metrics with data for remaining items
      # names(finalData) <- paste0("SimDataset_",1:length(finalData))
      output <- Biomics::getPCV( # PVC
        listData = finalData,
        grupos = unique(dm[,"Groups"]),
        dfGrupos = dm,
        grafico = T,
        interact = F
      )$grafico
    } else if (item =="item2"){ # Correlation
      allVectorsCorr <- lapply(
        finalData, 
        Biomics::getPooledCor, 
        dfGrupos = dm,
        metodo = "spearman")
      dfPlot <- data.frame(sapply(allVectorsCorr, "length<-", max(lengths(allVectorsCorr))))
      output <- Biomics::plotBoxMulti(
        base = dfPlot, 
        varResumen = colnames(dfPlot),
        tituloX = "Normalizations", interact = F,
        tituloY = "Spearman correlation")$grafico
    }
    
  } 
  
  return(output)
}



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
    
    # Correlación (intra > inter) vía factor correlacionado entre muestras
    rho_within  = 0.85,
    rho_between = 0.55,
    loading_sd  = 0.25,
    
    # “Cuña” MA: varianza residual depende de abundancia
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
    
    # Item 0 (shift global por muestra, afecta sumas)
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
  
  # --------- Item0: shift global por muestra ----------
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
                          Groups=factor(groups, levels=c("G1","G2"))),
    item_effects = list(b_shift=b_shift,
                        sigma_i_summary=summary(sigma_i)),
    miss_info = miss_info
  )
}



#.............................................................................
# Executions ####
#.............................................................................

#...............
## Standard ####
results <- simulate_proteomics_clean()
getResults(results)


#......................
## Item 0, 1, 5, 6 ####
# Probas iniciales #
# Opción 1: loading sd
results <- simulate_proteomics_clean(loading_sd = 0.25)
getResults(results)
results <- simulate_proteomics_clean(loading_sd = 1)
getResults(results)
# Opción 2: 
results <- simulate_proteomics_clean(
  sample_shift_sd = 0.1, sample_shift_cap = 0.15)
getResults(results)
results <- simulate_proteomics_clean(
  sample_shift_sd = 0.3, sample_shift_cap = 0.4)
getResults(results)


# Simulación final 
valores <- c(seq(0, 0.49, 0.1), 0.75, 1, 2)
tictoc::tic()
resByItem <- lapply(valores, function(x) simulate_proteomics_clean(
  semilla = 1000, 
  n_proteins = 1000,
  sample_shift_sd = x,
  sample_shift_cap = x+0.05))
tictoc::toc()
getResultsByItem(resByItem, item = "item0")
getResultsByItem(resByItem, item = "item6")
getResultsByItem(resByItem, item = "item5")
getResultsByItem(resByItem, item = "item1") # influyeeee
# getResultsByItem(resByItem, item = "item3") # algo influye si , pero creo que non vai ser a mellor forma de medilo
# getResultsByItem(resByItem, item = "item2") # Aumenta lixeiramente a correlación si... moi lixeiramente


#......................
## Item 2 ####
# Probas iniciales #
results <- simulate_proteomics_clean(
  rho_between = 0.1, 
  rho_within = 1) # Valores baixos de between e altos de within dan boa correlacion
getResults(results)
results <- simulate_proteomics_clean(
  rho_between = 1, 
  rho_within = 0.1)
getResults(results)

# Xa sabemos como funciona, ahora a usar varios valores #
values_rho_between <- c(seq(0, 1, 0.25), 1.5, 2, 3)
values_rho_within <- rev(c(0, 0.05, 0.2, seq(0.5, 1.5, 0.25)))
tictoc::tic()
resByItem <- lapply(1:length(values_rho_between), function(x) simulate_proteomics_clean(
  semilla = 10000, 
  n_proteins = 1000,
  rho_between = values_rho_between[x],
  rho_within = values_rho_within[x]))
tictoc::toc()
getResultsByItem(resByItem, item = "item2")
getResultsByItem(resByItem, item = "item1") # afecta para betw<0.5+within>0.5 respecto a betw>0.5+within<0.5
getResultsByItem(resByItem, item = "item0")
# getResultsByItem(resByItem, item = "item3")

# Alternativa #
sample_sd_strength <- c(seq(0, 1, 0.25), 1.5, 2, 3)
tictoc::tic()
resByItem <- lapply(1:length(values_rho_between), function(x) simulate_proteomics_clean(
  semilla = 10000, 
  n_proteins = 1000,
  sample_sd_strength = sample_sd_strength[x]))
tictoc::toc()
getResultsByItem(resByItem, item = "item2")


#......................
## Item 4 ####
# Probas iniciales #
results <- simulate_proteomics_clean(
  sample_shift_sd = 0.1,
  sample_sd_strength = 1,
  sample_sd_rho = 1.5, 
  sample_sd_cap = 0) # Aumentando este valor aumenta a pendiente das rectas
getResults(results)


# Xa sabemos como funciona, ahora a usar varios valores #
sample_sd_strength <- c(seq(0, 1, 0.25), 1.5, 2, 3)
sample_sd_cap <- c(seq(0, 0.5, 0.05), 1.5, 2, 3)
tictoc::tic()
resByItem <- lapply(1:length(sample_sd_strength), function(x) simulate_proteomics_clean(
  semilla = 10000, 
  n_proteins = 1000,
  sample_shift_sd = 0.5,  #sample_shift_cap = 0.15,
  sample_sd_strength = 1,
  sample_sd_rho = 1.5, 
  sample_sd_cap = sample_sd_cap[x]))
tictoc::toc()
getResultsByItem(resByItem, item = "item4")
getResultsByItem(resByItem, item = "item3")

# opcion 2 así flipas como cambia e mais aleatorio 
sample_sd_cap <- seq(1, 3, 1)
tictoc::tic()
resByItem <- lapply(1:length(sample_sd_cap), function(x) simulate_proteomics_clean(
  sample_shift_sd = 0.5,
  sample_sd_strength = 2,
  sample_sd_rho = 0, 
  sample_sd_cap = sample_sd_cap[x]))
tictoc::toc()
getResultsByItem(resByItem, item = "item4")
getResultsByItem(resByItem, item = "item3")
getResultsByItem(resByItem, item = "item0")
getResultsByItem(resByItem, item = "item1")





#......................
## Item 3 ####
# Probas iniciales #
results <- simulate_proteomics_clean(
  prop_de = 0.1,
  sigma_lo = 0.6, # poñendo valores crecientes inversos entre este argumento e o seguinte invírtese a forma de cuña
  sigma_hi = 2,
  gamma_sigma = 3) 
getResults(results)
results <- simulate_proteomics_clean(
  sample_shift_sd = 0.5,
  sample_sd_strength = 2,
  sample_sd_rho = 0, 
  sample_sd_cap = 1) # Aumentando este valor aumenta a pendiente das rectas (co resto de parámetros tal cual)
getResults(results)


# Xa sabemos como funciona, ahora a usar varios valores #
# Opcion 1 - cambio forma #
sigma_lo <- c(seq(0, 1, 0.1), 1.5, 2, 3)
sigma_hi <- rev(c(seq(0, 1, 0.1), 1.5, 2, 3))
tictoc::tic()
resByItem <- lapply(1:length(sigma_lo), function(x) simulate_proteomics_clean(
  semilla = 10000, 
  n_proteins = 1000,
  prop_de = 0.1,
  sigma_lo = sigma_lo[x],
  sigma_hi = sigma_hi[x],
  gamma_sigma = 3))
tictoc::toc()
getResultsByItem(resByItem, item = "item3")

# Opcion 2 - cambio pendiente #
sample_sd_cap <- -1*c(seq(0, 1, 0.25), 1.5, 2, 3) # Cambia o ancho dos puntos (canto más grande mais estreito)
sample_sd_cap <- c(seq(0, 1, 0.25), 1.5, 2, 3)
tictoc::tic()
resByItem <- lapply(1:length(sample_sd_cap), function(x) simulate_proteomics_clean(
  sample_shift_sd = 0.5,
  sample_sd_strength = 3,
  sample_sd_rho = 0, 
  sample_sd_cap = sample_sd_cap[x]))
tictoc::toc()
getResultsByItem(resByItem, item = "item3")






##
#......................
# other trials with different arguments
results <- simulate_proteomics_clean(
  n_proteins = 5000,
  enforce_rle = FALSE,
  add_missing = TRUE,
  
  gamma_sigma = 3.5,
  sigma_hi = 0.04,
  sigma_lo = 0.70,
  
  prop_de = 0.08,
  logFC_sd = 0.75,
  hetero_logFC = TRUE,
  fc_hi = 0.15,
  fc_lo = 2.5,
  gamma_fc = 7,
  
  loading_sd = 0.25
)
a <- 1.5
results <- simulate_proteomics_clean(
  n_proteins = 5000,
  enforce_rle = FALSE,
  add_missing = TRUE,
  
  gamma_sigma = 4.0*a,
  sigma_hi = 0.03*a,
  sigma_lo = 0.85*a,
  
  prop_de = 0.05*a,
  logFC_sd = 1.0*a,
  hetero_logFC = TRUE,
  fc_hi = 0.10*a,
  fc_lo = 2.5*a,
  gamma_fc = 7*a,
  
  loading_sd = 0.25*a
)








