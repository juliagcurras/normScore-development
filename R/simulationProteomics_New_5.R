
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
  dfRaw <- as.data.frame(results["rawData"])
  colnames(dfRaw) <- gsub(colnames(dfRaw), pattern = "rawData.", replacement = "")
  datos <- as.data.frame(results["logData"])
  colnames(datos) <- gsub(colnames(datos), pattern = "logData.", replacement = "")
  dm <- as.data.frame(results["metadata"])
  colnames(dm) <- c("Samples", "Groups")
  
  p1 <- Biomics::plotBarTI(data = dfRaw, interact = F)$grafico
  p11 <- Biomics::plotBoxMulti(base = datos, varResumen = colnames(datos),
                               interact = F, tituloX = "TI distribution")$grafico
  p2 <- Biomics::plotRLE(df = datos, normalizacion = "log", interact = F)$grafico
  p3 <- Biomics::plotMeanSD(df = datos, interact = F)$grafico
  p4 <- Biomics::plotMA(df = datos, dfGrupos = dm, gControl = "G1",
                        gCase = "G2", showR2 = F, interact = F
                        # limY = c(-2,2)
  )$grafico
  correlations <- data.frame(Biomics::getPooledCor(df = datos, dfGrupos = dm, metodo = "spearman"))
  colnames(correlations) <- "cor"
  p5 <- Biostatech::plotBox(base = correlations, tituloX = "Correlation (Spearman)",
                            varResumen = "cor", interact = F)$grafico
  
  ggpubr::ggarrange(p1, p11, p2, p3, p4, p5, ncol = 2, nrow = 3)
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
    n_proteins = 10000, # 1000,
    n_per_group = 20,
    
    # Medias por proteína (log2)
    mu_mean = 18,
    mu_sd   = 1.2,
    mu_clip = c(15, 25),
    
    # Correlación (intra > inter) vía factor correlacionado entre muestras
    rho_within  = 0.85,
    rho_between = 0.55,
    loading_sd  = 0.25,   # cuánto pesa el factor sobre proteínas (0.4)
    
    # “Cuña” MA: varianza residual depende de abundancia
    sigma_hi = 0.05,      # ruido en alta expresión (0.18)
    sigma_lo = 0.4,      # ruido en baja expresión (0.75)
    gamma_sigma = 2.5,    # mayor => cuña más marcada en MAplot (2.6)
    
    # DE (simétrica entre grupos)
    prop_de = 0.5, # 0.05
    logFC_sd = 1, # 2.0
    logFC_mean = 0,
    hetero_logFC = TRUE,  # DE más grandes en baja expresión
    fc_hi = 0.55, # 0.55
    fc_lo = 2.5, # 1.55
    gamma_fc = 7, # 2.0
    
    # Item 0
    sample_shift_sd = 0, 
    sample_shift_cap = 0.2,
    
    # Missing (por defecto OFF para “sin error”)
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
  
  # peso "baja expresión": 1 en baja (izqda), 0 en alta (dcha)
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
  
  # vector de aplicación simétrica
  gvec <- ifelse(groups == "G1", +0.5, -0.5)  # +logFC/2 y -logFC/2
  DE_mat <- outer(logFC, gvec)
  
  
  # --------- factor correlacionado para correlación entre muestras ----------
  # Construimos una correlación por bloques y sampleamos un score s_j
  R <- matrix(rho_between, m, m); diag(R) <- 1
  idx1 <- which(groups=="G1"); idx2 <- which(groups=="G2")
  R[idx1, idx1] <- rho_within; diag(R[idx1, idx1]) <- 1
  R[idx2, idx2] <- rho_within; diag(R[idx2, idx2]) <- 1
  
  ev <- eigen(R, symmetric=TRUE, only.values=TRUE)$values
  if (min(ev) <= 1e-8) R <- R + diag(abs(min(ev)) + 1e-6, m)
  
  s <- as.numeric(MASS::mvrnorm(1, mu = rep(0, m), Sigma = R))
  
  # CLAVE: centramos s dentro de cada grupo para que NO meta diferencia en logFC por grupos
  s[idx1] <- s[idx1] - mean(s[idx1])
  s[idx2] <- s[idx2] - mean(s[idx2])
  
  load <- rnorm(n_proteins, 0, loading_sd)
  FACT_mat <- outer(load, s)
  
  
  # --------- ruido residual heterocedástico (cuña MA) ----------
  sigma_i <- sigma_hi + (w_low^gamma_sigma) * (sigma_lo - sigma_hi)
  EPS <- matrix(rnorm(n_proteins*m), nrow=n_proteins, ncol=m) * sigma_i
  
  
  # --------- matriz final (log2) ----------
  X <- matrix(mu, nrow=n_proteins, ncol=m) + FACT_mat + DE_mat + EPS
  
  # ---- Item0: desajuste en suma por muestra (offset en log2) ----
  # if (sample_shift_sd > 0 || sample_shift_trend != 0) {
  #   j <- seq_len(m)
  #   b <- rnorm(m, 0, sample_shift_sd) + sample_shift_trend * scale(j)[,1]
  #   b <- as.numeric(b - mean(b))     # que el dataset global no se desplace
  #   X <- sweep(X, 2, b, "+")         # offset log2 por muestra
  # }
  
  # ---- Item0: desajuste en suma por muestra (offset log2 aleatorio) ----
  if (sample_shift_sd > 0) {
    b <- rnorm(m, mean = 0, sd = sample_shift_sd)
    
    # quitar media global para no desplazar el dataset entero
    b <- b - mean(b)
    
    # cap suave para evitar 1-2 muestras dominantes (winsorize)
    if (!is.null(sample_shift_cap) && is.finite(sample_shift_cap)) {
      b <- pmin(pmax(b, -sample_shift_cap), sample_shift_cap)
      b <- b - mean(b)  # re-centrar tras cap
    }
    
    X <- sweep(X, 2, b, "+")
  }
  
  # --------- Aesthetics ----------
  rownames(X) <- paste0("P", sprintf("%05d", 1:n_proteins))
  colnames(X) <- paste0(groups, "_", ave(seq_along(groups), groups, FUN = seq_along))
  
  
  # --------- missing MNAR opcional (self-contained) ----------
  miss_info <- NULL
  if (add_missing) {
    b <- rnorm(m, 0, missing_by_sample_sd)  # efecto por muestra (no por grupo)
    b <- b - mean(b)
    
    # calibrar intercepto a para cumplir target_missing
    # p_ij = sigmoid(a - k*X_ij + b_j)
    f <- function(a){
      p <- plogis(a - k_mnar * X + matrix(b, nrow=n_proteins, ncol=m, byrow=TRUE))
      mean(p, na.rm=TRUE) - target_missing
    }
    a_hat <- uniroot(f, interval=c(-50, 50))$root
    
    P <- plogis(a_hat - k_mnar * X + matrix(b, nrow=n_proteins, ncol=m, byrow=TRUE))
    M <- matrix(runif(n_proteins*m), nrow=n_proteins, ncol=m) < P
    X[M] <- NA_real_
    
    miss_info <- list(P = P, mask = M)
  }
  
  list(
    logData  = X,
    rawData  = 2^X,
    metadata = data.frame(Samples=colnames(X), Groups=factor(groups, levels=c("G1","G2")))
    # de_info  = data.frame(ProteinID=rownames(X)[de_idx], logFC_expected=logFC[de_idx]),
    # sim_info = list(rho_within=rho_within, rho_between=rho_between, sigma_summary=summary(sigma_i)),
    # miss_info = miss_info
  )
}



#.............................................................................
# Executions ####
#.............................................................................

#...............
## Standard ####
results <- simulate_proteomics_clean()
getResults(results)
dfRaw <- as.data.frame(results[["rawData"]])
datos <- as.data.frame(results[["logData"]])
dm <- as.data.frame(results[["metadata"]])


#......................
## Item 0, 1, 5, 6 ####
valores <- c(seq(0, 0.49, 0.1), 0.75, 1, 2)
tictoc::tic()
resByItem <- lapply(valores, function(x) simulate_proteomics_clean(
  semilla = 1000, 
  n_proteins = 100,
  sample_shift_sd = x,
  sample_shift_cap = x+0.05))
  # loading_sd = x, semilla = 1000))
tictoc::toc()
getResultsByItem(resByItem, item = "item0")
getResultsByItem(resByItem, item = "item6")
getResultsByItem(resByItem, item = "item5")
getResultsByItem(resByItem, item = "item3") # algo influye si , pero croe que non vai ser a mellor forma de medilo
getResultsByItem(resByItem, item = "item1") # influyeeee
getResultsByItem(resByItem, item = "item2") # Aumenta lixeiramente a correlación si... moi lixeiramente



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








