
# Julia G Curras - 07/10/2025
# Simulación de datos de proteómica (log2-scal e)
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/R")

# if(!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)
library(future.apply)
library(progressr)
library(dplyr)
source(file = "supportFunctions.R", encoding = "UTF-8")
source(file = "scoreFunction.R", encoding = "UTF-8")

pathToData <- "C:/Users/julia/Documents/GitHub/normScore/Simulations/"
norm <- "mean"
set.seed(9396)


#.........................................................................####
# Function ####
#.........................................................................####
simulate_mean <- 
  function(
    n_proteins = 1000, # total proteínas
    n_per_group = 10, # muestras por grupo
    k_groups = 2, # grupos
    p_high = 0.02, # proporción de proteínas "muy abundantes" (colas)
    mu_low_mean = 18, # centro de la masa principal
    mu_low_sd   = 1.2,
    mu_high_mean = 23, # centro de las pocas proteínas muy altas
    mu_high_sd   = 0.6,
    rho_samples = 0.6, # efectos por muestra para correlación
    sigma_sample = 1.0, # escala del efecto por muestra
    loading_sd = 0.6, # "carga" de cada proteína frente a los efectos de muestra (heterogeneidad entre proteínas)
    sigma_resid = 0.6, # añadir variación residual por proteína-muestra
    prop_de = 0.05, # proporción de proteínas DE
    logFC_mean = 2.0, # distribución de valores de logFC de proteínas DE
    logFC_sd   = 0.6, 
    sd_tech_effect = 1.5, # variabilidad del efecto técnico
    semilla = NULL
  ){
    if (!is.null(semilla)) set.seed(semilla)
    # ----- Novos parámetros -----
    # groups <- rep(c("G1","G2"), each = n_per_group)
    groups <- paste0("G", rep(1:k_groups, length.out = n_per_group * k_groups))
    m <- length(groups)       # nº total de muestras = 20

    # mezcla para medias proteicas (log2): muchas medias bajas/intermedias, pocas muy altas
    mu <- rnorm(n_proteins, mean = mu_low_mean, sd = mu_low_sd)
    high_idx <- sample(1:n_proteins, size = round(n_proteins * p_high))
    mu[high_idx] <- rnorm(length(high_idx), mean = mu_high_mean, sd = mu_high_sd)

    # truncar a rango razonable (15-25) para evitar extremos raros
    mu <- pmin(pmax(mu, 15), 25)

    # ----- estructura de correlación entre muestras -----
    # creamos un vector de efectos por muestra (común a todas las proteínas) para generar correlación
    Sigma <- matrix(rho_samples, nrow = m, ncol = m)
    diag(Sigma) <- 1
    sample_effects <- as.numeric(mvrnorm(n = 1, mu = rep(0, m), Sigma = Sigma)) * sigma_sample
    # sample_effects es un vector (length m) que da la estructura correlada entre muestras

    # Cada proteína tiene una "carga" frente a esos efectos de muestra (heterogeneidad entre proteínas)
    # generamos cargas y aseguramos que sean ortogonales a 'mu'
    loadings_raw <- rnorm(n_proteins, mean = 0, sd = loading_sd)
    # eliminar cualquier correlación lineal con 'mu' (conservar variabilidad)
    loadings <- resid(lm(loadings_raw ~ mu))

    # ----- definir proteínas diferencialmente abundantes (DE) -----
    # Seleccionamos preferentemente DE entre proteínas con medias NO MUY ALTAS
    candidate_idx <- setdiff(1:n_proteins, high_idx) # evitamos las muy altas
    n_de <- round(n_proteins * prop_de)
    de_idx <- sample(candidate_idx, size = n_de)

    # generar log2 fold-changes para DE (valor absoluto grande)
    # los DE tendrán logFC positivos o negativos con signo aleatorio
    de_logFC <- rnorm(n_de, mean = logFC_mean, sd = logFC_sd) * sample(c(-1,1), n_de, replace = TRUE)
    
    group_effect <- rep(0, n_proteins)
    group_effect[de_idx] <- de_logFC
    # group_effect es el efecto que se suma a G1 (si quieres G2 vs G1 invierte la asignación)
    
    # ----- construir la matriz (proteínas x muestras) -----
    mat <- matrix(NA, nrow = n_proteins, ncol = m)
    for(i in 1:n_proteins){
      # efecto base: la media de la proteína
      base <- mu[i]
      # efecto de muestra para esta proteína (loadings[i] * sample_effects)
      samp_eff_i <- loadings[i] * sample_effects
      # efecto de grupo: añadir group_effect[i] sólo para las muestras de G1
      grp_eff <- ifelse(groups == "G1", group_effect[i], 0)
      # ruido residual
      resid <- rnorm(m, mean = 0, sd = sigma_resid)
      mat[i, ] <- base + samp_eff_i + grp_eff + resid
    }
    rownames(mat) <- paste0("P", sprintf("%05d", 1:n_proteins))
    colnames(mat) <- paste0(groups, "_", rep(1:n_per_group, k_groups))
    
    
    # ------ Añadir efecto técnico por muestra -----
    # simulamos un sesgo técnico que afecte a cada muestra, (por ejemplo, distinta eficiencia de cuantificación)
    tech_effect <- rnorm(m, mean = 0, sd = sd_tech_effect)  # cambia sd para ajustar magnitud del sesgo
    # protein_sensitivity <- rnorm(n_proteins, 1, 0.1) # variación entre proteínas (más real, peor para mean)
    names(tech_effect) <- colnames(mat)
    
    # aplicamos el efecto técnico sumando un desplazamiento por muestra
    # mat <- mat + protein_sensitivity %*% t(tech_effect) # variación entre proteínas
    mat <- sweep(mat, 2, tech_effect, FUN = "+")
    
    # ----- Matrices finales -----
    meta <- data.frame(
      Samples = colnames(mat),
      Groups = factor(groups)
    )
    matRaw <- 2^mat
    
    # ----- Guardar -----
    # write.csv(matRaw, file = paste0(pathToData, norm, "/simulated_matrix.csv"), row.names = TRUE)
    # write.csv(meta, file = paste0(pathToData, norm, "/simulated_design_matrix.csv"), row.names = FALSE)
    
    return(list(
      rawData = matRaw,
      logData = as.data.frame(mat),
      metadata = meta
    ))
  }




#.........................................................................####
# Multiple simulations ####
#.........................................................................####

## Settings ####
n_proteins <- c(1000, 5000, 10000)
n_per_group <- c(5, 15, 20)
k_groups <- c(2, 3, 4)
sigma_resid <- c(0.2, 0.6, 1) # variación residual
sd_tech_effect <- c(0.75, 1.5, 2.5) # variación por efecto técnico a corregir con media

## All combinations of settings ####
grid <- expand.grid(
  n_proteins = n_proteins,
  n_per_group = n_per_group,
  k_groups = k_groups,
  sigma_resid = sigma_resid,
  sd_tech_effect = sd_tech_effect
)
nrow(grid)  # 243 combinacions
saveRDS(grid, file = paste0(pathToData, norm, "/optionsSim_243.rds"))

## Execution of simulations ####
set.seed(9396)
handlers(global = TRUE)
handlers("txtprogressbar")  
plan(multisession, workers = parallel::detectCores() - 16)
ordenNorm <- c("CyclicLoess", "GI", "Log","MAD", "Mean", "Median", "Quantile", "RLR", "VSN")
# datasets <- list()

with_progress({
  cat("\n\tIniciando simulacións e comparacións de normalizacións... \n\n")
  p <- progressor(steps = nrow(grid))
  results <- future_lapply(seq_len(nrow(grid)), function(i) {
    # results <- future_lapply(1:50, function(i) {
    # resultsOri <- sapply(c(1, 4:5), function(i) {
    # Simulate data
    params <- grid[i, ]
    sim <- simulate_mean(
      n_proteins = params$n_proteins,
      n_per_group =  params$n_per_group,
      k_groups = params$k_groups,
      sigma_resid = params$sigma_resid,
      sd_tech_effect = params$sd_tech_effect,
      semilla = 1000 + i
    )
    # datasets[[i]] <- sim
    
    # Retrieving data
    dfGrupos <- sim$metadata
    logIntensityMatrix <- sim$logData
    intensityMatrix <- sim$rawData
    grupos <- levels(as.factor(dfGrupos$Groups))
    
    # Normalization...
    mydata <- list(Log = logIntensityMatrix,
                   Mean = meanNorm(rawMatrix = intensityMatrix), 
                   Median = medianNorm(rawMatrix = intensityMatrix), 
                   GI = GINorm(rawMatrix = intensityMatrix),
                   Quantile = quantileNorm(log2Matrix = logIntensityMatrix),
                   VSN = VSNNorm(rawMatrix = intensityMatrix),
                   CyclicLoess = cyclicLoessNorm(log2Matrix = logIntensityMatrix),
                   RLR = RLRNorm(log2Matrix = logIntensityMatrix),
                   MAD = MADNormalization(log2Matrix = logIntensityMatrix))
    
    # Assessing normalization
    finalRank <- normScore( 
      designMatrix = dfGrupos, 
      normMatrixList = mydata, 
      dfRaw = intensityMatrix, 
      refGroup = "G2", altGroup = "G1")$finalRanking
    
    # finalRank <- finalRank[ordenNorm]
    p()  # avanza la barra
    return(finalRank)
  })
  # }, simplify = F)
  cat("\n\tSimulacións rematadas. \n\n")
})

# saveRDS(results, file = paste0(pathToData, norm, "/results_sim243_2.RDS"))
saveRDS(results, file = paste0(pathToData, norm, "/results_sim243_3.RDS"))
# results <- readRDS(file = paste0(pathToData, norm, "/results_sim243.RDS"))
# saveRDS(datasets, file = paste0(pathToData, norm, "/datasets_sim243.RDS"))

## Results for normalizations ####
dfRes <- sapply(results, function(i) i[ordenNorm])
# dfRes <- as.data.frame(dplyr::bind_cols(datasets))
colnames(dfRes) <- paste0("Sim ", 1:ncol(dfRes))
rownames(dfRes) <- ordenNorm

# Top 1 norm: mean is the second one
topNorm <- sapply(results, function(i) names(sort(i))[1])
table(topNorm)
# Top 2 norm: mean is the second one
top2Norm <- sapply(results, function(i) "Mean" %in% names(sort(i))[1:2])
table(top2Norm)
# Top 3 norm: mean is the second one
top3Norm <- sapply(results, function(i) "Mean" %in% names(sort(i))[1:3])
table(top3Norm)
# Top 5 norm: mean is the second one
top5Norm <- sapply(results, function(i) "Mean" %in% names(sort(i))[1:5])
table(top5Norm)


#
















# _ ####
# _ ####
# _ ####
#.........................................................................####
# COMPLETE AND ORIGINAL CODE ####
# ----- parámetros -----
n_proteins <- 1000        # número de proteínas (usé 5444 del ejemplo)
n_per_group <- 10         # 10 muestras por grupo
groups <- rep(c("G1","G2"), each = n_per_group)
m <- length(groups)       # nº total de muestras = 20

# mezcla para medias proteicas (log2): muchas medias bajas/intermedias, pocas muy altas
p_high <- 0.02            # proporción de proteínas "muy abundantes" (colas)
mu_low_mean <- 18         # centro de la masa principal
mu_low_sd   <- 1.2
mu_high_mean <- 23        # centro de las pocas muy altas
mu_high_sd   <- 0.6

mu <- rnorm(n_proteins, mean = mu_low_mean, sd = mu_low_sd)
high_idx <- sample(1:n_proteins, size = round(n_proteins * p_high))
mu[high_idx] <- rnorm(length(high_idx), mean = mu_high_mean, sd = mu_high_sd)

# truncar a rango razonable (15-25) para evitar extremos raros
mu <- pmin(pmax(mu, 15), 25)

# ----- estructura de correlación entre muestras -----
# creamos un vector de efectos por muestra (común a todas las proteínas) para generar correlación
rho_samples <- 0.6   # correlación entre muestras inducida por el factor común
sigma_sample <- 1.0  # escala del efecto por muestra

Sigma <- matrix(rho_samples, nrow = m, ncol = m)
diag(Sigma) <- 1
sample_effects <- as.numeric(mvrnorm(n = 1, mu = rep(0, m), Sigma = Sigma)) * sigma_sample
# sample_effects es un vector (length m) que da la estructura correlada entre muestras

# Cada proteína tiene una "carga" frente a esos efectos de muestra (heterogeneidad entre proteínas)
loading_sd <- 0.6
# generamos cargas y aseguramos que sean ortogonales a 'mu'
loadings_raw <- rnorm(n_proteins, mean = 0, sd = loading_sd)
# eliminar cualquier correlación lineal con 'mu' (conservar variabilidad)
loadings <- resid(lm(loadings_raw ~ mu))

# ----- añadir variación residual por proteína-muestra -----
sigma_resid <- 0.6  # ruido independiente

# ----- definir proteínas diferencialmente abundantes (DE) -----
prop_de <- 0.05  # 5% DE (ajusta si quieres más/menos)
# Seleccionamos preferentemente DE entre proteínas con medias NO MUY ALTAS
candidate_idx <- setdiff(1:n_proteins, high_idx) # evitamos las muy altas
n_de <- round(n_proteins * prop_de)
de_idx <- sample(candidate_idx, size = n_de)

# generar log2 fold-changes para DE (valor absoluto grande)
# los DE tendrán logFC positivos o negativos con signo aleatorio
logFC_mean <- 2.0
logFC_sd   <- 0.6
de_logFC <- rnorm(n_de, mean = logFC_mean, sd = logFC_sd) * sample(c(-1,1), n_de, replace = TRUE)

group_effect <- rep(0, n_proteins)
group_effect[de_idx] <- de_logFC
# group_effect es el efecto que se suma a G1 (si quieres G2 vs G1 invierte la asignación)

# ----- construir la matriz (proteínas x muestras) -----
mat <- matrix(NA, nrow = n_proteins, ncol = m)
for(i in 1:n_proteins){
  # efecto base: la media de la proteína
  base <- mu[i]
  # efecto de muestra para esta proteína (loadings[i] * sample_effects)
  samp_eff_i <- loadings[i] * sample_effects
  # efecto de grupo: añadir group_effect[i] sólo para las muestras de G1
  grp_eff <- ifelse(groups == "G1", group_effect[i], 0)
  # ruido residual
  resid <- rnorm(m, mean = 0, sd = sigma_resid)
  mat[i, ] <- base + samp_eff_i + grp_eff + resid
}

rownames(mat) <- paste0("P", sprintf("%05d", 1:n_proteins))
colnames(mat) <- paste0(groups, "_", rep(1:n_per_group, 2))


# ===== Añadir efecto técnico por muestra =====
# simulamos un sesgo técnico que afecte a cada muestra
# (por ejemplo, distinta eficiencia de cuantificación)
tech_effect <- rnorm(m, mean = 0, sd = 1.5)  # cambia sd para ajustar magnitud del sesgo
# protein_sensitivity <- rnorm(n_proteins, 1, 0.1)
names(tech_effect) <- colnames(mat)

# aplicamos el efecto técnico sumando un desplazamiento por muestra
# mat <- mat + protein_sensitivity %*% t(tech_effect)
mat <- sweep(mat, 2, tech_effect, FUN = "+")


# ===== Matrices finales =====
meta <- data.frame(
  Samples = colnames(mat),
  Groups = factor(groups)
)
matRaw <- 2^mat

# Opcional: guardar
write.csv(matRaw, file = paste0(pathToData, norm, "/simulated_matrix.csv"), row.names = TRUE)
write.csv(meta, file = paste0(pathToData, norm, "/simulated_design_matrix.csv"), row.names = FALSE)
# write.csv(data.frame(Protein=rownames(mat), mu=mu, isHigh = (1:n_proteins %in% high_idx),
#                      isDE = (1:n_proteins %in% de_idx), true_logFC = group_effect),
#           file = "simulated_proteins_annotation.csv", row.names = FALSE)

# Fin del script: devuelve objetos 'mat', 'mu', 'de_idx', 'high_idx', 'est_logFC'
