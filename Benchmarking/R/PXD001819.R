###############################################################################/
##################             UPS BENCHMARKING             ###################/
###############################################################################/

# Set directory
setwd("C:/Users/julia/Documents/GitHub/normScore-development/Benchmarking/R")
rm(list=ls())
graphics.off()

# Libraries ####
library(dplyr)
library(ggplot2)


# Functions ####

mae <- function(actual, predicted){
  mean(abs(actual - predicted))
}


getROC <- function(df, thresholds = seq(0, 1, by = 0.001)) {
  
  rocData <- do.call(rbind, lapply(thresholds, function(threshold) {
    
    predicted <- df$adj.P.Val <= threshold
    
    TP <- sum(predicted & df$truth == 1)
    FP <- sum(predicted & df$truth == 0)
    TN <- sum(!predicted & df$truth == 0)
    FN <- sum(!predicted & df$truth == 1)
    
    data.frame(
      threshold = threshold,
      TPR = TP / (TP + FN),
      FPR = FP / (FP + TN)
    )
  }))
  
  rocData
}



evalAUC <- function(dfDE, upsFC, yeastFC = 1){
  
  # Get initial dataset #
  dfDE$Protein.Ids <- rownames(dfDE)
  dfDE$Organism <-   grepl(pattern = "ups", x = rownames(dfDE))
  dfDE$Organism <- factor(dfDE$Organism, 
                          levels = c(TRUE, FALSE), 
                          labels = c("UPS", "YEAST"))
  df <- dfDE
  
  ## A) AUC and ROC curves 
  
  ## AUC sobre todos los puntos de corte relevantes
  # Estimate metrics #
  df$truth <- ifelse(is.na(df$Organism), NA, 
                     ifelse(df$Organism == "YEAST", 0, 1)) # Proteína de levadura non DE (0)
  
  rocObj <- pROC::roc(
    response = df$truth,
    predictor = df$adj.P.Val,
    levels = c(0, 1),
    direction = ">"
  )
  
  aucValue <- as.numeric(pROC::auc(rocObj))
  
  
  ## ROC data sobre una serie de puntos concretos para que sean comunes a todas 
  # las normalizaciones. 
  rocData <- getROC(df = df, thresholds = thresholdsROC)
  
  ## B) Regression 
  ### MAE 
  #### Expected values
  df$logFC_Real <- ifelse(is.na(df$Organism), NA,  
                          ifelse(df$Organism == "UPS", log(upsFC, 2), 
                                 ifelse(df$Organism == "YEAST", log(yeastFC, 2), "Error"
                                 )
                          )
  )
  df$logFC_Real <- as.numeric(df$logFC_Real)
  
  
  #### MAE global 
  maeG <- mean(abs(df$logFC_Real - df$logFC))
  
  #### MAE by class
  # df$Organism <- factor(df$Organism, levels = c("HUMAN", "YEAST", "ECOLI"))
  tablaC <- df %>% 
    group_by(Organism) %>% 
    summarise(MAE = mae(actual = logFC_Real, predicted = logFC)) %>%
    as.data.frame()
  maeGlobal <- maeG

  ## Returning 
  finalResult <- list(
    AUC = aucValue,
    dfROC = rocData, 
    maeGlobal = maeGlobal,
    tablaMaeOrg = tablaC
  )
  
  return(finalResult)
}




#............................................................................----
# Loading data ####

data <- readxl::read_excel(
  path = "../Data/PXD001819/MaxQuant.xlsx", 
  sheet = 1, 
  col_names = TRUE
)

# Comparaciones
table(data$`Original comparison`)
table(
  data$`Original comparison`,
  data$Species
)

dm <- data.frame(
  Samples = paste0(rep(LETTERS[1:2], each = 3), 1:3), 
  Groups = rep(LETTERS[1:2], each = 3)
)
dm

#................................................................................
# Exploring 50vs0.5 => FC = 100 ####
upsFC <- 100

df100 <- data |> 
  dplyr::filter(`Original comparison` == "50Vs0.5") |> 
  dplyr::select(Species, `Majority protein IDs`, `A1 (raw)`:`B3 (raw)`)
colnames(df100) <- gsub(colnames(df100), pattern = " (raw)", replacement = "", fixed = TRUE)
colnames(df100)
table(df100$Species)

df100[df100 == 0] <- NA

## Assesing missingness ####
data100 <- df100 |> dplyr::filter(Species != "UPS") |> dplyr::select(A1:B3)
data100 <- df100 |> dplyr::filter(Species == "UPS") |> dplyr::select(A1:B3)
data100 <- as.data.frame(data100)


Biomics::plotBarTI(data = data100, interact = FALSE, color = "red")
pNA <- Biomics::plotBarNA(data = data100, interact = FALSE)
pNA$dfMuestra
pNA$graficoMuestra
pNA$graficoProteina

pheatmap::pheatmap(mat = data100, show_rownames = FALSE, 
                   cluster_rows = FALSE, cluster_cols = FALSE)

# B1 mal, B2 y B3 con mucho missing para UPS como cabría de esperar. 


#................................................................................
# Exploring 25vs12.5 => FC = 2 ####
# upsFC <- 2
# 
# df2 <- data |> 
#   dplyr::filter(`Original comparison` == "25Vs12.5") |> 
#   dplyr::select(Species, `Majority protein IDs`, `A1 (raw)`:`B3 (raw)`)
# colnames(df2) <- gsub(colnames(df2), pattern = " (raw)", replacement = "", fixed = TRUE)
# colnames(df2)
# 
# df2[df2 == 0] <- NA
# 
# ## Assesing missingness ####
# # data2 <- df2 |> dplyr::filter(Species != "UPS") |> dplyr::select(A1:B3)
# data2 <- df2 |> dplyr::filter(Species == "UPS") |> dplyr::select(A1:B3)
# data2 <- as.data.frame(data2)
# 
# 
# Biomics::plotBarTI(data = data2, interact = FALSE)
# pNA <- Biomics::plotBarNA(data = data2, interact = FALSE)
# pNA$dfMuestra
# pNA$graficoMuestra
# pNA$graficoProteina
# 
# pheatmap::pheatmap(mat = data2, show_rownames = FALSE, 
#                    cluster_rows = FALSE, cluster_cols = FALSE)
# 
# ### Filtering ####
# 
# # Common processing
# data2 <- df2 |> dplyr::select(A1:B3)
# data2 <- as.data.frame(data2)
# 
# rownames(data2) <- gsub(
#   x = sapply(
#     strsplit(
#       x = df2$`Majority protein IDs`, split = "|", fixed = TRUE),
#     "[[", 1),
#   pattern = ">",
#   replacement = "",
#   fixed = TRUE)
# 
# Biomics::filterMissing(
#   df = data2, 
#   dfGrupos = dm, 
#   threshold = c(0, 1, 2, 3), 
#   mode = "minDataByGroup"
# )$tablaFormato
# 
# dataFinal <- Biomics::filterMissing(
#   df = data2, 
#   dfGrupos = dm, 
#   threshold = 2, 
#   mode = "minDataByGroup"
# )$tabla
# 
# rownames(dataFinal)[1:2]
# 
# dataLogFinal <- log(dataFinal, base = 2)
# 
# ## Log and visualization  ####
# pheatmap::pheatmap(dataFinal, show_rownames = F)
# apply(dataLogFinal, 2, min, na.rm = TRUE)
# pheatmap::pheatmap(dataLogFinal, show_rownames = F, scale = "row")
# pheatmap::pheatmap(cor(dataLogFinal, use = "complete.obs"), scale = "none")




#................................................................................
# Exploring 50vs5 => FC = 10 ####
upsFC <- 10

df10 <- data |> 
  dplyr::filter(`Original comparison` == "50Vs5") |> 
  dplyr::select(Species, `Majority protein IDs`, `A1 (raw)`:`B3 (raw)`)
colnames(df10) <- gsub(colnames(df10), pattern = " (raw)", replacement = "", fixed = TRUE)
colnames(df10)

df10[df10 == 0] <- NA


## Assesing missingness ####
data10 <- df10 |> dplyr::filter(Species != "UPS") |> dplyr::select(A1:B3)
# data10 <- df10 |> dplyr::filter(Species == "UPS") |> dplyr::select(A1:B3)
data10 <- as.data.frame(data10)


Biomics::plotBarTI(data = data10, interact = FALSE)
pNA <- Biomics::plotBarNA(data = data10, interact = FALSE)
pNA$dfMuestra
pNA$graficoMuestra
pNA$graficoProteina

pheatmap::pheatmap(mat = data10, show_rownames = FALSE, 
                   cluster_rows = FALSE, cluster_cols = FALSE)

### Filtering ####

data10 <- df10 |> dplyr::select(A1:B3)
data10 <- as.data.frame(data10)

rownames(data10) <- gsub(
  x = sapply(
    strsplit(
      x = df10$`Majority protein IDs`, split = "|", fixed = TRUE),
    "[[", 1),
  pattern = ">",
  replacement = "",
  fixed = TRUE)


Biomics::filterMissing(
  df = data10, 
  dfGrupos = dm, 
  # threshold = c(0, 0.3, 0.5, 0.7, 1), 
  threshold = c(0, 1, 2, 3), 
  mode = "minDataByGroup"
)$tablaFormato

dataFinal <- Biomics::filterMissing(
  df = data10, 
  dfGrupos = dm, 
  threshold = 2, 
  mode = "minDataByGroup"
)$tabla

rownames(dataFinal)[1:2]

dataLogFinal <- log(dataFinal, base = 2)

## Log and visualization  ####
pheatmap::pheatmap(dataFinal, show_rownames = F)
apply(dataLogFinal, 2, min, na.rm = TRUE)
pheatmap::pheatmap(dataLogFinal, show_rownames = F, scale = "row")
pheatmap::pheatmap(cor(dataLogFinal, use = "complete.obs"), scale = "none")


#..............................................................................----
## Common data processing  ####
#..............................................................................----

# Ejecutar uno de los 2 datasets de benchmarking y esta parte será común para todos

#....................#
## Normalization  ####
listaNorm <- Biomics::doNormalization(
  rawData = dataFinal, 
  logData = dataLogFinal,
  listaNorm = c("Mean", "Median", "TI", "VSN", 
                "Quantile", "CyclicLoess", "RLR"))

listaNorm <- lapply(listaNorm, function(df){
  rownames(df) <- rownames(dataFinal)
  df
})

#....................#
## DA analysis #### 
g1 <- "B"
g2 <- "A"

resultDE <- lapply(
  listaNorm, 
  Biomics::doTestT,
  dfGrupos = dm,
  g1 = g1,
  g2 = g2
)


#....................#
## NormScore #### 

resNS <- normScore::normScore(
  normalizedDataList = listaNorm,
  groupData = dm, 
  rawData = dataFinal, 
  refGroup = "B", 
  altGroup = "A", 
  returnDetails = T
)

resNS$finalRanking

pDiag <- normScore::plotNormScoreDiagnostics(
  normalizedDataList = listaNorm,
  groupData = dm,
  rawData = dataFinal,
  refGroup = "B",
  altGroup = "A"
)
# pDiag$item0
# pDiag$item1
# pDiag$item2
# pDiag$item3
# pDiag$item4
# pDiag$item5
# pDiag$item6



#....................#
## Gold Standard #### 

thresholdsROC <- sort(unique(c(
  seq(0, 0.01, by = 0.0001),
  seq(0.01, 0.10, by = 0.001),
  seq(0.10, 1, by = 0.01)
)))


evalResults <- lapply(
  resultDE,
  evalAUC, 
  upsFC  = upsFC
)


### All results together ####
dfRes <- data.frame(
  Normalizations = names(sapply(evalResults, "[[", 1)),
  NormScore = resNS$finalRanking[names(sapply(evalResults, "[[", 1))],
  MAE = sapply(evalResults, "[[", 3), 
  AUC = sapply(evalResults, "[[", 1))


# By organism MAE
tabsByOrg <- sapply(evalResults, "[[", 4, simplify = FALSE) # MAE
upsMAE <- sapply(tabsByOrg, function(df) df[1, 2])
dfRes$upsMAE <- upsMAE[dfRes$Normalizations]
yeastMAE <- sapply(tabsByOrg, function(df) df[2, 2])
dfRes$maeYeast <- yeastMAE[dfRes$Normalizations]

dfRes$RankNS <- unname(rank(resNS$finalRanking[dfRes$Normalizations]))
dfRes$RankMAE <- rank(dfRes$MAE)
dfRes$RankAUC <- rank(-dfRes$AUC)
dfRes$RankMaeYeast <- rank(dfRes$maeYeast)
dfRes$RankupsMAE <- rank(dfRes$upsMAE)


#............................................................................
# Correlations ####

# Kendall tau
dfRank <- dfRes |>  dplyr::select(starts_with("Rank"))
dfCor <- cor(dfRank, method = "kendall")
pheatmap::pheatmap(dfCor)

corrplot::corrplot.mixed(corr = dfCor, upper = "ellipse")



#............................................................................
# Saving ####

listaCompleta <- 
  list(
    dataRaw = dataFinal, 
    normData = listaNorm, 
    dm = dm, 
    resultDEA = resultDE, 
    resultNS = resNS, 
    dfRanks = dfRes, 
    dfCorKendall = dfCor
  )

# saveRDS(
#   listaCompleta,
#   file = "../Data/PXD001819/completeResult25vs12.5_MaxQuant.rds")
saveRDS(
  listaCompleta,
  file = "../Data/PXD001819/completeResult50vs5_MaxQuant.rds")


