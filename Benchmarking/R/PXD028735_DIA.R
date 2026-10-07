###############################################################################\
##################             DIA BENCHMARKING             ###################\
###############################################################################\

# Set directory
setwd("C:/Users/julia/Documents/GitHub/normScore-development/Benchmarking/R")

# Libraries ####
library(dplyr)
library(ggplot2)


#..............................................................................
# Loading data ####
## Raw data ####
dfRaw <- read.table(file = "../Data/PXD028735/report.pg_matrix.tsv", 
                    header = T, sep = "\t")
dfRaw$Specie <- sapply(strsplit(x = dfRaw$Protein.Names, split = "_", fixed = TRUE), "[[", 2)

table(dfRaw$Specie)

dfRaw$WeirdSpecie <- grepl(pattern = ";", x = dfRaw$Specie)

dfRawFilter <- dfRaw |> dplyr::filter(!WeirdSpecie)


## Selecting data ####
data <- dfRawFilter |>
  dplyr::select(
    C..Users.julia.Documents.Doctorado.benchmarking.DIA.LFQ_TTOF6600_SWATH_Condition_A_Sample_Alpha_01.wiff:C..Users.julia.Documents.Doctorado.benchmarking.DIA.LFQ_TTOF6600_SWATH_Condition_B_Sample_Gamma_03.wiff)

rownames(data) <- dfRawFilter$Protein.Names

colnames(data) <- gsub(
  pattern = "C..Users.julia.Documents.Doctorado.benchmarking.DIA.LFQ_TTOF6600_SWATH_Condition_", 
  replacement = "", 
  x = colnames(data))
colnames(data) <- gsub(
  pattern = ".wiff", 
  replacement = "", 
  x = colnames(data))
colnames(data) <- gsub(
  pattern = "_Sample", 
  replacement = "", 
  x = colnames(data))

colnames(data)

elementsSamples <- strsplit(x = colnames(data), split = "_", fixed = T)

## Initial group information ####
dfGroup <- data.frame(
  Samples = colnames(data), 
  Groups = sapply(elementsSamples, "[[", 1), 
  Replicates = paste0(
    sapply(elementsSamples, "[[", 1),  "_",
    sapply(elementsSamples, "[[", 2)
  )
)

dfGroup

## Log Transformations ####

dfLog <- log(data, base = 2)

head(dfLog)

#..............................................................................
# Initial Processing steps ####

## Distribution ####
Biomics::plotBarTI(data, color = "blue", interact = FALSE)$grafico
Biomics::plotBoxMulti(base = dfLog, varResumen = colnames(dfLog), interact = FALSE)$grafico

## Missing data ####
# Empty rows?
data$Miss <- rowSums(is.na(data))
table(data$Miss == (ncol(data)-1)) # hay filas vacias, fuera
data <- data[which(data$Miss < (ncol(data)-1)), -ncol(data)]
table(rowSums(is.na(data)) == (ncol(data)))

pNA <- Biomics::plotBarNA(data = data, interact = FALSE)
pNA$graficoMuestra
pNA$graficoProteina

dfLog <- log(data, base = 2)
# sum(is.nan(dfLog))
pheatmap::pheatmap(mat = as.matrix(dfLog), cluster_cols = FALSE, cluster_rows = FALSE)


#..............................................................................
# Join replicates ####
dataFilt <- Biomics::doJoinReplicates(
  dfQ = data, dm = dfGroup, sample_col = "Replicates", rep_col = "Samples"
)
dfLogFilt <- log(dataFilt, base = 2)

dfGroupFilt <- data.frame(
  Samples = colnames(dataFilt), 
  Groups = sapply(strsplit(x = colnames(dataFilt), split = "_", fixed = T), "[[", 1)
)


#..............................................................................
# Second Processing steps #### 
Biomics::plotBoxMulti(base = dataFilt, varResumen = colnames(dataFilt), interact = FALSE)$grafico
Biomics::plotBoxMulti(base = dfLogFilt, varResumen = colnames(dfLogFilt), interact = FALSE)$grafico

pNA <- Biomics::plotBarNA(data = dataFilt, interact = FALSE)
pNA$graficoMuestra
pNA$graficoProteina

## Filtering ####
Biomics::filterMissing(dataFilt, dfGroupFilt, threshold = c(0, 0.3, 0.6))$tablaFormato

dataFinal <- Biomics::filterMissing(dataFilt, dfGroupFilt, threshold = 0)$tabla
dataLogFinal <- log(dataFinal, base = 2)
Biomics::plotBoxMulti(base = dataLogFinal, varResumen = colnames(dataLogFinal), interact = F)$grafico

## Normalization ####
listaNorm <- Biomics::doNormalization(
  rawData = dataFinal, logData = dataLogFinal,
  listaNorm = c("Mean", "Median", "TI", "VSN", "Quantile", "CyclicLoess", "RLR"))

listaNorm <- lapply(listaNorm, function(df){
  rownames(df) <- rownames(dataFinal)
  df
})



#..............................................................................
# DAA #### 
g1 <- "A"
g2 <- "B"


# DE analysis ####
resultDE <- lapply(
  listaNorm, 
  Biomics::doTestT,
  dfGrupos = dfGroupFilt,
  g1 = g1,
  g2 = g2
)



#..............................................................................
# NormScore #### 

resNS <- normScore::normScore(
  normalizedDataList = listaNorm,
  groupData = dfGroupFilt, 
  rawData = dataFinal, 
  refGroup = "A", 
  altGroup = "B", 
  returnDetails = T
)

resNS$finalRanking
resNS$detailRanking

diagnosticPlots <- normScore::plotNormScoreDiagnostics(
  normalizedDataList = listaNorm,
  groupData = dfGroupFilt, 
  rawData = dataFinal, 
  refGroup = "A", 
  altGroup = "B"
)
diagnosticPlots$item0
diagnosticPlots$item1
diagnosticPlots$item2
diagnosticPlots$item3
diagnosticPlots$item4
diagnosticPlots$item5
diagnosticPlots$item6


#..............................................................................
# Gold Standard #### 

## Functions ####
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



evalAUC <- function(dfDE, dfIDs, minLogFC = 5){
  
  # Get initial dataset #
  dfDE$Protein.Ids <- rownames(dfDE)
  df <- base::merge(dfDE, dfIDs, by = "Protein.Ids")
  df$Organism <- factor(df$Organism, 
                        levels = c("HUMAN", "YEAST", "ECOLI"))
  
  ## Filtering 
  df <- df %>%
    filter(logFC < 5) %>%
    filter(logFC > -5) %>%
    filter(AveExpr <= Biostatech::getOutlier(df$AveExpr)$max_sinout)
  
  
  ## A) AUC and ROC curves 
  
  ## AUC sobre todos los puntos de corte relevantes
  # Estimate metrics #
  df$truth <- ifelse(is.na(df$Organism), NA, 
                     ifelse(df$Organism == "HUMAN", 0, 1)) # Proteína de humano non DE (0)
  
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
                          ifelse(df$Organism == "ECOLI", 2, 
                                 ifelse(df$Organism == "YEAST", -1,
                                        ifelse(df$Organism == "HUMAN", 0, "Error"
                                        )
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


## Analyzing ####
dfIDs <- dfRawFilter |> dplyr::select(Protein.Names, Specie) |> 
  dplyr::rename(Organism = Specie, Protein.Ids = Protein.Names)

thresholdsROC <- sort(unique(c(
  seq(0, 0.01, by = 0.0001),
  seq(0.01, 0.10, by = 0.001),
  seq(0.10, 1, by = 0.01)
)))

## EXECUTION ####

evalResults <- lapply(
  resultDE,
  evalAUC,
  dfIDs = dfIDs
)

dfRes <- data.frame(
  Normalizations = names(sapply(evalResults, "[[", 1)),
  NormScore = resNS$finalRanking[names(sapply(evalResults, "[[", 1))],
  MAE = sapply(evalResults, "[[", 3), 
  AUC = sapply(evalResults, "[[", 1)
  )

# By organism MAE
tabsByOrg <- sapply(evalResults, "[[", 4, simplify = FALSE) # MAE

humanMAE <- sapply(tabsByOrg, function(df) df[1, 2])
sort(humanMAE)
dfRes$maeHUMAN <- humanMAE[dfRes$Normalizations]

yeastMAE <- sapply(tabsByOrg, function(df) df[2, 2])
sort(yeastMAE)
dfRes$maeYeast <- yeastMAE[dfRes$Normalizations]

ecoliMae <- sapply(tabsByOrg, function(df) df[3, 2])
sort(ecoliMae)
dfRes$ecoliMae <- ecoliMae[dfRes$Normalizations]


# Ranks
dfRes$RankNS <-  rank(dfRes$NormScore, ties.method = "average")
dfRes$RankMAE  <- rank(dfRes$MAE, ties.method = "average")
dfRes$RankAUC  <- rank(-dfRes$AUC, ties.method = "average")  # mayor AUC = mejor
dfRes$RankHuman <- rank(dfRes$maeHUMAN)
dfRes$RankYeast <- rank(dfRes$maeYeast)
dfRes$RankEcoli <- rank(dfRes$ecoliMae)

dfRes

#............................................................................
# Correlations ####
# Kendal
dfRank <- dfRes |>  dplyr::select(starts_with("Rank"))
dfCor <- cor(dfRank, method = "kendall")
pheatmap::pheatmap(dfCor)

corrplot::corrplot.mixed(corr = dfCor, upper = "ellipse")



#..............................................................................
# SAVING #### 

listaCompleta <- list(
  dataRaw = dataFinal, 
  dataLog = dataLogFinal, 
  groupData = dfGroupFilt, 
  normList = listaNorm, 
  resultsDAA = resultDE, 
  gsMetrics = evalResults, 
  resNS = resNS, 
  allMetrics = dfRes, 
  correlationResults = dfCor, 
  diagnosticPlots = diagnosticPlots
)


saveRDS(
  listaCompleta, 
  file = "../Data/PXD028735/PXD028735_DIA_Complete_Analysis.rds"
)


