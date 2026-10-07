###############################################################################\
##################             DDA BENCHMARKING             ###################\
###############################################################################\

# Set directory
setwd("C:/Users/julia/Documents/GitHub/normScore-development/Benchmarking/R")
rm(list=ls())
graphics.off()

# Libraries ####
library(dplyr)
library(ggplot2)

#.............................................................................
# Initial processing steps####

## Loading data ####
dfRaw <- read.table(file = "../Data/PXD028735/proteinGroups.txt", 
                    header = T, sep = "\t")
dfRaw <- dfRaw |> 
  dplyr::filter(!(Fasta.headers %in% c(";", ";;", "")))

dfRaw$namesProtein <- Biomics::getProteinName(protNames = dfRaw$Fasta.headers)


## Selecting interesting information ####
# Selecting interesting cols #
df0 <- dfRaw %>% 
  # dplyr::filter(Number.of.proteins <10) |>
  dplyr::select(
    Protein.IDs, 
    starts_with("Intensity"), #"LFQ.Intensity
    Reverse, Only.identified.by.site,
    Potential.contaminant, 
    Taxonomy.names, 
    -Intensity)

# Set 0 to NA #
df0[df0 == 0] <- NA

# Filters #
dfFilter <- df0 %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

## Annotation data frame ####
table(dfFilter$Taxonomy.names)
dfFilter$Specie <- sapply(strsplit(x = dfFilter$Taxonomy.names, ";", fixed = TRUE), "[[", 1)
table(dfFilter$Specie)

dfAnnotation <- dfFilter |>  
  dplyr::select(Protein.IDs, Specie) |> 
  dplyr::mutate(
    Organism = if_else(Specie == "Homo sapiens", "HUMAN",
                       if_else(Specie == "Escherichia coli K-12", "ECOLI", "YEAST"))
  )
table(dfAnnotation$Organism)
dfFilter$Taxonomy.names <- NULL
dfFilter$Specie <- NULL

# Adjusting rownames #
rownames(dfFilter) <- dfFilter$Protein.IDs
dfFilter$Protein.IDs <- NULL
table(table(rownames(dfFilter)>1))

colnames(dfFilter) <- gsub(x= colnames(dfFilter), pattern = "Intensity.", replacement = "")
colnames(dfFilter) <- gsub(x= colnames(dfFilter), pattern = "Sample_", replacement = "")
colnames(dfFilter)

## Log transforation ####
dfLogFilter <- log(dfFilter,  base = 2)


## Distribution ####
Biomics::plotBarTI(dfFilter, color = "blue", interact = FALSE)$grafico
Biomics::plotBoxMulti(base = dfLogFilter, varResumen = colnames(dfLogFilter), interact = FALSE)$grafico
pheatmap::pheatmap(
  mat = as.matrix(dfLogFilter), 
  cluster_cols = TRUE, cluster_rows = FALSE, show_rownames = FALSE, scale = "row")
pheatmap::pheatmap(cor(dfLogFilter, use = "complete.obs"),show_rownames = TRUE, scale = "none")


pNA <- Biomics::plotBarNA(data = dfFilter, interact = FALSE)
pNA$graficoMuestra
pNA$graficoProteina

## Join replicates ####
elementsSamples <- strsplit(x = colnames(dfFilter), split = "_", fixed = TRUE)
dfReplicates <- data.frame(
  Samples = colnames(dfFilter), 
  Groups = sapply(elementsSamples, "[[", 1), 
  Replicates = paste0(
    sapply(elementsSamples, "[[", 1),  "_",
    sapply(elementsSamples, "[[", 2)
  )
)
dfReplicates

# Dealing with technical replicates for the sample sample #
dataFilter <- Biomics::doJoinReplicates(
  dfQ = dfFilter, dm = dfReplicates,
  sample_col = "Replicates", 
  rep_col = "Samples"
  )
dataLogFilter <- log(dataFilter, 2)

## Final design matrix ####
dm <- data.frame(
  Samples = colnames(dataFilter),
  Groups = sapply(strsplit(x = colnames(dataFilter), split = "_", fixed = TRUE), "[[", 1)
)
dm




#..............................................................................
# Second Processing steps #### 
Biomics::plotBarTI(data = dataFilter, interact = FALSE)$grafico
Biomics::plotBoxMulti(base = dataLogFilter, varResumen = colnames(dataLogFilter), interact = FALSE)$grafico

## Missing data ####
# Empty rows?
dataFilter$Miss <- rowSums(is.na(dataFilter))
table(dataFilter$Miss == (ncol(dataFilter)-1)) # hay filas vacias, fuera
data <- dataFilter[which(dataFilter$Miss < (ncol(dataFilter)-1)), -ncol(dataFilter)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(
  df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
dataFinal <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLogFinal <- log(dataFinal, base = 2)
Biomics::plotBoxMulti(
  base = dataLogFinal, 
  varResumen = colnames(dataLogFinal), interact = F)$grafico


## Log and visualization  ####
pheatmap::pheatmap(dataFinal, show_rownames = F)
apply(dataLogFinal, 2, min)
min(dataLogFinal)
pheatmap::pheatmap(dataLogFinal, show_rownames = F, scale = "row")
pheatmap::pheatmap(cor(dataLogFinal), scale = "none")


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



#..............................................................................
# DAA #### 
g1 <- "A"
g2 <- "B"


# DE analysis ####
resultDE <- lapply(
  listaNorm, 
  Biomics::doTestT,
  dfGrupos = dm,
  g1 = g1,
  g2 = g2
)


#..............................................................................
# NormScore #### 

resNS <- normScore::normScore(
  normalizedDataList = listaNorm,
  groupData = dm, 
  rawData = dataFinal, 
  refGroup = "A", 
  altGroup = "B", 
  returnDetails = T
)

resNS$finalRanking
resNS$detailRanking

diagnosticPlots <- normScore::plotNormScoreDiagnostics(
  normalizedDataList = listaNorm,
  groupData = dm,
  rawData = dataFinal,
  refGroup = "A",
  altGroup = "B"
)




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
dfIDs <- dfAnnotation |> dplyr::select(Protein.IDs, Organism) |> 
  dplyr::rename(Protein.Ids = Protein.IDs)

thresholdsROC <- sort(unique(c(
  seq(0, 0.01, by = 0.0001),
  seq(0.01, 0.10, by = 0.001),
  seq(0.10, 1, by = 0.01)
)))

## EXECUTION ####

### with log data ####
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
  groupData = dm, 
  annotationData = dfAnnotation,
  normList = listaNorm, 
  resultsDAA = resultDE, 
  gsMetrics = evalResults, 
  resNS = resNS, 
  allMetrics = dfRes, 
  correlationResults = dfCor
  # diagnosticPlots = diagnosticPlots
)


saveRDS(
  listaCompleta, 
  file = "../Data/PXD028735/PXD028735_DDA_Complete_Analysis.rds"
)

