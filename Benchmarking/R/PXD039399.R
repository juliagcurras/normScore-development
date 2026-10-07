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

doTestTMod <- function (df, dfGrupos, g1 = "Control", g2 = "Case") 
{
  df <- as.data.frame(df)
  dfGrupos <- as.data.frame(dfGrupos)
  colnames(dfGrupos) <- c("Samples", "Groups")
  if ((!(g1 %in% dfGrupos$Groups)) | (!(g2 %in% dfGrupos$Groups))) {
    return(NULL)
  }
  samplesG1 <- as.character(dfGrupos[dfGrupos$Groups == g1, 
                                     "Samples"])
  samplesG2 <- as.character(dfGrupos[dfGrupos$Groups == g2, 
                                     "Samples"])
  df <- df[, c(samplesG1, samplesG2)]
  df$logFC <- apply(df, 1, function(x) {
    mean(x[samplesG2], na.rm = T) - mean(x[samplesG1], na.rm = T)
  })
  df$AveExpr <- apply(df, 1, function(x) {
    (mean(x[samplesG2], na.rm = T) + mean(x[samplesG1], 
                                          na.rm = T))/2
  })
  df[, c("t", "P.Val")] <- t(apply(df, 1, function(x) {
    evalG1NA <- sum(is.na(x[samplesG1]))
    evalG2NA <- sum(is.na(x[samplesG2]))
    testDone <- try(stats::t.test(x[samplesG2], x[samplesG1], 
                                  alternative = "two.sided", var.equal = TRUE), T)
    if (!(class(testDone) == "try-error")) {
      res <- testDone
      return(c(res$statistic, res$p.value))
    }
    else {
      return(c(NA, NA))
    }
  }))
  df$adj.P.Val <- stats::p.adjust(df$P.Val, method = "BH")
  # df <- df[, c("logFC", "AveExpr", "t", "P.Val", "adj.P.Val")]
  df <- df[, c("logFC", "adj.P.Val")]
  df$Comparison <- paste0(g1, "vs", g2)
  return(df)
}


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



evalAUC <- function(dfDE, comparison = NULL){
  
  # Get initial dataset #
  if (!is.null(comparison) && !(comparison %in% colnames(dfDE))){
    dfDE$Comparison <- comparison
  }
  # table(dfDE$Comparison)
  dfDE$Comparison <- factor(dfDE$Comparison, 
                            levels = c("Y2HvsY1H",  "Y3HvsY1H", "Y4HvsY1H"))
  df <- na.omit(dfDE)
  
  ## B) Regression 
  ### MAE 
  #### Expected values
  df$logFC_Real <- ifelse(is.na(df$Comparison), NA,  
                          ifelse(df$Comparison == "Y2HvsY1H", log(200/100, 2), 
                                 ifelse(df$Comparison == "Y3HvsY1H", log(200/50, 2), 
                                        ifelse(df$Comparison == "Y4HvsY1H", log(200/25, 2),"Error"
                                )
                              )
                            )
                          )
  df$logFC_Real <- as.numeric(df$logFC_Real)
  table(df$logFC_Real)
  
  
  #### MAE global 
  maeG <- mean(abs(df$logFC_Real - df$logFC))
  
  #### MAE by class
  tablaC <- df %>% 
    group_by(Comparison) %>% 
    summarise(MAE = mae(actual = logFC_Real, predicted = logFC)) %>%
    as.data.frame()
  maeGlobal <- maeG
  
  ## Returning 
  finalResult <- list(
    maeGlobal = maeGlobal,
    tablaMaeOrg = tablaC

  )
  
  return(finalResult)
}



#............................................................................----
# Loading data ####
dfRaw <- read.table(file = "../Data/PXD039399/report.pg_matrix.tsv", 
                    header = T, sep = "\t")

table(duplicated(dfRaw$Protein.Group))

## Selecting data ####
data <- dfRaw |>
  dplyr::select(C..Users.dell.Desktop.14.0930Y1H_R1.raw:C..Users.dell.Desktop.14.0930Y4H_R3.raw)

rownames(data) <- dfRaw$Protein.Group

colnames(data) <- gsub(
  pattern = "C..Users.dell.Desktop.14.0930", 
  replacement = "", 
  x = colnames(data))
colnames(data) <- gsub(
  pattern = ".raw", 
  replacement = "", 
  x = colnames(data))

colnames(data)

elementsSamples <- strsplit(x = colnames(data), split = "_", fixed = T)

## Initial group information ####
dm <- data.frame(
  Samples = colnames(data), 
  Groups = sapply(elementsSamples, "[[", 1)
)

table(dm$Groups)

## Log Transformations ####

dfLog <- log(data, base = 2)

head(dfLog)


## Missingness ####
data$Miss <- rowSums(is.na(data))
table(data$Miss == (ncol(data)-1)) # hay filas vacias, fuera
data <- data[which(data$Miss < (ncol(data)-1)), -ncol(data)]
table(rowSums(is.na(data)) == (ncol(data)))


Biomics::plotBarTI(data = data, interact = FALSE, color = "red")
pNA <- Biomics::plotBarNA(data = data, interact = FALSE)
pNA$dfMuestra
pNA$graficoMuestra
pNA$graficoProteina

pheatmap::pheatmap(mat = data, show_rownames = FALSE, 
                   cluster_rows = FALSE, cluster_cols = FALSE)

### Filtering ####
Biomics::filterMissing(
  df = data, 
  dfGrupos = dm, 
  threshold = c(0, 0.3, 0.5, 0.7, 1)
)$tablaFormato

dataFinal <- data

rownames(dataFinal)[1:2]

dataLogFinal <- log(dataFinal, base = 2)

## Log and visualization  ####
apply(dataLogFinal, 2, min, na.rm = TRUE)
pheatmap::pheatmap(dataLogFinal, show_rownames = F, scale = "row", 
                   cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap::pheatmap(cor(dataLogFinal, use = "complete.obs"), scale = "none")


## Only 2 groups ####
# dataFinal <- dataFinal |> select(Y1H_R1:Y1H_R3, Y2H_R1:Y2H_R3)
# dataLogFinal <- log(dataFinal, base = 2)

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


#..............................................................................----
## Differential abundance analysis  ####

g2 <- "Y1H"
grupos1 <- c("Y2H", "Y3H", "Y4H")

resultDE <- lapply(listaNorm, function(df){
  resDE <- lapply(
    grupos1, 
    doTestTMod,
    df = df,
    dfGrupos = dm,
    # g1 = g1,
    g2 = g2
  )
  do.call(rbind, resDE)
})



#..............................................................................----
## NormScore #### 

## All groups ####

resNS <- normScore::normScore(
  normalizedDataList = listaNorm,
  # groupData = dmSub, 
  groupData = dm,
  rawData = dataFinal, 
  returnDetails = T
)

resNS$finalRanking

pDiag <- normScore::plotNormScoreDiagnostics(
  normalizedDataList = listaNorm,
  groupData = dmSub,
  rawData = dataFinal,
  refGroup = "Y1H", 
  altGroup = "Y2H"
)
pDiag$item0
pDiag$item1
pDiag$item2
pDiag$item3
pDiag$item4
pDiag$item5
pDiag$item6

## By groups ####
nsByNorm <- lapply(
  unique(dm$Groups), 
  function(grupo){
    dmSub <- dm |> filter(Groups %in% grupo)
    dataFinalSub <- dataFinal[, dmSub$Samples]
    normListSub <- lapply(listaNorm, function(df) df[, dmSub$Samples])
    
    normScore::normScore(
      normalizedDataList = normListSub,
      groupData = dmSub, 
      rawData = dataFinalSub
    )
  })

names(nsByNorm) <- unique(dm$Groups)



#........................................................................####
## Gold Standard #### 

evalResults <- lapply(
  resultDE,
  evalAUC,
  comparison = NULL
)

names(evalResults$Log)

### All results together ####
dfRes <- data.frame(
  Normalizations = names(sapply(evalResults, "[[", 1)),
  normScore = resNS$finalRanking[names(sapply(evalResults, "[[", 1))],
  MAE = sapply(evalResults, "[[", 1)) 


dfRes$NSY1vsY2 <- (nsByNorm$Y1H$finalRanking[dfRes$Normalizations] + nsByNorm$Y2H$finalRanking[dfRes$Normalizations])/2
dfRes$NSY1vsY3 <- (nsByNorm$Y1H$finalRanking[dfRes$Normalizations] + nsByNorm$Y3H$finalRanking[dfRes$Normalizations])/2
dfRes$NSY1vsY3 <- (nsByNorm$Y1H$finalRanking[dfRes$Normalizations] + nsByNorm$Y3H$finalRanking[dfRes$Normalizations])/2
dfRes$NSY1vsY4 <- (nsByNorm$Y1H$finalRanking[dfRes$Normalizations] + nsByNorm$Y4H$finalRanking[dfRes$Normalizations])/2


# By organism MAE
tabsByOrg <- sapply(evalResults, "[[", 2, simplify = FALSE) # MAE

compHMAE <- sapply(tabsByOrg, function(df) df[1, 2])
dfRes$Y2Hmae <- compHMAE[dfRes$Normalizations]

compHMAE <- sapply(tabsByOrg, function(df) df[2, 2])
dfRes$Y3Hmae <- compHMAE[dfRes$Normalizations]

compHMAE <- sapply(tabsByOrg, function(df) df[3, 2])
dfRes$Y4Hmae <- compHMAE[dfRes$Normalizations]


dfRes$RankNS <- unname(rank(resNS$finalRanking[dfRes$Normalizations]))
dfRes$RankMAE <- rank(dfRes$MAE)

dfRes$RankNSNSY1vsY2 <- rank(dfRes$NSY1vsY2)
dfRes$RankNSNSY1vsY3 <- rank(dfRes$NSY1vsY3)
dfRes$RankNSNSY1vsY3 <- rank(dfRes$NSY1vsY3)
dfRes$RankNSNSY1vsY4 <- rank(dfRes$NSY1vsY4)

dfRes$RankY2Hmae <- rank(dfRes$Y2Hmae)
dfRes$RankY3Hmae <- rank(dfRes$Y3Hmae)
dfRes$RankY4Hmae <- rank(dfRes$Y4Hmae)

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
    resultNSGlobal = resNS, 
    resultNSBygroup = nsByNorm, 
    dfRanks = dfRes, 
    dfCorKendall = dfCor
  )

saveRDS(
  listaCompleta,
    file = "../Data/PXD039399//completeResult.rds")




