##############################################################################_

############       Non-redundancy analysis - normScore       #################-

###############################################################################_


# 2026/09/14 - Julia G Currás
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets")


#...........................................................................####
# Libraries & global objectss ####
library(ggplot2)
library(ggpubr)

outDir <-  "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/ProcessedDatasets/"


#...........................................................................####
# Functions ####
source(file = "../R/scoreFunction.R", encoding = "UTF-8")

compareNSvsLOO <- function(
    rankNS, 
    rankLOO, 
    normalization = "Median",
    nNormExcluded = 1,
    maxColor = "darkblue"){
  
  
  # 1) Comparing with kendall
  kendallNormX <- sapply(seq_along(rankNS), function(i) {
    # cat("\t", i)
    
    original <- rankNS[[i]]
    original <- original[!(names(original) %in% normalization)]
    loo <- rankLOO[[i]]
    
    # Posición de cada normalización en el ranking original
    posOriginal <- seq_along(original)
    names(posOriginal) <- names(original)
    
    # Posición de esas mismas normalizaciones en el ranking LOO
    posLOO <- match(names(original), names(loo))
    
    stopifnot(setequal(names(original), names(loo)))
    
    cor(posOriginal, posLOO, method = "kendall")
  })
  
  # 1.1) Just adjusting parameters
  if (nNormExcluded == 1){
    normLabel <- normalization
    nComparisons <- 7
    }
  if (nNormExcluded == 2){
    normLabel <- paste0(normalization, collapse = ", ")
    nComparisons <- 21
  } 
  
  # 2) Plotting Kendall
  dfKendall <- data.frame(
    dataset = seq_along(kendallNormX),
    normalization = normLabel,
    kendall = kendallNormX
  )
  
  pKendall <- ggplot(dfKendall, aes(x = normalization, y = kendall)) +
    geom_boxplot(outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.08, alpha = 0.4) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Kendall's tau"
    ) + ylim(c(-1, 1))
  
  # 3) Heatmaps
  dfRanks <- do.call(rbind, lapply(seq_along(rankNS), function(i) {
    
    original <- rankNS[[i]]
    original <- original[!(names(original) %in% normalization)]
    loo <- rankLOO[[i]]
    
    # Posición de cada normalización en el ranking original
    posOriginal <- seq_along(original)
    names(posOriginal) <- names(original)
    
    # Posición de esas mismas normalizaciones en el ranking LOO
    posLOO <- match(names(original), names(loo))
    
    
    data.frame(
      dataset = i,
      normalization = names(original),
      rankOriginal = posOriginal,
      rankLOO = posLOO
    )
  }))
  
  dfHeatmap <- as.data.frame(
    table(dfRanks$rankOriginal, dfRanks$rankLOO)
  )
  
  names(dfHeatmap) <- c("rankOriginal", "rankLOO", "frequency")
  
  dfHeatmap$rankOriginal <- as.numeric(as.character(dfHeatmap$rankOriginal))
  dfHeatmap$rankLOO <- as.numeric(as.character(dfHeatmap$rankLOO))
  
  pHeatmap <- ggplot(dfHeatmap, aes(x = rankOriginal, y = rankLOO, fill = frequency)) +
    geom_tile() +
    geom_text(aes(label = frequency)) +
    scale_x_continuous(breaks = 1:nComparisons) +
    scale_y_continuous(breaks = 1:nComparisons) +
    scale_fill_gradient(
      low = "white",
      high = maxColor, limits = c(0, 100)
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "bottom") +
    labs(
      x = "Original rank",
      y = paste0("LOO ", normLabel, " rank"),
      fill = "Frequency"
    )
  
  
  # 4) Outputs
  return(
    list(
      kendalDF = kendallNormX, 
      kendalPloDF = dfKendall,
      pKendall = pKendall, 
      pHeatmap = pHeatmap
    )
  )
  
}




#...........................................................................####
# Loading data ####
dataGS <- readxl::read_excel(path = "_Assessment.xlsx", sheet = "Main", col_names = T)[-1,]
allFiles <- dataGS[, "ID", drop = T]

# normScore results over 100 dataset - 8 normalization methods
resNS <- readRDS(file = "AssessmentFiles/normScore_dataGS_All_Item0_x4_Item2_01.rds")
rankingNS <- lapply(resNS, function(df){
  scoreDF_norm <- df |> dplyr::arrange(TotalCorrected)
  finalRank <- stats::setNames(scoreDF_norm$TotalCorrected, rownames(scoreDF_norm))
  # finalRank[names(finalRank)[!(names(finalRank) == "Log")]]
})


#...........................................................................####
# LOO for normalization methods  ####
allNorms <- c(
  "Median", 
  "Mean", 
  "TI",
  "Quantile", 
  "CyclicLoess",
  "RLR",
  "VSN" 
)
normToRemove <- "CyclicLoess"

## Removing a normalization method + normScore ####
nsListLOO <- lapply(allNorms, function(normToRemove){
  cat(normToRemove, "      .................................................\n")
  normScoreList <- sapply(allFiles, function(i){
    cat("\t* ", i , "\n")
    reducedNorms <- c("Log", allNorms[!(normToRemove == allNorms)])
    output <- readRDS(file = paste0(outDir, i, ".rds"))
    return(
      abc <- normScore(
        normMatrixList = output$listaNorm[reducedNorms],
        designMatrix = output$dm,
        dfRaw = output$data,
        corrected = F,
        onlyFinalRank = F,
        onlyDetailRanking = T
      )$detailRanking
      )
  }, simplify = FALSE, USE.NAMES = TRUE)
  cat("\n\n")
  normScoreList
})
names(nsListLOO) <- allNorms

saveRDS(object = nsListLOO, file = "AssessmentFiles/rankStability_LOO_normalization_Item0_x4_Item2_01.rds")

## Just ranking by dataset ####
rankingNsLoo <- lapply(nsListLOO, function(resNsLOO){
  lapply(resNsLOO, function(df){
    scoreDF_norm <- df |> dplyr::arrange(TotalCorrected)
    finalRank <- stats::setNames(scoreDF_norm$TotalCorrected, rownames(scoreDF_norm))
  })
})



## Global comparison  ####
# abc <- compareNSvsLOO(rankingNS, rankLOO = rankingNsLoo$Median, normalization = "Median")
res <- lapply(names(rankingNsLoo), function(i) 
  compareNSvsLOO(rankNS = rankingNS, rankLOO = rankingNsLoo[[i]], normalization = i))
names(res) <- names(rankingNsLoo)
ggpubr::ggarrange(plotlist = sapply(res, "[[", 3), nrow = 2, ncol = 4)
ggpubr::ggarrange(plotlist = sapply(res, "[[", 4), nrow = 2, ncol = 4, common.legend = TRUE)


# Testing...
dfKendall <- do.call(rbind, sapply(res, "[[", 2, USE.NAMES = T, simplify = F))

# Test first
dfTest <- dfKendall |>
  rstatix::pairwise_wilcox_test(
    kendall ~ normalization,
    paired = TRUE,
    p.adjust.method = "holm"
  )

# Plot without test

ggplot(dfKendall, aes(x = normalization, y = kendall)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Kendall's tau"
  )

# plot with test result

dfTest <- dfTest %>%
  rstatix::add_xy_position(x = "normalization")
dfTest$y.position <- 1+seq(0.1, 1.1, 0.05)

ggplot(dfKendall, aes(x = normalization, y = kendall)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25) +
  stat_pvalue_manual(
    dfTest,
    label = "p.adj.signif",
    hide.ns = TRUE
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Kendall's tau"
  )



#...........................................................................####
# LOO for pairs of normalization methods  ####
dfPairs <- as.data.frame(t(combn(allNorms, 2)))


## Removing a normalization method + normScore ####
# nsListLOO <- lapply(dfPairs, function(normToRemove){
listPairsLOO <- apply(dfPairs, 1, function(normToRemove){
  cat(normToRemove[1]," - ", normToRemove[2], "      .................................................\n")
  reducedNorms <- c("Log", allNorms[!(allNorms %in% normToRemove)])
  normScoreList <- sapply(allFiles, function(i){
    cat("\t* ", i , "\n")
    output <- readRDS(file = paste0(outDir, i, ".rds"))
    return(
      normScore(
        normMatrixList = output$listaNorm[reducedNorms],
        designMatrix = output$dm,
        dfRaw = output$data,
        corrected = F,
        onlyFinalRank = F,
        onlyDetailRanking = T
      )$detailRanking
    )
  }, simplify = FALSE, USE.NAMES = TRUE)
  cat("\n\n")
  normScoreList
})
dfPairs$Name <- paste0(dfPairs$V1, "_", dfPairs$V2)
names(listPairsLOO) <- dfPairs$Name

saveRDS(object = listPairsLOO, file = "AssessmentFiles/rankStability_LOO_Pairs_normalization_Item0_x4_Item2_01.rds")
listPairsLOO <- readRDS(file = "AssessmentFiles/rankStability_LOO_Pairs_normalization_Item0_x4_Item2_01.rds")

## Just ranking by dataset ####
rankingNsPairsLoo <- lapply(listPairsLOO, function(resNsLOO){
  lapply(resNsLOO, function(df){
    scoreDF_norm <- df |> dplyr::arrange(TotalCorrected)
    finalRank <- stats::setNames(scoreDF_norm$TotalCorrected, rownames(scoreDF_norm))
  })
})




## Global comparison  ####
dfPairs <- as.data.frame(dfPairs)
resPairs <- lapply(dfPairs$Name, function(i){
  cat(i, "\n")
  compareNSvsLOO(
    rankNS = rankingNS, 
    rankLOO = rankingNsPairsLoo[[i]], 
    normalization = unname(unlist(as.vector(dfPairs[dfPairs$Name == i, 1:2]))), 
    nNormExcluded = 2
    )
  }
)
names(resPairs) <- names(rankingNsPairsLoo)
ggpubr::ggarrange(
  plotlist = sapply(resPairs, "[[", 3), nrow = 7, ncol = 3)
ggpubr::ggarrange(
  plotlist = sapply(resPairs, "[[", 4), nrow = 7, ncol = 3, 
  common.legend = TRUE)


# Testing...
dfKendall <- do.call(rbind, sapply(resPairs, "[[", 2, USE.NAMES = T, simplify = F))

# Test first
dfTest <- dfKendall |>
  rstatix::pairwise_wilcox_test(
    kendall ~ normalization,
    paired = TRUE,
    p.adjust.method = "holm"
  )

# Plot without test

ggplot(dfKendall, aes(x = normalization, y = kendall)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  geom_jitter(width = 0.12, alpha = 0.25) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Kendall's tau"
  )

# plot with test result

dfTest <- dfTest %>%
  rstatix::add_xy_position(x = "normalization")
dfTest$y.position <- 1+seq(0.1, 2.105, 0.009569)

ggplot(dfKendall, aes(x = normalization, y = kendall)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25) +
  stat_pvalue_manual(
    dfTest,
    label = "p.adj.signif",
    hide.ns = TRUE
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Kendall's tau"
  )
