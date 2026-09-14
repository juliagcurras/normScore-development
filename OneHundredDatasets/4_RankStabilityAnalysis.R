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

compareNSvsLOO <- function(rankNS, rankLOO, normalization = "Median", maxColor = "darkblue"){
  
  # 1) Comparing with kendall
  kendallNormX <- sapply(seq_along(rankNS), function(i) {
    
    original <- rankNS[[i]]
    original <- original[names(original) != normalization]
    loo <- rankLOO[[i]]
    
    # Posición de cada normalización en el ranking original
    posOriginal <- seq_along(original)
    names(posOriginal) <- names(original)
    
    # Posición de esas mismas normalizaciones en el ranking LOO
    posLOO <- match(names(original), names(loo))
    
    stopifnot(setequal(names(original), names(loo)))
    
    cor(posOriginal, posLOO, method = "kendall")
  })
  
  # 2) Plotting Kendall
  dfKendall <- data.frame(
    dataset = seq_along(kendallNormX),
    normalization = normalization,
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
    original <- original[names(original) != normalization]
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
    scale_x_continuous(breaks = 1:6) +
    scale_y_continuous(breaks = 1:6) +
    scale_fill_gradient(
      low = "white",
      high = maxColor, limits = c(0, 100)
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "bottom") +
    labs(
      x = "Original rank",
      y = paste0("LOO ", normalization, " rank"),
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
  finalRank[names(finalRank)[!(names(finalRank) == "Log")]]
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
normToRemove <- "Median"

## Removing a normalization method + normScore ####
nsListLOO <- lapply(allNorms, function(normToRemove){
    cat("\t* ", normToRemove, "\n")
  normScoreList <- sapply(allFiles, function(i){
    cat("\t* ", i , "\n")
    reducedNorms <- c("Log", allNorms[!(normToRemove == allNorms)])
    output <- readRDS(file = paste0(outDir, i, ".rds"))
    return(
      normScore(
        normMatrixList = output$listaNorm[reducedNorms],
        designMatrix = output$dm,
        dfRaw = output$data,
        corrected = F,
        onlyFinalRank = F,
        onlyDetailRanking = T
      )$detailRanking)
  }, simplify = FALSE, USE.NAMES = TRUE)
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



#...........................................................................####
# Global comparison  ####
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




# CyclicLoess exploration ####

resCyclicloess <- nsListLOO$CyclicLoess

common <- intersect(names(itemOriginal), names(itemLOO))

cor(
  itemOriginal[common],
  itemLOO[common],
  method = "spearman"
)

