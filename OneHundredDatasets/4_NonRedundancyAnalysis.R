##############################################################################_

############       Non-redundancy analysis - normScore       #################-

###############################################################################_


# 2026/09/14 - Julia G Currás
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets")


#...........................................................................####
# Libraries & functions ####
library(ggplot2)
library(ggpubr)


# Function to rank 
recalculateNS <- function(scoreDF_norm, itemToRemove){
  
  # Selecting cols
  allItems <- 1:6
  reducedItems <- allItems[-itemToRemove]
  
  # Obtained original item 0 
  totalCorrected <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "TotalCorrected"]
  totalVal <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "Total"]
  item0 <- totalVal/totalCorrected
  scoreDF_norm <- scoreDF_norm[,reducedItems]
  
  # Now, recalculate final score
  scores_matrix <- t(scoreDF_norm)
  scoreDF_norm$Total <- rowSums(scoreDF_norm)
  scoreDF_norm$TotalCorrected <- scoreDF_norm$Total
  scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "TotalCorrected"] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "TotalCorrected"]*item0
  
  # Sort for final ranking
  scoreDF_norm <- scoreDF_norm |> dplyr::arrange(TotalCorrected)
  finalRank <- stats::setNames(scoreDF_norm$TotalCorrected, rownames(scoreDF_norm))
  
  return(finalRank)
}

# Function to compare 
compareNSvsLOO <- function(rankNS, rankLOO, item = "Item 1", maxColor = "darkred"){
  
  # 1) Comparing with kendall
  kendallItemX <- sapply(seq_along(rankNS), function(i) {
    
    original <- rankNS[[i]]
    loo <- rankLOO[[i]]
    
    # Posición de cada normalización en el ranking original
    posOriginal <- seq_along(original)
    names(posOriginal) <- names(original)
    
    # Posición de esas mismas normalizaciones en el ranking LOO
    posLOO <- match(names(original), names(loo))
    
    cor(posOriginal, posLOO, method = "kendall")
  })
  
  # 2) Plotting Kendall
  dfKendall <- data.frame(
    dataset = seq_along(kendallItemX),
    item = item,
    kendall = kendallItemX
  )
  
  pKendall <- ggplot(dfKendall, aes(x = item, y = kendall)) +
    geom_boxplot(outlier.shape = NA, width = 0.4) +
    geom_jitter(width = 0.08, alpha = 0.4) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Kendall's tau"
    ) + ylim(c(0, 1))
  
  # 3) Heatmaps
  dfRanks <- do.call(rbind, lapply(seq_along(rankNS), function(i) {
    
    original <- rankNS[[i]]
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
    scale_x_continuous(breaks = 1:8) +
    scale_y_continuous(breaks = 1:8) +
    scale_fill_gradient(
      low = "white",
      high = maxColor, limits = c(0, 100)
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "bottom") +
    labs(
      x = "Original rank",
      y = paste0("LOO ", item, " rank"),
      fill = "Frequency"
    )
  
  
  # 4) Outputs
  return(
    list(
      kendalDF = kendallItemX, 
      kendalPloDF = dfKendall,
      pKendall = pKendall, 
      pHeatmap = pHeatmap
    )
  )
  
}



#...........................................................................####
# Loading data ####

# normScore results over 100 dataset
resNS <- readRDS(file = "AssessmentFiles/normScore_dataGS_All_Item0_x4_Item2_01.rds")

# what will happen if we reset item2 to 0-1?
resNsNoWeight <- lapply(resNS, function(scoreDF_norm){
  scoreDF_norm[,2] <- scoreDF_norm[,2]/0.1
  
  # Obtained original item 0 
  totalCorrected <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "TotalCorrected"]
  totalVal <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "Total"]
  item0 <- totalVal/totalCorrected
  
  # Now, recalculate final score
  scores_matrix <- t(scoreDF_norm)
  scoreDF_norm$Total <- rowSums(scoreDF_norm)
  scoreDF_norm$TotalCorrected <- scoreDF_norm$Total
  scoreDF_norm[which(rownames(scoreDF_norm) == "Log"), "TotalCorrected"] <- scoreDF_norm[which(rownames(scoreDF_norm) == "Log"),  "TotalCorrected"]*item0
  
  # Sort for final ranking
  scoreDF_norm <- scoreDF_norm |> dplyr::arrange(TotalCorrected)
  scoreDF_norm
})




#...........................................................................####
# Correlation item by item ####

# No change between item 2 weithged and non-weighted because the weight does not
# affect the ranking, only the magnitude

## Median matrix ####
# Wit item 2*0.1
corListWeight <- lapply(resNS, function(df){
  df <- df[,1:6]
  cor(df, method = "spearman")
})

corArray <- simplify2array(corListWeight)
apply(corArray, c(1, 2), median, na.rm = TRUE)

# With item 2 (without weights)
corListNoWeight <- lapply(resNsNoWeight, function(df){
  df <- df[,1:6]
  cor(df, method = "spearman")
})

corArrayNoW <- simplify2array(corListNoWeight)
apply(corArrayNoW, c(1, 2), median, na.rm = TRUE)



## Boxplots ####
corList <- corListNoWeight
corList <- corListWeight

# Long format by each matrix
corLong <- do.call(rbind, lapply(seq_along(corList), function(i) {
  
  mat <- as.matrix(corList[[i]])
  
  # Only unique pairs
  ind <- which(upper.tri(mat), arr.ind = TRUE)
  
  # All pairs
  ind <- which(row(mat) != col(mat), arr.ind = TRUE)
  
  data.frame(
    dataset = i,
    var1 = rownames(mat)[ind[, 1]],
    var2 = colnames(mat)[ind[, 2]],
    correlation = mat[ind]
  )
}))

# Pair names
corLong$pair <- paste(corLong$var1, corLong$var2, sep = " - ")

### Distribution by item ####
ggplot(corLong, aes(x = var1, y = correlation, group = var1)) +
  geom_boxplot(aes(fill = var1), alpha = 0.5) +
  theme_bw() +
  labs(
    x = "Pairs within each item",
    y = "Correlation"
  ) +
  scale_fill_manual(values = viridis::viridis(n = 6)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### Distribution by pairs ####
ggplot(corLong, aes(x = pair, y = correlation, color = var1)) +
  geom_boxplot(aes(fill = var1, colour = var1), alpha = 0.5) +
  theme_bw() +
  labs(
    x = "All pairs",
    y = "Correlation"
  ) +
  scale_fill_manual(values = viridis::viridis(n = 6)) +
  scale_colour_manual(values = viridis::viridis(n = 6)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )





#...........................................................................####
# LOO ####

## RANKING WITH WEITGHTED ITEM 2 ####

### Ranking with all items ####
rankingNS <- lapply(resNS, function(df){
  scoreDF_norm <- df |> dplyr::arrange(TotalCorrected)
  finalRank <- stats::setNames(scoreDF_norm$TotalCorrected, rownames(scoreDF_norm))
})
rankingNoItem0 <- lapply(resNS, function(df){
  scoreDF_norm <- df |> dplyr::arrange(Total)
  finalRank <- stats::setNames(scoreDF_norm$Total, rownames(scoreDF_norm))
})


### Ranking removing one item ####
rankingLOO <- lapply(1:6, function(i) {
  lapply(resNS, function(df){
  recalculateNS(df, itemToRemove = i)
  })
})
names(rankingLOO) <- paste0("Item ", 1:6)

rankingLOO$`Item 0` <- rankingNoItem0

### Comparison with kendall cor and graphs ####
res <- lapply(names(rankingLOO), function(i) 
  compareNSvsLOO(rankNS = rankingNS, rankLOO = rankingLOO[[i]], item = i))
names(res) <- names(rankingLOO)
ggpubr::ggarrange(plotlist = sapply(res, "[[", 3), nrow = 2, ncol = 4)
ggpubr::ggarrange(plotlist = sapply(res, "[[", 4), nrow = 2, ncol = 4, common.legend = TRUE)


# Testing...
dfKendall <- do.call(rbind, sapply(res, "[[", 2, USE.NAMES = T, simplify = F))

# Test first
dfTest <- dfKendall |>
  rstatix::pairwise_wilcox_test(
    kendall ~ item,
    paired = TRUE,
    p.adjust.method = "holm"
  )

# Plot without test

ggplot(dfKendall, aes(x = item, y = kendall)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Kendall's tau"
  )

# plot with test result

dfTest <- dfTest %>%
  rstatix::add_xy_position(x = "item")
dfTest$y.position <- 1+seq(0.1, 1.1, 0.05)

ggplot(dfKendall, aes(x = item, y = kendall)) +
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





## RANKING WITH NON-WEITGHTED ITEM 2 ####

### Ranking with all items ####
rankingNS <- lapply(resNsNoWeight, function(df){
  scoreDF_norm <- df |> dplyr::arrange(TotalCorrected)
  finalRank <- stats::setNames(scoreDF_norm$TotalCorrected, rownames(scoreDF_norm))
})
rankingNoItem0 <- lapply(resNsNoWeight, function(df){
  scoreDF_norm <- df |> dplyr::arrange(Total)
  finalRank <- stats::setNames(scoreDF_norm$Total, rownames(scoreDF_norm))
})


### Ranking removing one item ####
rankingLOO <- lapply(1:6, function(i) {
  lapply(resNsNoWeight, function(df){
  recalculateNS(df, itemToRemove = i)
  })
})
names(rankingLOO) <- paste0("Item ", 1:6)

rankingLOO$`Item 0` <- rankingNoItem0

### Comparison with kendall cor and graphs ####
res <- lapply(names(rankingLOO), function(i) 
  compareNSvsLOO(rankNS = rankingNS, rankLOO = rankingLOO[[i]], 
                 item = i, maxColor = "darkgreen"))
names(res) <- names(rankingLOO)
ggpubr::ggarrange(plotlist = sapply(res, "[[", 3), nrow = 2, ncol = 4)
ggpubr::ggarrange(plotlist = sapply(res, "[[", 4), nrow = 2, ncol = 4, common.legend = TRUE)


# Testing...
dfKendall <- do.call(rbind, sapply(res, "[[", 2, USE.NAMES = T, simplify = F))

# Test first
dfTest <- dfKendall |>
  rstatix::pairwise_wilcox_test(
    kendall ~ item,
    paired = TRUE,
    p.adjust.method = "holm"
  )

# Plot without test

ggplot(dfKendall, aes(x = item, y = kendall)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Kendall's tau"
  )

# plot with test result
library(ggpubr)

dfTest <- dfTest %>%
  rstatix::add_xy_position(x = "item")
dfTest$y.position <- 1+seq(0.1, 1.1, 0.05)

ggplot(dfKendall, aes(x = item, y = kendall)) +
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



### Additional result: discrimination ability ####
ordenNormalizations <- c(
  "Log", 
  "Median", 
  "Mean", 
  "TI",
  "Quantile", 
  "CyclicLoess",
  "RLR",
  "VSN" 
)
rankList <- lapply(seq_len(6), function(item) {
  
  dfRank <- do.call(rbind, lapply(resNsNoWeight, function(df) {
    df[ordenNormalizations, item]
  }))
  
  dfRank <- as.data.frame(dfRank)
  
  # Nombres de las normalizaciones
  colnames(dfRank) <- ordenNormalizations
  
  dfRank
})

names(rankList) <- paste0("Item", 1:6)
dim(rankList[[1]])
head(rankList[[1]])


heatmapList <- lapply(seq_along(rankList), function(i) {
  
  df <- rankList[[i]]
  
  dfLong <- data.frame(
    dataset = rep(seq_len(nrow(df)), times = ncol(df)),
    normalization = rep(colnames(df), each = nrow(df)),
    rank = unlist(df, use.names = FALSE)
  )
  
  ggplot(dfLong, aes(x = normalization, y = dataset, fill = rank)) +
    geom_tile() +
    scale_fill_gradient(
      # low = "#543005",
      low = "#F5F5F5",
      high = "#003C30",
      limits = c(0, 1),
      name = "Scaled item-score"
    ) +
    # scale_fill_gradientn(
    #   colours = RColorBrewer::brewer.pal(9, "PuOr"), 
    #   limits = c(0, 1), 
    #   name = "Scaled item-score",
    #   values = scales::rescale(c(0, 0.5, 1)))+
    labs(
      title = paste("Item", i),
      x = NULL,
      y = "Dataset"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), 
      panel.grid = element_blank()
    )
})

ggpubr::ggarrange(plotlist = heatmapList, nrow = 2, ncol = 3, common.legend = TRUE)
