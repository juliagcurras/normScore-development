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
resNS <- readRDS(file = "AssessmentFiles/3_normScore_results_100datasets.rds")

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
dfCorWeight <- apply(corArray, c(1, 2), median, na.rm = TRUE)

saveRDS(
  list(corList = corListWeight, dfCor = dfCorWeight), 
  "AssessmentFiles/4_nonRedundancyResults.rds")


pH1 <- pheatmap::pheatmap(
  dfCorWeight, 
  cluster_cols = TRUE,  
  cluster_rows = TRUE,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 14,
  border_color = "white",
  angle_col = 45,
  # clustering_callback = callback, 
  # color = grDevices::colorRampPalette(c(rev(RColorBrewer::brewer.pal(n = 9, name = "PuBu")), "white", "white", "white"))(319), # Top 1
  color = grDevices::colorRampPalette(rev(viridis::inferno(n =5)))(22), # Top 1
  # color = RColoseq(-1, 1, 0.5)rBrewer::brewer.pal(n = 6, name = "PuBu"), # Average -->
  # color = grDevices::colorRampPalette(colores)(351), # Average -->
  breaks = seq(0, 1, 0.05),
  legend_breaks = seq(0, 1, 0.25),
  legend_labels = seq(0, 1, 0.25)
)
svg(filename = paste0(outputDir, "NewFigure4A.svg"), width = 6, height = 6)
pH1
dev.off()


# With item 2 (without weights)
corListNoWeight <- lapply(resNsNoWeight, function(df){
  df <- df[,1:6]
  cor(df, method = "spearman")
})

corArrayNoW <- simplify2array(corListNoWeight)
dfCorNoWeight <- apply(corArrayNoW, c(1, 2), median, na.rm = TRUE)

pheatmap::pheatmap(dfCorNoWeight)

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
corLong |> 
  rename(kendall = correlation) |> 
  summarise(
    mean = mean(kendall), 
    sd = sd(kendall), 
    median = median(kendall), 
    IQR = IQR(kendall), 
    min = min(kendall), 
    max = max(kendall), 
    q1 = quantile(kendall)[[2]],
    q3 = quantile(kendall)[[4]])

# mean        sd    median       IQR       min max
# 0.4036418 0.4277955 0.4910268 0.5584513 -0.952381   1
# mean        sd    median       IQR       min max        q1        q3
# 1 0.4036418 0.4277955 0.4910268 0.5584513 -0.952381   1 0.1796439 0.7380952

### Distribution by item ####
ggplot(corLong, aes(x = var1, y = correlation, group = var1)) +
  geom_boxplot(aes(fill = var1), alpha = 0.75) +
  theme_bw() +
  labs(
    x = "Pairs within each item",
    y = "Correlation"
  ) +
  scale_fill_manual(values = viridis::inferno(n =6)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### Distribution by pairs ####
outputDir <- "G:/Mi unidad/_General/DoctoradoIndustrial2022/Publicaciones/NORMALIZACIONES/Paper_normScore/Figures/"

pByPairs <- ggplot(corLong, aes(x = pair, y = correlation, color = var1)) +
  # geom_boxplot(aes(fill = var1, colour = var1), alpha = 0.8) +
  geom_boxplot(aes(fill = var1), colour = "black", alpha = 0.8) +
  theme_bw(base_size = 18) +
  labs(
    x = "All pairs of items",
    y = "Spearman's correlation", 
    fill = "Item"
  ) +
  scale_fill_manual(values = viridis::inferno(n =6)) +
  # scale_colour_manual(values = viridis::inferno(n =6)) +
  theme(
    # panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pByPairs

ggsave(filename = paste0(outputDir, "NewFigure4B.svg"), 
       plot = pByPairs, width = 15, height = 6)

# jpeg(filename = paste0(outputDir, "Figure6.jpeg"), 
#      width = 6, height = 6, units = "in", res = 1000)
# pC
# dev.off()




