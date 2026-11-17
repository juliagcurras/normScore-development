
# Images for papers
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore-development/R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


# Global variables ####
outputDir <- "G:/Mi unidad/_General/DoctoradoIndustrial2022/Publicaciones/NORMALIZACIONES/Paper_normScore/FiguresNew/"
normOrder <- c("Log", "Mean", "TI", "Median", "Quantile", "CyclicLoess", "VSN", "RLR")

# Functions ####
callback = function(hc, mat){
    # d_rows <- dist(mat)
  d_rows <- as.dist(1 - cor(t(mat), use = "pairwise.complete.obs"))
  ord_rows <- seriation::get_order(seriation::seriate(d_rows, method = "OLO"))
  dend_rows <- as.dendrogram(hc)
  dend_rows <- reorder(dend_rows, wts = ord_rows)
  as.hclust(dend_rows)
}

dfToMat <- function(df, mode = "top1"){
  colnames(df)[4] <- "Ranking"
  if (mode == "top1"){
    matRanking <- df %>%
      filter(Ranking == 1) %>%
      group_by(Normalization, Item) %>%
      summarise(sumRanking = sum(Ranking, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = Item,
        values_from = sumRanking,
        values_fill = 0
      ) %>%
      tibble::column_to_rownames("Normalization") %>%
      as.matrix() %>%
      t()
  } else if (mode == "average"){
    matRanking <- df %>%
      mutate(rankScore = 9 - Ranking) %>%  # 8 = mejor, 1 = peor
      mutate(rankScore = Ranking) %>% 
      group_by(Normalization, Item) %>%
      summarise(meanScore = mean(rankScore, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = Item,
        values_from = meanScore
      ) %>%
      tibble::column_to_rownames("Normalization") %>%
      as.matrix() %>%
      t()
  }
  
  # matRanking <- matRanking[normOrder, paste0("Item", 1:6)]
  matRanking <- matRanking[paste0("Item", 1:6), normOrder]
  return(matRanking)
}



#.........................................................................####
# Figure 2B ####

## Required functions ####
source(file = "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/3_AssessmentFunctions.R")
inputDir <- "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/AssessmentFiles/"

## Load data ####
resAllGS <- readRDS(file = paste0(inputDir, "3_Reviewer1_results_100dataset.rds")) # Manual ranking
resAllGSR2 <- readRDS(file = paste0(inputDir, "3_Reviewer2_results_100dataset.rds")) # Manual ranking
resNS <- readRDS(file = paste0(inputDir, "3_normScore_results_100datasets.rds")) # Originales
resAllNS <- formatNormScoreResults(finalList = resNS, item0 = T)
resAllNS$resByItem <- resAllNS$resByItem %>%
  group_by(ID, Item) %>%
  mutate(Ranking = min_rank(Ranking)) %>%
  ungroup() #%>%
colores <- c("#A50026","#F4A582", "lightyellow", "#92C5DE", "#313695","#171945")
# colores <- RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")

## Figures By Item ####

# Width: 400 | Height: 309 => saving options

### R1 ####
mat <- dfToMat(resAllGS$byItem, mode = "top1")
mat <- dfToMat(resAllGS$byItem, mode = "average")
pH1 <- pheatmap::pheatmap(
  mat = mat, 
  cluster_cols = F,  
  cluster_rows = F,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 14,
  border_color = "white",
  angle_col = 45,
  clustering_callback = callback, 
  color = grDevices::colorRampPalette(colores)(351), # Average -->
  breaks = seq(1, 8, 0.02),
  legend_breaks = 1:8,
  legend_labels = as.character(1:8)
)
pH1

### R2 ####
mat <- dfToMat(resAllGSR2$byItem, mode = "top1")
mat <- dfToMat(resAllGSR2$byItem, mode = "average")
pH2 <- pheatmap::pheatmap(
  mat = mat, 
  cluster_cols = F,  
  cluster_rows = F,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 14,
  angle_col = 45,
  border_color = "white",
  clustering_callback = callback, 
  # color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PuBu"))(200)# Top 1
  color = grDevices::colorRampPalette(colores)(351), # Average -->
  breaks = seq(1, 8, 0.02),
  legend_breaks = 1:8,
  legend_labels = as.character(1:8)
)



### NS ####
mat <- dfToMat(resAllNS$resByItem, mode = "top1")
mat <- dfToMat(resAllNS$resByItem, mode = "average")
min(mat)
max(mat)
pHns <- pheatmap::pheatmap(
  mat = mat, 
  cluster_cols = F,  
  cluster_rows = F,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 14,
  border_color = "white",
  angle_col = 45,
  clustering_callback = callback, 
  # color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PuBu"))(200) # Top 1
  color = grDevices::colorRampPalette(colores)(351), # Average -->
  breaks = seq(1, 8, 0.02),
  legend_breaks = 1:8,
  legend_labels = as.character(1:8)
)
pHns





## Figures Global ####

### NS ####

#### Ranking ####
metodos <- unique(unlist(resAllNS$resGlobal))

matRanking <- t(sapply(resAllNS$resGlobal, function(x) {
  match(metodos, x)
}))

colnames(matRanking) <- metodos
rownames(matRanking) <- names(resAllNS$resGlobal)

matRanking <- 9-matRanking

pG <- pheatmap::pheatmap(
  mat = matRanking, 
  cluster_cols = T,  
  cluster_rows = T,
  show_rownames  = F,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 14,
  clustering_callback = callback,
  border_color = NA, 
  angle_col = 45,
  # treeheight_row = 0,
  # treeheight_col = 0,
  # color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "PuBu"))(200) # Top 1
  color = RColorBrewer::brewer.pal(n = 9, name = "PuBu"), # Average -->
  breaks = 0:8,
  legend_breaks = 1:8,
  legend_labels = as.character(8:1)
)
pG

# svg(filename = paste0(outputDir, "Figure3/Hetamap_NS_Global.svg"), width = 6, height = 5)
# pG
# dev.off()


#### Puntuaciones ####
listaScores <- resNS
dfScores <- bind_rows(
  lapply(listaScores, function(x) {
    x %>%
      tibble::rownames_to_column("Metodo") %>%
      rename(Score = last_col()) %>%
      select(Metodo, Score)
  }),
  .id = "ID"
) %>%
  pivot_wider(
    names_from = Metodo,
    values_from = Score
  ) %>%
  tibble::column_to_rownames("ID") %>%
  as.matrix
max(dfScores)
min(dfScores)
summary(dfScores)

colores <- c("#6C001A", "#A50026","#F4A582", "lightyellow", "#92C5DE", "#313695", "#171945")
# colores <- c("aliceblue", "#92C5DE", "#313695","#171945")
colorsBase <- rev(colorRampPalette(
  rev(colores)
  # RColorBrewer::brewer.pal(8, "PuBu")
)(100))

colors <- c(colorsBase, tail(colorsBase, 1))

breaks <- c(
  seq(0, 5, length.out = length(colorsBase) + 1)
)


pG <- pheatmap::pheatmap(
  mat = dfScores, 
  cluster_cols = T,  
  cluster_rows = T,
  show_rownames  = F,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 14,
  clustering_callback = callback,
  angle_col = 45,
  color = colors,
  breaks = breaks,
  treeheight_row = 20,
  treeheight_col = 20,
  legend_breaks = 0:5,
  legend_labels = c("0", "1", "2", "3", "4", ">5"),
  border_color = NA
)
pG


### Reviewer 1 ####
metodos <- normOrder
dfConteo <- bind_rows(
  imap_dfr(resAllGS$globalBests, ~ tibble(ID = .y, Metodo = .x, Tipo = "Best")),
  imap_dfr(resAllGS$globalWorst, ~ tibble(ID = .y, Metodo = .x, Tipo = "Worst"))
) %>%
  distinct(ID, Metodo, Tipo) %>%
  count(Tipo, Metodo, name = "nID") %>%
  complete(
    Tipo = c("Best", "Worst"),
    Metodo = metodos,
    fill = list(nID = 0)
  )

dfConteo$Metodo <- factor(dfConteo$Metodo, levels =  normOrder) 

p1 <- ggplot(dfConteo, aes(x = Metodo, y = nID, fill = Tipo)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.95) +
  geom_text(
    aes(label = nID),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 5
  ) +
  scale_fill_manual(values = c("#A50026", "#313695")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    limits = c(0, 100),
    breaks = seq(0, 100, 25)
  ) +
  labs(x = "", y = "No. of datasets", fill = "Selected as...") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 18), legend.position = "top")
p1




### Reviewer 2 ####
dfConteo <- bind_rows(
  imap_dfr(resAllGSR2$globalBests, ~ tibble(ID = .y, Metodo = .x, Tipo = "Best")),
  imap_dfr(resAllGSR2$globalWorst, ~ tibble(ID = .y, Metodo = .x, Tipo = "Worst"))
) %>%
  distinct(ID, Metodo, Tipo) %>%
  count(Tipo, Metodo, name = "nID") %>%
  complete(
    Tipo = c("Best", "Worst"),
    Metodo = metodos,
    fill = list(nID = 0)
  )

dfConteo$Metodo <- factor(dfConteo$Metodo, levels =  normOrder) 

p2 <- ggplot(dfConteo, aes(x = Metodo, y = nID, fill = Tipo)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.95) +
  geom_text(
    aes(label = nID),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 5
  ) +
  scale_fill_manual(values = c("#A50026", "#313695")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    limits = c(0, 100),
    breaks = seq(0, 100, 25)
  ) +
  labs(x = "", y = "No. of datasets", fill = "Selected as...") +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 18), legend.position = "top")
p2





# Only best #
listaMetodos <- resAllGSR2$globalBests
dfConteo <- data.frame(
  Metodo = c(names(table(unlist(listaMetodos))), "TI", "RLR"),
  nID = c(as.integer(table(unlist(listaMetodos))), 0, 0)
)
rownames(dfConteo) <- dfConteo$Metodo
dfConteo <- dfConteo[normOrder, ]
dfConteo$Metodo <- factor(dfConteo$Metodo, levels =  normOrder) 

p2 <- ggplot(dfConteo, aes(x = Metodo, y = nID)) +
  geom_col(colour = "#313695", fill = "#313695") + 
  geom_text(aes(label = nID), vjust = -0.3, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "", y = "No. of datasets") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2
# svg(filename = paste0(outputDir, "Figure3/Barplot_R2_Global.svg"), width = 5.5, height = 5.7)
# p2
# dev.off()







# Material supp
listaScores <- resNS
itemsKeep <- c(
  "Item1", "Item2", "Item3", "Item4", "Item5", "Item6",
  "TotalCorrected"
)

listaPorItem <- map(itemsKeep, function(item) {
  
  map_dfr(listaScores, function(df) {
    df %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Normalizacion") %>%
      select(Normalizacion, Score = all_of(item)) %>%
      tidyr::pivot_wider(
        names_from = Normalizacion,
        values_from = Score
      ) %>%
      select(all_of(ordenNorm))
    
  }, .id = "ID")
  
})

names(listaPorItem) <- itemsKeep

writexl::write_xlsx(x = listaPorItem, path = paste0(outputDir, "SupplTableNormScore.xlsx"))






#.........................................................................####
# Figure 1. normScore Scheme 3 ####

## Requiered functions ####
source(file = "simulationFunction.R")
# lista_plots: list de objetos ggplot
getYrange <- function(p) {
  b <- ggplot2::ggplot_build(p)
  ys <- unlist(lapply(b$data, function(d) d$y), use.names = FALSE)
  range(ys[is.finite(ys)], na.rm = TRUE)
}


getResultsByItem <- function(dfRaw, listaNorm, dm, item, colores = "#003C72", 
                             numCols = 3, numRows = 1){
  # just a common object
  # dm <- as.data.frame(lista[[1]][["metadata"]])
  # cat(item, "\n")
  colnames(dm) <- c("Samples", "Groups")
  grupos <- unique(dm$Groups)
  # numRows <- ceiling(length(listaNorm)/2)
  
  # Starting with data extraction and metrics/graphs estimation... 
  if (item == "item0"){ # only item 0
    output <- Biomics::plotBarTI(data = dfRaw, interact = F, color = colores)$grafico +
    # output <- plotBarTI1(data = dfRaw, interact = F, color = colores)$grafico + 
                   theme(axis.text.x = element_text(vjust = 0.5))
    # output <-  ggpubr::ggarrange(plotlist = empties, ncol = 1, nrow = numRows-2)
  } else { # other items
    # First extrating individual data and generating individual plots
    a <<- 0
    finalData <- lapply(listaNorm, function(datos){
      a <<- a+1
      cores <- rep(colores[a], nrow(dm))
      # cat(a, "\n")
      # cat(names(listaNorm)[a], "\n")
      nome <- names(listaNorm)[a]
      if (item == "item6"){
        Biomics::plotBoxMulti(
          base = datos, varResumen = colnames(datos), color = cores, 
          interact = F, titulo = nome, tituloX = "Samples",
          tituloY = ifelse(a == 1, "Intensity distribution", "")
        )$grafico + 
          theme_classic() +
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      } else if (item == "item5"){
        Biomics::plotRLE(
          df = datos, normalizacion = nome, interact = F, 
          titulo = nome, tituloX = "Samples", color = cores, 
          tituloY = ifelse(a == 1, "RLE", ""))$grafico  + 
          theme_classic() +
          ylim(c(-3, 3)) +
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      } else if (item == "item4"){
        Biomics::plotMeanSD(df = datos, interact = F, titulo = nome,
                            color = cores)$grafico  + 
          theme_classic() +
          ylim(c(1, 1.5)) +
          ylab(ifelse(a == 1, "SD", ""))+
          xlab("Sorted mean")+
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      } else if (item == "item3"){
        Biomics::plotMA(df = datos, dfGrupos = dm, gControl = grupos[1], 
                        gCase = grupos[2], titulo = nome, limY = c(-2, 2),
                        showR2 = F, interact = F, color = cores)$grafico +
          theme_classic()+
          ylab(ifelse(a == 1, "Mean ratio", ""))+
          xlab("Average mean")+
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      }
    })
    names(finalData) <- paste0("SimDataset_",1:length(finalData))
    
    # Select ylim for MAplot after calculations
    if (item == "item3"){
      rangos <- lapply(finalData, getYrange)
      ylimGlobal <- ceiling(max(abs(range(unlist(rangos), na.rm = TRUE)*1.1)))
      # Aplicar el mismo límite Y a todos (sin recortar datos del todo)
      finalData <- lapply(finalData, function(p) {
        p  + coord_cartesian(ylim = c(-ylimGlobal, ylimGlobal))
      })
    }
    
    # Mix plots into a single image or...
    if (item %in% c("item3", "item4", "item5", "item6")){ # plot graphs together: MAplot, RLEplot, meanSDplot, TIboxplot
      output <- ggpubr::ggarrange(plotlist = finalData, 
                                  ncol = numCols, nrow = numRows)
    } else if (item == "item1"){ # ...or generate other metrics with data for remaining items
      output <- Biomics::getPCV( # PVC
        listData = listaNorm,
        grupos = grupos,
        dfGrupos = dm,
        grafico = T,
        interact = F
      )$grafico
      # output <-  ggpubr::ggarrange(plotlist = empties, ncol = 1, nrow = numRows-2)
    } else if (item =="item2"){ # Correlation
      allVectorsCorr <- lapply(
        listaNorm, 
        Biomics::getPooledCor, 
        dfGrupos = dm,
        method = "spearman")
      dfPlot <- data.frame(sapply(allVectorsCorr, "length<-", max(lengths(allVectorsCorr))))
      output <- Biomics::plotBoxMulti(
        base = dfPlot, 
        varResumen = colnames(dfPlot),
        tituloX = "Normalizations", interact = F,
        tituloY = "Spearman correlation")$grafico
      # output <- ggpubr::ggarrange(plotlist = empties, ncol = 1, nrow = numRows-2)
    }
    
  } 
  res <- list(output)
  names(res) <- item
  return(res)
}


## Data for items 0, 1, 5, 6 ####

# Parameters
item0156_sample_shift_sd <- seq(0.5, 1.5, 0.5)
item0156_sample_shift_sd <- seq(0, 1, 0.5)
item0156_sample_shift_cap <- item0156_sample_shift_sd + 0.05
data <- list()

# Simulate data for example
for(i in 1:length(item0156_sample_shift_sd)){
  data[[i]] <- simulate_proteomics_clean(
    semilla = 9693, 
    n_proteins = 1000,
    n_per_group = 5, # 15 for Item 1
    sample_shift_sd = item0156_sample_shift_sd[i],
    sample_shift_cap = item0156_sample_shift_cap[i]
  )
}

listaNorm <- sapply(data, "[[", 1, simplify = F)
names(listaNorm) <- c("Norm2", "Norm1", "Norm3")
listaNorm <- listaNorm[c("Norm1", "Norm2", "Norm3")]
dfRaw <- data[[3]][["rawData"]]
dm <- data[[3]][["metadata"]]
colores <- c("#003C72", "#1786A3", "#2FB2AD")

### Item 0 ####
p0 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listNorm, 
  dm = dm, 
  item = "item0")$item0 + 
  theme_classic() +
  ylab("Total intensity\n (raw data)")+
  xlab("Samples")+
  scale_fill_manual(values = "#8076B9")+
  ggplot2::theme(
    legend.position = "none", 
    axis.ticks = element_blank(),
    axis.title = element_text(size = 26), 
    axis.text = element_blank(), 
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))
ggsave(filename = paste0(outputDir, "normScoreScheme/item0.svg"), plot = p0, 
       width = 7.5, height = 3.2, bg = "transparent")

### Item 1 ####
dmAlt <- dm # Hey! Reexecute simulated data with n=15 per group
dmAlt$Groups <- rep(paste0("G", 1:5), each = 6)
pI1 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listaNorm, 
  dm = dmAlt, 
  item = "item1")$item1 + 
  theme_classic()+
  ylab("PCV")+
  ggplot2::theme(
    legend.position = "none", 
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 22), 
    axis.text.y = element_blank(), 
    axis.text.x = element_text(hjust = 0.5,size = 18))
pI1
ggsave(filename = paste0(outputDir, "normScoreScheme/item1.svg"), plot = pI1, 
       width = 5, height = 3)


### Item 5 ####
pI5 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listaNorm, 
  dm = dm, 
  item = "item5", 
  colores = colores)$item5
pI5
ggsave(filename = paste0(outputDir, "normScoreScheme/item5.svg"), plot = pI5, 
       width =8, height = 3.5)

### Item 6 ####
pI6 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listaNorm, 
  dm = dm, 
  item = "item6", 
  colores = colores)$item6
pI6
ggsave(filename = paste0(outputDir, "normScoreScheme/item6.svg"), plot = pI6, 
       width =8, height = 3.5)


## Data for item 2 ####
# Parameters
item2_values_rho_between <- c(0.15, 0.3, 0.6)
item2_values_rho_within <- c(0.45, 0.45, 1)
data <- list()

# Simulate data for example
for(i in 1:length(item0156_sample_shift_sd)){
  data[[i]] <- simulate_proteomics_clean(
    semilla = 9693, 
    n_proteins = 1000,
    n_per_group = 15, # 15 for Item 1
    rho_between = item2_values_rho_between[i],
    rho_within = item2_values_rho_within[i]
  )
}
listaNorm <- sapply(data, "[[", 1, simplify = F)
names(listaNorm) <- c("Norm1", "Norm3", "Norm2")
listaNorm <- listaNorm[c("Norm1", "Norm2", "Norm3")]
dfRaw <- data[[3]][["rawData"]]
dm <- data[[3]][["metadata"]]
colores <- c("#003C72", "#1786A3", "#2FB2AD")

### Item 2 ####
dmAlt <- dm # Hey! Reexecute simulated data with n=15 per group
dmAlt$Groups <- rep(paste0("G", 1:5), each = 6)
pI2 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listaNorm, 
  dm = dmAlt, 
  item = "item2")$item2 + 
  theme_classic()+
  ylim(c(0.85, 0.99))+
  ggplot2::theme(
    legend.position = "none", 
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.text.x = element_text(hjust = 0.5,size = 18))
pI2
ggsave(filename = paste0(outputDir, "normScoreScheme/item2.svg"), plot = pI2, 
       width = 5.8, height = 3.1)
#






## Data for item 4 ####
# Parameters
item4_sample_sd_cap <- c(0, 0.05, 0.03)
item4_sample_sd_strength <- c(1, 1, -0.5)
item4_sample_sd_rho <- c(1.5, 1.5, 1.2)


# Simulate data for example
for(i in 1:length(item4_sample_sd_cap)){
  data[[i]] <- simulate_proteomics_clean(
    semilla = 9693, 
    n_proteins = 1000,
    n_per_group = 15, # 15 for Item 1
    sample_shift_sd = 0.5,
    sample_sd_cap = item4_sample_sd_cap[[i]],
    sample_sd_strength = item4_sample_sd_strength[[i]],
    sample_sd_rho = item4_sample_sd_rho[[i]]
    )
}
listaNorm <- sapply(data, "[[", 1, simplify = F)
names(listaNorm) <- c("Norm1", "Norm3", "Norm2")
listaNorm <- listaNorm[c("Norm1", "Norm2", "Norm3")]
dfRaw <- data[[3]][["rawData"]]
dm <- data[[3]][["metadata"]]
colores <- c("#003C72", "#1786A3", "#2FB2AD")

### Item 4 ####
pI4 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listaNorm, 
  dm = dm, 
  item = "item4", 
  colores = colores)$item4
pI4

ggsave(filename = paste0(outputDir, "normScoreScheme/item4.svg"), plot = pI4, 
       width =8, height = 3.5)





## Data for item 3 ####

item3_sigma_lo <- c(1, 1, 1) #*1.5
item3_sigma_hi <- c(0.1, 0.1, 0.1)
item3_sample_sd_cap <- c(0.35, 0.5, 0.6)
item3_sample_sd_strength = c(0, 2.5, 3)
item3_sample_sd_rho = c(0.8, 0.2, 0)

# Simulate data for example
for(i in 1:length(item3_sample_sd_cap)){
  data[[i]] <- simulate_proteomics_clean(
    semilla = 9693, 
    n_proteins = 1500,
    n_per_group = 15, # 15 for Item 1
    sample_shift_sd = 0.5,
    prop_de = 0.2,
    sigma_lo = item3_sigma_lo[[i]],
    sigma_hi = item3_sigma_hi[[i]], 
    sample_sd_strength = item3_sample_sd_strength[[i]],
    sample_sd_rho = item3_sample_sd_rho[[i]],
    sample_sd_cap = item3_sample_sd_cap[[i]]
  )
}
listaNorm <- sapply(data, "[[", 1, simplify = F)
names(listaNorm) <- c("Norm1", "Norm2", "Norm3")
listaNorm <- listaNorm[c("Norm1", "Norm2", "Norm3")]
dfRaw <- data[[3]][["rawData"]]
dm <- data[[3]][["metadata"]]
colores <- c("#003C72", "#1786A3", "#2FB2AD")

## Item 3 ####
pI3 <- getResultsByItem(
  dfRaw = dfRaw, 
  listaNorm = listaNorm, 
  dm = dm, 
  item = "item3", 
  colores = colores)$item3
pI3

ggsave(filename = paste0(outputDir, "normScoreScheme/item3.svg"), plot = pI3, 
       width =8, height = 3.5)





#.........................................................................####
# Supplementary Figure 1. Simulation Scheme ####

## Requiered functions ####
source(file = "simulationFunction.R")


getResultsByItem <- function(dfRaw, listaNorm, dm, item, colores = "#003C72", 
                             numCols = 3, numRows = 1){
  # just a common object
  # dm <- as.data.frame(lista[[1]][["metadata"]])
  # cat(item, "\n")
  colnames(dm) <- c("Samples", "Groups")
  grupos <- unique(dm$Groups)
  # numRows <- ceiling(length(listaNorm)/2)
  
  # Starting with data extraction and metrics/graphs estimation... 
  if (item == "item0"){ # only item 0
    output <- Biomics::plotBarTI(data = dfRaw, interact = F, color = colores)$grafico +
      # output <- plotBarTI1(data = dfRaw, interact = F, color = colores)$grafico + 
      theme(axis.text.x = element_text(vjust = 0.5))
    # output <-  ggpubr::ggarrange(plotlist = empties, ncol = 1, nrow = numRows-2)
  } else { # other items
    # First extrating individual data and generating individual plots
    a <<- 0
    finalData <- lapply(listaNorm, function(datos){
      a <<- a+1
      cores <- rep(colores, nrow(dm))
      # cat(a, "\n")
      # cat(names(listaNorm)[a], "\n")
      nome <- names(listaNorm)[a]
      if (item == "item6"){
        Biomics::plotBoxMulti(
          base = datos, varResumen = colnames(datos), color = cores, 
          interact = F, titulo = nome, tituloX = "Samples",
          tituloY = ifelse(a == 1, "Intensity distribution", "")
        )$grafico + 
          theme_classic() +
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      } else if (item == "item5"){
        Biomics::plotRLE(
          df = datos, normalizacion = nome, interact = F, 
          titulo = nome, tituloX = "Samples", color = cores, 
          tituloY = ifelse(a == 1, "RLE", ""))$grafico  + 
          theme_classic() +
          ylim(c(-3, 3)) +
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      } else if (item == "item4"){
        Biomics::plotMeanSD(df = datos, interact = F, titulo = nome,
                            color = cores)$grafico  + 
          theme_classic() +
          ylim(c(1, 1.5)) +
          ylab(ifelse(a == 1, "SD", ""))+
          xlab("Sorted mean")+
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      } else if (item == "item3"){
        Biomics::plotMA(df = datos, dfGrupos = dm, gControl = grupos[1], 
                        gCase = grupos[2], titulo = nome, limY = c(-2, 2),
                        showR2 = F, interact = F, color = cores)$grafico +
          theme_classic()+
          ylab(ifelse(a == 1, "Mean ratio", ""))+
          xlab("Average mean")+
          ggplot2::theme(
            legend.position = "none",
            title = element_text(size = 20),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 18))
      }
    })
    names(finalData) <- paste0("SimDataset_",1:length(finalData))
    
    # Select ylim for MAplot after calculations
    if (item == "item3"){
      rangos <- lapply(finalData, getYrange)
      ylimGlobal <- ceiling(max(abs(range(unlist(rangos), na.rm = TRUE)*1.1)))
      # Aplicar el mismo límite Y a todos (sin recortar datos del todo)
      finalData <- lapply(finalData, function(p) {
        p  + coord_cartesian(ylim = c(-ylimGlobal, ylimGlobal))
      })
    }
    
    # Mix plots into a single image or...
    if (item %in% c("item3", "item4", "item5", "item6")){ # plot graphs together: MAplot, RLEplot, meanSDplot, TIboxplot
      output <- ggpubr::ggarrange(plotlist = finalData, 
                                  ncol = numCols, nrow = numRows)
    } else if (item == "item1"){ # ...or generate other metrics with data for remaining items
      output <- Biomics::getPCV( # PVC
        listData = listaNorm,
        grupos = grupos,
        dfGrupos = dm,
        grafico = T,
        interact = F
      )$grafico
      # output <-  ggpubr::ggarrange(plotlist = empties, ncol = 1, nrow = numRows-2)
    } else if (item =="item2"){ # Correlation
      allVectorsCorr <- lapply(
        listaNorm, 
        Biomics::getPooledCor, 
        dfGrupos = dm,
        method = "spearman")
      dfPlot <- data.frame(sapply(allVectorsCorr, "length<-", max(lengths(allVectorsCorr))))
      output <- Biomics::plotBoxMulti(
        base = dfPlot, 
        varResumen = colnames(dfPlot),
        tituloX = "Normalizations", interact = F,
        tituloY = "Spearman correlation")$grafico
      # output <- ggpubr::ggarrange(plotlist = empties, ncol = 1, nrow = numRows-2)
    }
    
  } 
  res <- list(output)
  names(res) <- item
  return(res)
}



goldStandard <- paste0("Dataset ", 1:9)

### Item 0, 1, 5, 6 ####
item0156_sample_shift_sd <- list(
  rep(0, 9), 
  seq(0, 0.5, 0.0625), 
  # seq(0.5, 1.5, 0.125), 
  seq(1, 2.5, 0.1875),
  seq(0, 2.5, 0.3125))
names(item0156_sample_shift_sd) <- paste0("Sim", 1:length(item0156_sample_shift_sd))

offset <- c(0.2, 0.05, 0.05, -0.05)
item0156_sample_shift_cap <- sapply(1:length(item0156_sample_shift_sd), function(x){
  item0156_sample_shift_sd[[x]] + offset[x]
}, simplify = F, USE.NAMES = T)
names(item0156_sample_shift_cap) <- paste0("Sim", 1:length(item0156_sample_shift_cap))

simulacionN <- names(item0156_sample_shift_cap)
coloresSim <- c("darkgrey", "#1786A3", "#2FB2AD", "#003C72")
names(coloresSim) <- simulacionN
res <- sapply(simulacionN, function(simu) {
  # simu <- 1
  valores <- length(item0156_sample_shift_sd[[simu]])
  resByItem <- lapply(1:valores, function(val) simulate_proteomics_clean(
    semilla = 9396, 
    n_proteins = 2000,
    n_per_group = 5,
    sample_shift_sd = item0156_sample_shift_sd[[simu]][[val]],
    sample_shift_cap = item0156_sample_shift_cap[[simu]][[val]])
  )
  names(resByItem) <- goldStandard
  
  #--- Cheking item 5 ---
  listaNorm <- sapply(resByItem, "[[", 1, simplify = F)
  dm <- resByItem[[1]][["metadata"]]
  dfRaw <- resByItem[[1]][["rawData"]]
  getResultsByItem(
    dfRaw = dfRaw, 
    listaNorm = listaNorm, 
    dm = dm, 
    item = "item6", 
    numCols = 9, 
    numRows = 1, 
    # colores = rep(coloresSim[simu], 9)
    colores = coloresSim[[simu]]
    )
  }
)

ggsave(filename = paste0(outputDir, "simulationScheme/error1.svg"), 
       plot = res$Sim1.item6, 
       width = 25, height = 3)
ggsave(filename = paste0(outputDir, "simulationScheme/error2.svg"), 
       plot = res$Sim2.item6, 
       width = 25, height = 3)
ggsave(filename = paste0(outputDir, "simulationScheme/error3.svg"), 
       plot = res$Sim3.item6, 
       width = 25, height = 3)
ggsave(filename = paste0(outputDir, "simulationScheme/error4.svg"), 
       plot = res$Sim4.item6, 
       width = 25, height = 3)







#.........................................................................####
# Figure 3 A to F ####

plotBoxChulo <- function(
    data, 
    varResumen = "tau", 
    varGrupo = "Item", 
    tituloY = "Kendall's Tau",
    color = "#003C72", 
    colorMean = "#C02A4E",
    titulo = "Reviewer 1 vs normScore",
    tituloX = NULL, 
    sizeLetra = 10,
    semilla = 969,
    anguloX = 45, 
    sizeDots = 1,
    widthDots = 0.05,
    faceX = "italic"
    ){
  
  df <- data[, c(varResumen, varGrupo)]
  colnames(df) <- c("Metric", "Item")
  set.seed(semilla)
  ggplot(df, aes(x = Item, y = Metric)) +
    
    # Boxplot
    geom_boxplot(
      width = 0.6,
      fill = color,
      # color = color,
      alpha = 0.4,
      box.linewidth = 0.15,
      staple.linewidth = 0.15,
      outlier.shape = NA
    ) +
    
    # Puntos individuales
    ggrastr::rasterise(
      geom_jitter(
        width = widthDots,
        size = sizeDots,
        alpha = 0.6,
        color = color
      ), 
      dpi = 400) +
    
    # Línea de la mediana
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3,
      fill = colorMean,
      color = colorMean
    ) +
    
    # Etiquetas
    labs(
      x = tituloX,
      y = tituloY,
      title = titulo
    ) +
    ylim(c(-1, 1)) +
    # Tema
    theme_minimal(base_size = sizeLetra) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = sizeLetra, face = , 
                                 colour = "black", angle = anguloX, 
                                 hjust = ifelse(anguloX != 0, 1, 0)),
      # axis.ticks.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black")
    )
}

## Loading results ####
inputDir <- "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/AssessmentFiles/"
resAllGS <- readRDS(file = paste0(inputDir, "3_Reviewer1_results_100dataset.rds")) # Manual ranking
resNS <- readRDS(file = paste0(inputDir, "3_normScore_results_100datasets.rds")) # Originales
resAllGSR2 <- readRDS(file = paste0(inputDir, "3_Reviewer2_results_100dataset.rds")) # Manual ranking

source(file = "../OneHundredDatasets/3_AssessmentFunctions.R")
resAllNS <- formatNormScoreResults(finalList = resNS, item0 = T)
resAllNSR2 <- formatNormScoreResults(finalList = resNSR2, item0 = T)


listPlots <- list()
## Reviewer 1 vs normScore ####
summaryByItem <- resultsByItem(resNormScore = resAllNS$resByItem, 
                               resGS = resAllGS$byItem)
data <- summaryByItem$results
colnames(data)

col1 <- "#31688EFF" 
col2 <- "#DAA600"

# Plot #
listPlots[["p1"]] <- plotBoxChulo(
  data, varResumen = "spearman", varGrupo = "Item",
  tituloY = "Spearman's coefficient", color = col1,# "#005B9A",
  titulo = "Reviewer 1 vs normScore")
listPlots[["p4"]] <- plotBoxChulo(
  data, varResumen = "pairAccuracy", varGrupo = "Item",
  tituloY = "Pairwise Accuracy", color = col2, # "#46AC8A",
  titulo = "Reviewer 1 vs normScore")



## Reviewer 2 vs normScore ####
summaryByItemR2 <- resultsByItem(resNormScore = resAllNSR2$resByItem, 
                                 resGS = resAllGSR2$byItem)
data <- summaryByItemR2$results
colnames(data)

listPlots[["p2"]] <- plotBoxChulo(
  data, varResumen = "spearman", varGrupo = "Item",
  tituloY = "Spearman's coefficient", color = col1, # "#005B9A",
  titulo = "Reviewer 2 vs normScore")
listPlots[["p5"]] <- plotBoxChulo(
  data, varResumen = "pairAccuracy", varGrupo = "Item",
  tituloY = "Pairwise Accuracy", color = col2, #"#46AC8A",
  titulo = "Reviewer 2 vs normScore")


## Reviewer 1 vs Reviewer 2 ####
resAllGS$byItem <- resAllGS$byItem %>% rename(Ranking = RankGS)
summaryByItemR1R2 <- resultsByItem(resNormScore = resAllGS$byItem, 
                                   resGS = resAllGSR2$byItem)
data <- summaryByItemR1R2$results
colnames(data)

listPlots[["p3"]] <- plotBoxChulo(
  data, varResumen = "spearman", varGrupo = "Item",
  tituloY = "Spearman's coefficient", color = col1, # "#005B9A",
  titulo = "Reviewer 1 vs Reviewer 2")
listPlots[["p6"]] <- plotBoxChulo(
  data, varResumen = "pairAccuracy", varGrupo = "Item",
  tituloY = "Pairwise Accuracy", color =  col2,# "#46AC8A",
  titulo = "Reviewer 1 vs Reviewer 2")


p1 <- ggpubr::ggarrange(plotlist = listPlots, ncol = 2, nrow = 3, labels = "AUTO")
p1

ggsave(filename = paste0(outputDir, "Figure3_AtoF_.svg"), plot = p1, width = 5.3, height = 8)
# jpeg(filename = paste0(outputDir, "Figure4.jpeg"), 
#      width = 8, height = 12, units = "in", res = 600)
# p1
# dev.off()



listPlots <- listPlots[paste0("p", 1:6)]
p1 <- ggpubr::ggarrange(plotlist = listPlots, ncol = 3, nrow = 2, labels = "AUTO")
p1
ggsave(filename = paste0(outputDir, "Figure4Longer.svg"), plot = p1, width = 11, height = 7)




#.........................................................................####
# Figure 2A ####

df <- readRDS(file = "C:/Users/julia/Documents/GitHub/normScore-development/Simulations/allSimulations.rds")
df <- df %>% filter(!(especificos %in% c(1,4,8))) # quitamos o 1 que ten control negativo en item 3, o 4 pq é moi reiterativo co 3 e non aporta nada, e o 8 porque pa valores extremos en item3 a situación ponse rara.
df$especificos <- factor(df$especificos, levels= c(0, 2:3, 5:7), labels = LETTERS[1:6])


## Plot simulation results ####
colores <- Biostatech::colorPalette(n = 100)
colores <- c("#003C72","#00467F", "#005B9A","#0A6E9D", "#117CA1", "#249EA8")
names(colores) <- unique(df$especificos)
graficos <- sapply(unique(df$especificos), function(i){
  dfAux <- df %>% 
    filter(especificos == i) %>% 
    dplyr::select(Item0:normScore) %>% 
    tidyr::pivot_longer(cols = colnames(.), names_to = "Metrics", values_to = "Cor")
  plotBoxChulo(
    data = dfAux, 
    varResumen = "Cor", 
    varGrupo = "Metrics", 
    tituloY = "Kendall's Correlation",
    # color = colores[i], 
    anguloX = 65,
    sizeDots = 0.6,
    colorMean = "#C02A4E",
    titulo = paste0("Error magnitude ", i),
    sizeLetra = 10)

})


p1 <- ggpubr::ggarrange(plotlist = graficos, ncol = 3, nrow = 2)
ggsave(filename = paste0(outputDir, "Figure2A.svg"), plot = p1, width = 8, height = 5.5)

p1 <- ggpubr::ggarrange(plotlist = graficos, ncol = 2, nrow = 3)
p1
ggsave(filename = paste0(outputDir, "Figure2.svg"), plot = p1, width = 8, height = 12)

jpeg(filename = paste0(outputDir, "Figure2.jpeg"), 
     width = 8, height = 12, units = "in", res = 600)
p1
dev.off()



#.........................................................................####
# Supplementary Figure 1 ####
df <- readRDS(file = "C:/Users/julia/Documents/GitHub/normScore-development/Simulations/allSimulations.rds")
df <- df %>% filter(!(especificos %in% c(1,4,8))) # quitamos o 1 que ten control negativo en item 3, o 4 pq é moi reiterativo co 3 e non aporta nada, e o 8 porque pa valores extremos en item3 a situación ponse rara.
df$especificos <- factor(df$especificos, levels= c(0, 2:3, 5:7), labels = LETTERS[1:6])

## Opcion 1 ####
graficos <- sapply(unique(df$n_proteins), function(i){
  dfAux <- df %>% 
    filter(n_proteins == i) %>% 
    filter(especificos != "A") %>% # negative control
    # dplyr::select(Item0:normScore) %>%
    dplyr::select(normScore) %>%
    tidyr::pivot_longer(cols = colnames(.), names_to = "Metrics", values_to = "Cor")
  plotBoxChulo(
    data = dfAux, 
    varResumen = "Cor", 
    varGrupo = "Metrics", 
    tituloY = "Kendall's Correlation",
    color = "#2FB2AD",
    anguloX = 65,
    sizeDots = 0.6,
    colorMean = "#C02A4E",
    titulo = paste0("No. of proteins: ", i),
    sizeLetra = 10)
  
})

graficos2 <- sapply(unique(df$n_per_group), function(i){
  dfAux <- df %>% 
    filter(n_per_group == i) %>% 
    filter(especificos != "A") %>% # negative control
    dplyr::select(Item0:normScore) %>%
    tidyr::pivot_longer(cols = colnames(.), names_to = "Metrics", values_to = "Cor")
  plotBoxChulo(
    data = dfAux, 
    varResumen = "Cor", 
    varGrupo = "Metrics", 
    tituloY = "Kendall's Correlation",
    color = "#2E2753",
    anguloX = 65,
    sizeDots = 0.6,
    colorMean = "#C02A4E",
    titulo = paste0("N per group: ", i),
    sizeLetra = 10)
  
})

## Opcion 2 ####

### All error magnitudes ####
dfAux <- df %>% 
  filter(especificos != "A") %>% # negative control
  dplyr::select(normScore, n_proteins, n_per_group) 
dfAux$n_per_group <- factor(dfAux$n_per_group)
dfAux$n_proteins <- factor(dfAux$n_proteins)

pN <- plotBoxChulo(
  data = dfAux, 
  varResumen = "normScore", 
  varGrupo = "n_per_group", 
  tituloY = "Kendall's Correlation",
  color = "#2E2753",
  anguloX = 0,
  widthDots = 0.1,
  sizeDots = 1, 
  colorMean = "#C02A4E",
  tituloX = "Number of samples per group",
  titulo = "NormScore (global)",
  # titulo = "NormScore (global) - All error magnitudes",
  sizeLetra = 14, faceX = "bold") + 
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(hjust = 0.5), 
        plot.margin = margin(t = 30, r = 10, b = 10, l = 40), 
        plot.title = element_text(vjust = 8, face = "bold"))
pN
pP <- plotBoxChulo(
  data = dfAux, 
  varResumen = "normScore", 
  varGrupo = "n_proteins", 
  tituloY = "Kendall's Correlation",
  color = "#2FB2AD",
  anguloX = 0,
  widthDots = 0.1,
  sizeDots = 1, 
  colorMean = "#C02A4E",
  tituloX = "Number of proteins",
  titulo = "NormScore (global)",
  # titulo = "NormScore (global) - All error magnitudes",
  sizeLetra = 14, faceX = "bold") + 
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(hjust = 0.5), 
        plot.title = element_text(vjust = 8, face = "bold"),
        plot.margin = margin(t = 30, r = 10, b = 10, l = 40))
pP

ggpubr::ggarrange(plotlist = list(pN, pP), ncol = 2, nrow = 1, labels = "AUTO")


### All error magnitudes ####
dfAux <- df %>% 
  filter(especificos == "B") %>%
  filter(especificos != "A") %>% # negative control
  dplyr::select(normScore, n_proteins, n_per_group) 
dfAux$n_per_group <- factor(dfAux$n_per_group)
dfAux$n_proteins <- factor(dfAux$n_proteins)

pN_B <- plotBoxChulo(
  data = dfAux, 
  varResumen = "normScore", 
  varGrupo = "n_per_group", 
  tituloY = "Kendall's Correlation",
  color = "#2E2753",
  anguloX = 0,
  widthDots = 0.1,
  sizeDots = 1, 
  colorMean = "#C02A4E",
  tituloX = "Number of samples per group",
  # titulo = "NormScore (global) - Error magnitude B",
  titulo = "",
  # titulo = "NormScore (global)",
  sizeLetra = 14, faceX = "bold") + 
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(hjust = 0.5), 
        plot.margin = margin(t = 25, r = 10, b = 10, l = 40))
pN_B
pP_B <- plotBoxChulo(
  data = dfAux, 
  varResumen = "normScore", 
  varGrupo = "n_proteins", 
  tituloY = "Kendall's Correlation",
  color = "#2FB2AD",
  anguloX = 0,
  widthDots = 0.1,
  sizeDots = 1, 
  colorMean = "#C02A4E",
  tituloX = "Number of proteins",
  titulo = "",
  # titulo = "NormScore (global) - Error magnitude B",
  sizeLetra = 14, faceX = "bold") + 
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(hjust = 0.5), 
        plot.margin = margin(t = 25, r = 10, b = 10, l = 40))
pP_B

ggpubr::ggarrange(plotlist = list(pN_B, pP_B), ncol = 2, nrow = 1, labels = "AUTO")


### All error magnitudes ####
dfAux <- df %>% 
  filter(especificos == "F") %>%
  filter(especificos != "A") %>% # negative control
  dplyr::select(normScore, n_proteins, n_per_group) 
dfAux$n_per_group <- factor(dfAux$n_per_group)
dfAux$n_proteins <- factor(dfAux$n_proteins)

pN_F <- plotBoxChulo(
  data = dfAux, 
  varResumen = "normScore", 
  varGrupo = "n_per_group", 
  tituloY = "Kendall's Correlation",
  color = "#2E2753",
  anguloX = 0,
  widthDots = 0.1,
  sizeDots = 1, 
  colorMean = "#C02A4E",
  tituloX = "Number of samples per group",
  titulo = "",
  # titulo = "NormScore (global) - Error magnitude F",
  sizeLetra = 14, faceX = "bold") + 
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(hjust = 0.5), 
        plot.margin = margin(t = 25, r = 10, b = 10, l = 40))
pN_F
pP_F <- plotBoxChulo(
  data = dfAux, 
  varResumen = "normScore", 
  varGrupo = "n_proteins", 
  tituloY = "Kendall's Correlation",
  color = "#2FB2AD",
  anguloX = 0,
  widthDots = 0.1,
  sizeDots = 1, 
  colorMean = "#C02A4E",
  tituloX = "Number of proteins",
  titulo = "",
  # titulo = "NormScore (global) - Error magnitude F",
  sizeLetra = 14, faceX = "bold") + 
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(hjust = 0.5), 
        plot.margin = margin(t = 25, r = 10, b = 10, l = 40))
pP_F


### All together ####
listaPlots <- list(pN, pP, pN_B, pP_B, pN_F, pP_F)
labels <- c("Error magnitudes B to F", "", "Only error magnitude B", "", 
            "Only error magnitude F", "")
ps1 <- ggpubr::ggarrange(plotlist = listaPlots, ncol = 2, nrow = 3, labels = labels, 
                         label.x = -0.1, label.y = 0.925)
ps1

ggsave(filename = paste0(outputDir, "SupplementaryFigure1.svg"), 
       plot = ps1, width = 10, height = 12)




#.........................................................................####
# Figure 3G ####
## Loading results ####
inputDir <- "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/AssessmentFiles/"
resAllGS <- readRDS(file = paste0(inputDir, "3_Reviewer1_results_100dataset.rds")) # Manual ranking
resNS <- readRDS(file = paste0(inputDir, "3_normScore_results_100datasets.rds")) # Originales
resAllGSR2 <- readRDS(file = paste0(inputDir, "3_Reviewer2_results_100dataset.rds")) # Manual ranking

source(file = "../OneHundredDatasets/3_AssessmentFunctions.R")
resAllNS <- formatNormScoreResults(finalList = resNS, item0 = T)
resAllNSR2 <- resAllNS


## Estimating data ####
# R1 vs normScore
jacIndR1NS <- sapply(names(resAllGS$globalBests), function(i) 
  jaccardTopSet(
    rankedMethods = resAllNS$resGlobal[[i]], 
    referenceSet = resAllGS$globalBests[[i]]
  )
)

# R2 vs normScore
jacIndR2NS <- sapply(names(resAllGSR2$globalBests), function(i) 
  jaccardTopSet(
    rankedMethods = resAllNS$resGlobal[[i]], 
    referenceSet = resAllGSR2$globalBests[[i]]
  )
)

# R1 vs R2
jacIndR1R2 <- sapply(names(resAllGSR2$globalBests), function(i) 
  jaccardIndex(
    setA = resAllGS$globalBests[[i]], 
    setB = resAllGSR2$globalBests[[i]]
  )
)

## All together
dfJI <- data.frame(
  Datasets = names(jacIndR1NS), 
  R1NS = jacIndR1NS[names(jacIndR1NS)],
  R2NS = jacIndR2NS[names(jacIndR1NS)],
  R1R2 = jacIndR1R2[names(jacIndR1NS)]
)

# Format longer
dfJaccard <- dfJI %>%
  tidyr::pivot_longer(
    cols = -Datasets,
    names_to = "Comparison", 
    values_to = "Jaccard index"
  )

dfJaccard$Comparison <- factor(
  dfJaccard$Comparison,
  levels = c(
    "R1NS", 
    "R2NS", 
    "R1R2"
  ),
  labels = c(
    "NormScore vs Review 1",
    "NormScore vs Review 2",
    "Review 1 vs Review 2"
  )
)

# Count by jaccard value

dfJaccardCount <- dfJaccard %>%
  count(Comparison, `Jaccard index`) %>%
  mutate(`Jaccard index` = round(`Jaccard index`, 2)) %>%
  tidyr::complete(
    Comparison,
    `Jaccard index`,
    # `Jaccard index` = jaccardLevels,
    fill = list(n = 0)
  )


## Figure ####

### Option A ####
pA <- ggplot(dfJaccardCount, aes(x = factor(`Jaccard index`), y = n, 
                                 fill = Comparison, colour = Comparison)) +
  geom_col(width = 0.75, alpha = 1, linewidth = 0.5) +
  facet_wrap(~ Comparison, 
             ncol = 1,
             # axes = "all",
             strip.position = "right"
             # nrow = 3, 
             # axes = "all"
  ) +
  scale_fill_manual(values = c("darkgrey","royalblue4", "#A80729"))+
  scale_colour_manual(values = c("darkgrey","royalblue4",  "#A80729"))+
  labs(
    x = "Jaccard index",
    y = "Number of datasets",
    title = NULL
  ) +
  # ylim(c(0, 60))+
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "bold"),
    strip.background = element_rect(
      fill = "#F9F9F9",
      colour = "black",
      linewidth = 0.8),
    strip.placement = "outside",
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5) 
    # axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
    # axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black")
  )
pA
ggsave(filename = paste0(outputDir, "Figure6A.svg"), 
       plot = pA, width = 5, height = 10)


### Option B ####
pB <- ggplot(dfJaccardCount, aes(x = factor(`Jaccard index`), y = n, 
                                 fill = Comparison, colour = Comparison)) +
  geom_col(width = 0.75, alpha = 0.6) +
  facet_wrap(~ Comparison, 
             nrow = 3,
             axes = "all"
  ) +
  scale_fill_manual(values = c("darkgrey", "darkorange", "royalblue4"))+
  scale_colour_manual(values = c("darkgrey", "darkorange", "royalblue4"))+
  labs(
    x = "Jaccard index",
    y = "Number of datasets",
    title = NULL
  ) +
  # ylim(c(0, 60))+
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(
      size = 14,
      face = "bold"),
    # plot.title = element_text(face = "bold", hjust = 0.5),
    axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black")
  )
pB
ggsave(filename = paste0(outputDir, "Figure6B.svg"), 
       plot = pB, width = 6, height = 10)



### Option C ####

comparisonsList <- list(
  c("Review 1 vs Review 2", "NormScore vs Review 1"),
  c("Review 1 vs Review 2", "NormScore vs Review 2"),
  c("NormScore vs Review 1", "NormScore vs Review 2")
)


id_col <- "Datasets"
stopifnot(id_col %in% names(dfJaccard))

testsJaccard <- do.call(rbind, lapply(comparisonsList, function(grupos) {
  g1 <- dfJaccard[dfJaccard$Comparison %in% grupos[1], ]
  g2 <- dfJaccard[dfJaccard$Comparison %in% grupos[2], ]
  
  # Comprueba que hay un único valor por dataset en cada grupo.
  stopifnot(
    !anyNA(g1[[id_col]]), !anyNA(g2[[id_col]]),
    !anyDuplicated(g1[[id_col]]), !anyDuplicated(g2[[id_col]])
  )
  
  # Alinea los valores de ambos grupos por identificador de dataset.
  ids <- intersect(g1[[id_col]], g2[[id_col]])
  x <- g1[["Jaccard index"]][match(ids, g1[[id_col]])]
  y <- g2[["Jaccard index"]][match(ids, g2[[id_col]])]
  completos <- complete.cases(x, y)
  
  data.frame(
    group1 = grupos[1],
    group2 = grupos[2],
    n_pairs = sum(completos),
    p = stats::wilcox.test(x[completos], y[completos], paired = TRUE)$p.value
  )
}))

testsJaccard$p.adj <- stats::p.adjust(testsJaccard$p, method = "BH")
testsJaccard$p.adj.signif <- as.character(stats::symnum(
  testsJaccard$p.adj,
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  symbols = c("****", "***", "**", "*", "ns")
))
testsJaccard$y.position <- seq(1.06, 1.28, length.out = nrow(testsJaccard))

testsJaccard


coloresBP <- c("#00204DFF" , "#707173FF", "#F7D03CFF")

pC <- ggplot(dfJaccard, aes(x = Comparison, y = `Jaccard index`, 
                      fill = Comparison, colour = Comparison)) +
  geom_boxplot(
    width = 0.35,
    outlier.shape = NA,
    alpha = 0.3, 
    linewidth = 0.4, whisker.linewidth = 0.4, median.linewidth = 0.7
  ) +
  scale_fill_manual(values = coloresBP) + #c("darkgrey","royalblue4", "#A80729"))+
  scale_colour_manual(values = coloresBP) +#c("darkgrey","royalblue4",  "#A80729"))+
  ggrastr::rasterize(
    geom_jitter(
      width = 0.10,
      size = 1.8,
      alpha = 0.65
    ), dpi = 400) +
  # ggpubr::stat_compare_means(
  #   comparisons = comparisonsList,
  #   method = "wilcox.test",
  #   paired = TRUE,
  #   label = "p.signif"
  # ) +
  ggpubr::stat_pvalue_manual(
    testsJaccard,
    label = "p.adj.signif",
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1.4),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = NULL,
    y = "Jaccard index",
    title = "Global assessment"
  ) +
  ggpubr::theme_pubr(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
pC


ggsave(filename = paste0(outputDir, "Figure6G.svg"), 
       plot = pC, width = 6.5, height = 5)





#.........................................................................####
# Figure 4  ####

## Figure 4 B ####
lista <- readRDS(
  "../OneHundredDatasets/AssessmentFiles/4_nonRedundancyResults.rds"
)
pH1 <- pheatmap::pheatmap(
  lista$dfCor, 
  cluster_cols = TRUE,  
  cluster_rows = TRUE,
  fontsize = 14,
  fontsize_row = 16,
  fontsize_col = 16,
  border_color = "white",
  angle_col = 45,
  treeheight_row = 20,
  treeheight_col = 20,
  color = grDevices::colorRampPalette(rev(viridis::inferno(n =5)))(22), # Top 1
  breaks = seq(0, 1, 0.05),
  legend_breaks = seq(0, 1, 0.5),
  legend_labels = seq(0, 1, 0.5)
)


## Figure 4 A ####
# Long format by each matrix
corList <- lista$corList
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


outputDir <- "G:/Mi unidad/_General/DoctoradoIndustrial2022/Publicaciones/NORMALIZACIONES/Paper_normScore/FiguresNew/"

pByPairs <- ggplot(corLong, aes(x = pair, y = correlation, color = var1)) +
  geom_boxplot(aes(fill = var1), colour = "black", alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    x = "All pairs of items",
    y = "Kendall's correlation", 
    fill = NULL
  ) +
  scale_fill_manual(values = viridis::inferno(n =6)) +
  guides(fill = guide_legend(nrow = 1)) +
  # scale_colour_manual(values = viridis::inferno(n =6)) +
  theme(
    # panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom", 
    legend.text = element_text(size = 14)
  )
pByPairs
# ggsave(filename = paste0(outputDir, "NewFigure4A.svg"), 
#        plot = pByPairs, width = 15, height = 6)

## Figure 4 C ####
res <- readRDS(file = "../OneHundredDatasets/AssessmentFiles/4_rankStability_LOO_normalization_plots.rds")

pmatrix <- ggpubr::ggarrange(plotlist = sapply(res, "[[", 4), nrow = 2, ncol = 4, 
                             common.legend = TRUE)
pmatrix

ggsave(filename = paste0(outputDir, "Figure4C.svg"), 
       plot = pmatrix, width = 9, height = 5.5)

plegend <- sapply(res, "[[", 4)[[1]] + theme(legend.position = "right")
ggsave(filename = paste0(outputDir, "Figure4C_Legend.svg"), 
       plot = plegend, width = 8, height = 8)




#.........................................................................####
#  Figure 5 ####

# colorVector <- c("#762A83", "#F7F7F7", "#1B7837")
corPlot <- function(
    corMatrix, 
    varNames,
    colorVector =  c("firebrick", "white", "cadetblue4"),
    sizeLetra = 18, 
    sizeInnerText = 5, 
    titulo = NULL, 
    legPos = "bottom"
    ){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  dfCor <- as.data.frame(corMatrix) |>
    mutate(var1 = rownames(corMatrix)) |>
    pivot_longer(
      -var1,
      names_to = "var2",
      values_to = "correlation"
    )
  dfCor <- as.data.frame(dfCor)
  
  dfCor$var1 <- gsub(pattern = "RankNSNS", replacement = "normScore ", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "Rank", replacement = "", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "NS", replacement = "normScore", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "Human", replacement = "MAE Human", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "MaeYeast", replacement = "Yeast", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "Yeast", replacement = "MAE Yeast", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "Ecoli", replacement = "MAE Ecoli", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "upsMAE", replacement = "MAE UPS", x = dfCor$var1)
  dfCor$var1 <- gsub(pattern = "Y1vsY4", replacement = "Y1-Y4", x = dfCor$var1, fixed = T)
  dfCor$var1 <- gsub(pattern = "Y1vsY2", replacement = "Y1-Y2", x = dfCor$var1, fixed = T)
  dfCor$var1 <- gsub(pattern = "Y1vsY3", replacement = "Y1-Y3", x = dfCor$var1, fixed = T)
  
  dfCor$var2 <- gsub(pattern = "RankNSNS", replacement = "normScore ", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "Rank", replacement = "", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "NS", replacement = "normScore", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "Human", replacement = "MAE Human", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "MaeYeast", replacement = "Yeast", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "Yeast", replacement = "MAE Yeast", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "Ecoli", replacement = "MAE Ecoli", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "upsMAE", replacement = "MAE UPS", x = dfCor$var2)
  dfCor$var2 <- gsub(pattern = "Y1vsY4", replacement = "Y1-Y4", x = dfCor$var2, fixed = T)
  dfCor$var2 <- gsub(pattern = "Y1vsY2", replacement = "Y1-Y2", x = dfCor$var2, fixed = T)
  dfCor$var2 <- gsub(pattern = "Y1vsY3", replacement = "Y1-Y3", x = dfCor$var2, fixed = T)
  

  
  dfCor <- dfCor |>
    mutate(
      var1 = factor(var1, levels = rev(varNames)),
      var2 = factor(var2, levels = varNames),
      rowId = match(as.character(var1), rev(varNames)),
      colId = match(as.character(var2), varNames)
    ) |>
    filter(rowId + colId <= length(varNames))
  
  
  
  ggplot(
    dfCor,
    aes(
      x = var2,
      y = var1,
      fill = correlation
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 1
    ) +
    geom_text(
      aes(label = sprintf("%.2f", correlation)),
      size = sizeInnerText
    ) +
    scale_fill_gradient2(
      low = colorVector[1],
      mid = colorVector[2],
      high = colorVector[3],
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL, 
      subtitle = titulo,
      fill = "Kendall correlation"
    ) +
    theme_minimal(base_size = sizeLetra) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        colour = "black"),
      axis.text.y = element_text(colour = "black"),
      legend.position = legPos, 
      fill = guide_colorbar(
        barwidth = unit(10, "cm"),
        barheight = unit(0.6, "cm")
      )
    )
}

## A - DIA PXD028735 ####
varNames <- c("normScore", "MAE", "AUC", "MAE Human", "MAE Yeast", "MAE Ecoli")
colorVector <- c("firebrick", "white", "cadetblue4")
allResults <- readRDS(file = "../Benchmarking/Data/PXD028735/PXD028735_DIA_Complete_Analysis.rds")
corrplot::corrplot.mixed(
  allResults$correlationResults, upper = "ellipse")
corMatrix <- allResults$correlationResults

pCor <- corPlot(corMatrix = corMatrix, colorVector = colorVector, sizeInnerText = 5, 
        titulo = "PXD028735 - DIA", legPos = "none", varNames = varNames)
ggsave(
  filename = paste0(outputDir, "Figure5A.svg"), 
  plot = pCor, width = 4.5, height = 4.3)

## B - DDA PXD028735 ####
allResults <- readRDS(file = "../Benchmarking/Data/PXD028735/PXD028735_DDA_Complete_Analysis.rds")
corMatrix <- allResults$correlationResults
pCor <- corPlot(corMatrix = corMatrix, colorVector = colorVector, sizeInnerText = 5, 
        titulo = "PXD028735 - DDA", legPos = "right", varNames = varNames)
ggsave(
  filename = paste0(outputDir, "Figure5b.svg"), 
  plot = pCor, width = 7.5, height = 4.3)


## C - FC=10 PXD001819 ####
varNames <- c("normScore", "MAE", "AUC", "MAE Yeast", "MAE UPS")
res <- readRDS( file = "../Benchmarking/Data/PXD001819/completeResult50vs5_MaxQuant.rds")
dfCor <- res$dfCorKendall

colorVector <-  c("firebrick", "white", "cadetblue4") #c("salmon2", "white", "royalblue3")
pCor <- corPlot(corMatrix = dfCor, colorVector = colorVector, sizeInnerText = 5, 
                titulo = "PXD001819 - 50vs5 (FC =10)", legPos = "none", 
                varNames = varNames)
pCor
ggsave(
  filename = paste0(outputDir, "Figure5C_.svg"), 
  plot = pCor, width = 4, height = 3.7)

## D - FC=2 PXD001819 ####
res <- NULL
res <- readRDS( file = "../Benchmarking/Data/PXD001819/completeResult25vs12.5_MaxQuant.rds")
dfCor <- res$dfCorKendall
pCor <- corPlot(corMatrix = dfCor, colorVector = colorVector, sizeInnerText = 5, 
                titulo = "PXD001819 - 25vs12.5 (FC =12.5)", legPos = "right", 
                varNames = varNames)
pCor
ggsave(
  filename = paste0(outputDir, "Figure5D_.svg"), 
  plot = pCor, width = 7, height = 3.7)


## F - PXD039399 ####
varNames <- c("normScore", "MAE", "normScore Y1-Y2", "normScore Y1-Y3", "normScore Y1-Y4")
colorVector <-  c("firebrick", "white", "cadetblue4") #c("deeppink4", "white", "deepskyblue3")
res <- NULL
listaCompleta <- readRDS(file = "../Benchmarking/Data/PXD039399/completeResult.rds")
dfCor <- listaCompleta$dfCorKendall
dfCor <- dfCor[c(2, 1, 3:5), c(2, 1, 3:5)]
pCor <- corPlot(corMatrix = dfCor, colorVector = colorVector, sizeInnerText = 5, 
                titulo = "PXD039399", legPos = "bottom", 
                varNames = varNames)
pCor
ggsave(
  filename = paste0(outputDir, "Figure5F_LEGEND.svg"), 
  plot = pCor, width = 7.4, height = 4)

## E - PXD039399 ####

data <- listaCompleta$normData$Log

df <- as.data.frame(apply(data, 2, sum, na.rm = T))
colnames(df) <- "Intensity"
df$Samples <- factor(rownames(df), levels = rownames(df), ordered = T)
df$Groups <- rep(c("Y1H", "Y2H", "Y3H", "Y4H"), each = 3)
rownames(df) <- 1:nrow(df)
rownames(df) <- df$Samples
TIplot <- ggplot2::ggplot(
  df, 
  ggplot2::aes(x = Samples, y = Intensity, fill = Groups)) + 
  ggplot2::geom_bar(stat = "identity") +
  scale_fill_manual(values = viridis::mako(n=5)[1:4])+
  ggplot2::theme_minimal(base_size = 14) + 
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
    axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
    axis.text.x = ggplot2::element_text(angle = 75, hjust = 1)) + 
  ggplot2::ylab("Total Intensity") + 
  ggplot2::xlab("Samples") +
  labs(title = "Sum of total intensity (item 0)")
TIplot




df_long <- data %>% 
  as.data.frame() |> 
  tidyr::pivot_longer(cols = tidyr::everything(), names_to = "Variable", 
                      values_to = "Valor")
df_long$Variable <- factor(df_long$Variable, levels = varResumen)
df_long$Groups <- sapply(strsplit(x = df_long$Variable, split = "_", fixed = T), "[[", 1)
pBMulti <- ggplot2::ggplot(
  df_long, 
  ggplot2::aes(
    x = Variable, 
    y = Valor, 
    fill = Groups)) + 
  ggplot2::geom_boxplot(
    width = 0.5, alpha = 0.75, 
  ) + 
  ggplot2::scale_fill_manual(values = viridis::mako(n=5)[1:4])+
  ggplot2::labs(title = "Distribution of log intensity (item 6)", 
                x = NULL, y = "Log(Intensity)") +
  ggplot2::theme_minimal(base_size = 12) + 
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.5,  colour = "black"), 
    axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = ggplot2::unit(0.2,"cm"), 
    axis.text.x = ggplot2::element_text(hjust = 1, angle = 65), 
    panel.grid.major.x = ggplot2::element_blank(), 
    panel.grid.minor.x = ggplot2::element_blank(), 
    plot.title = ggplot2::element_text(hjust = 0.5))
pBMulti


#.........................................................................####
# Supplementary Figure 2 ####
allOut <- readRDS(
  file = "C:/Users/julia/Documents/GitHub/normScore-development/OneHundredDatasets/AssessmentFiles/3_NormScoreAssesment_models_Covars.rds")
# names(allOut) <- c(names(allOut)[1:5, "Species"])
# allOut$Specie
covariables <- c("Acquisition", "Instrument", "CellCulture", "Software", 
                 "Species")
etiquetaP <- allOut$Etiquetas
allOut$Etiquetas <- NULL
listDfMedias <- sapply(allOut, "[[", 2, simplify = F)


graficosCovariables <- purrr::map(
  covariables,
  function(covariable) {
    dfMedias <- listDfMedias[[covariable]]
    ggplot(
      dfMedias,
      aes(
        x = Comparison,
        y = emmean,
        color = .data[[covariable]],
        shape = .data[[covariable]],
        group = .data[[covariable]]
      )
    ) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(
        aes(
          ymin = lower.CL,
          ymax = upper.CL
        ),
        width = 0.1
      ) +
      scale_color_manual(values = RColorBrewer::brewer.pal(
        n = length(unique(dfMedias[,covariable])), 
        name = "Set1")) +
      # scale_color_manual(values = colores[1:length(unique(dfComp[,covariable]))]) +
      labs(
        title = covariable,
        subtitle = etiquetaP[covariable],
        x = NULL,
        y = "Jaccard index",
        color = covariable
      ) +
      ylim(c(0, 1))+
      theme_minimal(base_size = 12) + 
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "top"
      )
  }
)

pAll <- ggpubr::ggarrange(plotlist = graficosCovariables, ncol = 2, nrow = 3)
pAll
ggsave(filename = paste0(outputDir, "SupplementaryFigure2.svg"), 
       plot = pAll, width = 10, height = 15)





#.........................................................................####
# Supplementary Table 3 ####


## Items 0, 1, 5, 6 ####
# errors: 2, 3, 5:7
### sample_shift_sd ####
item0156_sample_shift_sd <- list(
  # seq(0, 0.05, 0.00625),
  seq(0, 0.1, 0.0125),
  seq(0, 0.5, 0.0625), 
  # seq(0, 0.5, 0.0625),
  seq(0.5, 1.5, 0.125),
  seq(0.5, 1.5, 0.125),
  seq(1.5, 3.5, 0.25))
  # seq(1.5, 3.5, 0.25))
names(item0156_sample_shift_sd) <- LETTERS[2:6]

### sample_shift_cap ####
offset <- c(0.05, 0.05, 0.05, -0.05, 0.05)
item0156_sample_shift_cap <- sapply(1:length(item0156_sample_shift_sd), function(x){
  item0156_sample_shift_sd[[x]] + offset[x]
}, simplify = F, USE.NAMES = T)
names(item0156_sample_shift_cap) <- LETTERS[2:6]



df <- as.data.frame(do.call(rbind, item0156_sample_shift_sd))
rownames(df) <- paste0("Item0156_shiftSd_", rownames(df))
df2 <- as.data.frame(do.call(rbind, item0156_sample_shift_cap))
rownames(df2) <- paste0("Item0156_shiftCap_", rownames(df2))

df <- rbind(df, df2)
colnames(df) <- paste0("Perturbation ", 1:9)






## Item 2 ####

### rho_between ####
item2_values_rho_between <- list(
  # seq(0, 0.05, 0.00625),
  seq(0, 1.25, 0.15625), 
  seq(0, 1, 0.125), 
  # seq(0, 0.75, 0.09375), 
  seq(0, 0.5, 0.0625),
  seq(0, 0.5, 0.0625),
  seq(0, 0.25, 0.03125))
  # seq(0, 0.05, 0.00625)) 
names(item2_values_rho_between) <- LETTERS[2:6]

### rho_within ####
item2_values_rho_within <- list(
  # seq(3.5, 1.5, -0.25),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  # rev(seq(0, 1, 0.125)),
  rev(seq(0, 1, 0.125)),
  rev(seq(0, 0.75, 0.09375)),
  rev(seq(0, 1, 0.125)))
  # rev(seq(0, 1, 0.125)))
names(item2_values_rho_within) <- LETTERS[2:6]


df1 <- as.data.frame(do.call(rbind, item2_values_rho_between))
rownames(df1) <- paste0("Item2_rhoBetween_", rownames(df1))
df2 <- as.data.frame(do.call(rbind, item2_values_rho_within))
rownames(df2) <- paste0("Item2_rhoWithin_", rownames(df2))

df3 <- rbind(df1, df2)
colnames(df3) <- paste0("Perturbation ", 1:9)
df <- rbind(df, df3)
df




## Item 4 ####

### sample_sd_cap ####
item4_sample_sd_cap <- list(
  # seq(0, 0.5, 0.0625),
  seq(0, 1, 0.125),
  seq(0, 1.5, 0.1875),
  # seq(0.5, 2, 0.1875),
  rep(1, 9),
  rep(1, 9),
  rep(1, 9)
  # rep(1, 9)
)
names(item4_sample_sd_cap) <-LETTERS[2:6]

### sample_sd_cstrength ####
item4_sample_sd_strength <- list( # rep(c(1, 2), each = 4) 
  # rep(1, 9),
  rep(1, 9),
  rep(1, 9),
  # rep(1, 9),
  seq(0, 0.25, 0.03125)*(-1),
  seq(0, 0.5, 0.0625)*(-1),
  seq(0, 1, 0.125)*(-1)
  # seq(0, 1.5, 0.1875)*(-1)
) # 2
names(item4_sample_sd_strength) <- LETTERS[2:6]


### sample_sd_rho ####
item4_sample_sd_rho <- c(1.5, 1.5, 0.9, 0.9, 0.9, 0.9)
df4 <- t(data.frame(
  Item4_sdRho_B = rep(1.5, 9),
  Item4_sdRho_C = rep(1.5, 9),
  Item4_sdRho_D = rep(0.9, 9),
  Item4_sdRho_E = rep(0.9, 9),
  Item4_sdRho_F = rep(0.9, 9)
))
colnames(df4) <- paste0("Perturbation ", 1:9)



df1 <- as.data.frame(do.call(rbind, item4_sample_sd_cap))
rownames(df1) <- paste0("Item4_sdCap_", rownames(df1))
df2 <- as.data.frame(do.call(rbind, item4_sample_sd_strength))
rownames(df2) <- paste0("Item4_sdStrength_", rownames(df2))

df3 <- rbind(df1, df2)
colnames(df3) <- paste0("Perturbation ", 1:9)
df <- rbind(df, df3)
df <- rbind(df, df4)



## Item 3 ####

### sigma_lo ####
item3_sigma_lo <- list(
  # rep(0.4, 9), # Control (-1)
  rep(0.4, 9), 
  rep(0.4, 9)*1.5, 
  # rep(0.4, 9), 
  rep(0.4, 9), 
  rev(seq(0, 1.5, 0.1875)), 
  rev(seq(0.4, 2, 0.2))
  # c(3, 2, 1, 0.5, 0.3, 0.3, 0.3, 0.3, 0.1) # solo forma
)

names(item3_sigma_lo) <- LETTERS[2:6]

### sigma_hi ####
item3_sigma_hi <- list(
  # rep(0.05, 9), # Control
  rep(0.05, 9), 
  rep(0.05, 9)*1.5, 
  # rep(0.05, 9), 
  rep(0.05, 9), 
  seq(0, 1.5, 0.1875)*(1.5),
  seq(0.4, 2, 0.2)*(1.5)
  # c(0.05, 0.1, 0.2, 0.4, 1.2, 1.2, 1.2, 1.9, 2.5)
)
names(item3_sigma_hi) <- LETTERS[2:6]

df1 <- as.data.frame(do.call(rbind, item3_sigma_hi))
rownames(df1) <- paste0("item3_sigmaHi_", rownames(df1))
df2 <- as.data.frame(do.call(rbind, item3_sigma_lo))
rownames(df2) <- paste0("item3_sigmaLo_", rownames(df2))

df3 <- rbind(df1, df2)
colnames(df3) <- paste0("Perturbation ", 1:9)
df <- rbind(df, df3)

### sample_sd_cap ####
item3_sample_sd_cap <- list(
  # seq(2.5, 1, -0.1875), # Control - orden inverso, ten que dar 0 
  seq(1, 2.5, 0.1875), 
  seq(1, 2.5, 0.1875),  
  # seq(0.5, 2, 0.1875),
  seq(0, 1, 0.125),
  seq(0, 1.5, 0.1875), 
  seq(0, 1.5, 0.1875) 
  # rep(0.2, 9)
  # 8 
)
names(item3_sample_sd_cap) <- LETTERS[2:6]

### sample_sd_strength ####
item3_sample_sd_strength = list(
  # rep(3, 9), 
  rep(3, 9), 
  rep(3, 9), 
  # rep(3, 9), 
  rep(3, 9), 
  rep(3, 9), 
  rep(2.5, 9) 
  # rep(0.8, 9)
)
names(item3_sample_sd_strength) <- LETTERS[2:6]

df1 <- as.data.frame(do.call(rbind, item3_sample_sd_cap))
rownames(df1) <- paste0("item3_sdCap_", rownames(df1))
df2 <- as.data.frame(do.call(rbind, item3_sample_sd_strength))
rownames(df2) <- paste0("item3_sdStrength_", rownames(df2))

df3 <- rbind(df1, df2)
colnames(df3) <- paste0("Perturbation ", 1:9)
df <- rbind(df, df3)


### sample_sd_rho ####
item3_sample_sd_rho = list(
  # rep(0, 9), 
  rep(0, 9), 
  rep(0, 9), 
  # rep(0, 9), 
  rep(0, 9), 
  rep(0, 9), 
  rep(0.2, 9)
  # rep(0.35, 9)
)
names(item3_sample_sd_rho) <- LETTERS[2:6]


df1 <- as.data.frame(do.call(rbind, item3_sample_sd_rho))
rownames(df1) <- paste0("item3_sdRho_", rownames(df1))
colnames(df1) <- paste0("Perturbation ", 1:9)
df <- rbind(df, df1)

metadata <- strsplit(rownames(df), split = "_", fixed = T)
df$Item <- sapply(metadata, "[[", 1)
df$Parameter <- sapply(metadata, "[[", 2)
df$ErrorMagnitude <- sapply(metadata, "[[", 3)
#


# ALL SCORE ####


sevList <- list(
  seq(0, 0.1, 0.0125),
  seq(0, 0.25, 0.03125),
  seq(0.25, 0.5, 0.03125),
  seq(0, 0.5, 0.0625),
  seq(0.25, 1, 0.09375),
  seq(0, 1, 0.125),
  seq(0.5, 1.5, 0.125),
  seq(1, 1.5, 0.0625)
)

sevList <- sevList[c(2:3, 5:7)]

df3 <- as.data.frame(do.call(rbind, sevList))
colnames(df3) <- paste0("Perturbation ", 1:9)
df3$Item <- "Global normScore"
df3$Parameter <- "Multiple"
df3$ErrorMagnitude <- LETTERS[2:6]

df <- rbind(df, df3)


writexl::write_xlsx(x = df, path = paste0(outputDir, "SupplementaryTable2.xlsx"))





