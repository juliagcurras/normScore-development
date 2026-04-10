################################################################################-
##################     Extract data for 100 dataset          ###################
################################################################################-


# Julia G Currás - 08/04/2026
rm(list=ls())
graphics.off()
setwd("H:/Mi unidad/Doctorado/SCORE/R")

# Librerías ####
library(dplyr)
library(ggplot2)


# Global vars ####
inDir <-  "C:/Users/julia/Documents/GitHub/normScore-development/goldStandard/ProcessedDatasets/"


# Functions ####

getkappaMatrix <- function(df){
  # parámetro "df": data frame con las variables categóricas de interés
  KappaMatrix <- apply(df, 2, function(i) apply(df, 2, function(j) { # "iterar" por columnas dos veces
    vcd::Kappa(table(i, j))[["Unweighted"]]["value"]
  }))
  return(KappaMatrix)
}


getSignProts <- function(
    listaResults, 
    dfRaw, # just for protein names
    limPval = 0.05,
    grupo1 = "Case", 
    grupo2 = "Control"){
  
  #----- Output 1. Get significant proteins for each method
  listaProtsDF <- lapply(listaResults, function(dfDE) {
    rownames(dfDE) <- rownames(dfRaw)
    datos <- dfDE %>% 
      # dplyr::filter(adj.P.Val < limPval) %>% 
      dplyr::filter(P.Val < limPval) %>% 
      tibble::rownames_to_column("Protein") %>% 
      as.data.frame
    datos
  })
  
  #List of DE proteins
  listaProts2 <- sapply(listaProtsDF, "[[", 1, simplify = F)
  
  #----- Output 2. Total number of DE by normalization
  totalProts <- sapply(listaProts2, length)
  percDAProts <- totalProts/nrow(dfRaw)
  
  # List of unique proteins
  listaProts3 <-  sort(unique(Reduce(c, listaProts2))) # vector of total unique prots
  
  #----- Output 3. Prots by normalization method - get df 0 and 1
  dftt <- data.frame(matrix(nrow = length(listaProts3), ncol = length(listaResults)))
  rownames(dftt) <- listaProts3
  colnames(dftt) <- names(listaResults)
  for (nome in names(listaProts2)){
    varBody <- vector()
    for (el in listaProts3){
      if (el %in% listaProts2[[nome]]){
        varBody <- c(varBody, "Si")
      } else {
        varBody <- c(varBody, "No")
      }
    }
    dftt[nome] <- varBody
  }
  dftt <- as.data.frame(dftt)
  
  #----- Output 4.  Common Proteins
  percDECommonProts <- sum(rowSums(dftt == "Si") == ncol(dftt))/nrow(dfRaw)
  # commonProts <- names(which(rowSums(dftt == "Si") == ncol(dftt)))
  
  #----- Output 5. Concordance (kappa index)
  prots <- rownames(dftt)
  protsTotal <- rownames(dfRaw)
  protsNosign <- protsTotal[!(protsTotal %in% prots)]
  
  # Joint information 
  dfAux <- data.frame(matrix(data = "No", nrow = length(protsNosign), ncol = ncol(dftt)))
  colnames(dfAux) <- colnames(dftt)
  dfTotal <- rbind(dftt, dfAux)
  rownames(dfTotal) <- c(rownames(dftt), protsNosign)
  dfTotal[,1:ncol(dfTotal)] <- lapply(colnames(dfTotal), 
                                      function(i) factor(dfTotal[,i], 
                                                         levels = c("Si", "No")))
  kappaIndex <- getkappaMatrix(dfTotal) # kappa estimation
  
  #---- Final: Return info
  return(list(
    percDAProts = percDAProts, 
    percDACommonProts = percDECommonProts, 
    dfProtsBymethod = dftt, 
    kappas = kappaIndex,
    dfFC = listaProtsDF
    )
    )
}




#..........................................................................####
# Checkpoint: GROUPS ####
# Selected Control and Case groups always appear in DM matrix?
df <- readxl::read_xlsx(
  path = "C:/Users/julia/Documents/GitHub/normScore-development/goldStandard/_Index.xlsx", 
  sheet = 1, 
  col_names = T)
allFiles <- df %>% 
  pull(ID) %>% paste0(., ".rds")

i <- "PXD056771"
resCheck <- sapply(allFiles, function(i){
  # dm <- as.data.frame(readRDS(file = paste0(inDir, i, ".rds"))[["dm"]])
  dm <- as.data.frame(readRDS(file = paste0(inDir, i))[["dm"]])
  # table(dm$Groups)
  # paste0(unique(dm$Groups), collapse = ";")
  gCaso <- df %>% 
    dplyr::filter(ID == gsub(pattern = ".rds", replacement = "", x = i)) %>% 
    dplyr::pull(Case)
  gControl <- df %>% 
    dplyr::filter(ID == gsub(pattern = ".rds", replacement = "", x = i)) %>% 
    dplyr::pull(Control)
  
  if (all(gCaso %in% unique(dm$Groups), gControl %in% unique(dm$Groups))){
    return("OK")
  } else if (all(!(gCaso %in% unique(dm$Groups)), gControl %in% unique(dm$Groups))){
    return("Lacking case")
  } else if (all(gCaso %in% unique(dm$Groups), !(gControl %in% unique(dm$Groups)))){
    return("Lacking control")
  } else {
    return("Any")
  }
}, simplify = T, USE.NAMES = T)
table(resCheck)
# resCheck[which(resCheck == "Any")]
# resCheck[which(resCheck == "Lacking case")]
# resCheck[which(resCheck == "Lacking control")]



#..........................................................................####
# Obter obxectos necesarios ####

## Cargar metadata ####
dfGrupo <- readxl::read_xlsx(
  path = "C:/Users/julia/Documents/GitHub/normScore-development/goldStandard/_Index.xlsx", 
  sheet = 1,
  col_names = T)
dfGrupo <- dfGrupo %>% 
  tibble::column_to_rownames("ID") %>%
  as.data.frame() %>%
  dplyr::select(Control, Case)

## Ejecutar toooodo ####
todo <- sapply(rownames(dfGrupo), function(idDataset){
  cat("Empezaaaando con dataset ", idDataset, "\n")
  # Data
  out <- readRDS(file = paste0(inDir, idDataset, ".rds"))
  dm <- out$dm
  dfRaw <- out$data
  gControl <- dfGrupo[idDataset, "Control"]
  gCaso <- dfGrupo[idDataset, "Case"]

  # DE analysis
  resComp <- lapply(out$listaNorm, Biomics::doTestT, dfGrupos = dm,
                    g1 = gControl, g2 = gCaso)
  
  # Compare DE results across normalizations
  resFinal <- getSignProts(
    listaResults = resComp, 
    dfRaw = dfRaw,
    limPval = 0.05, 
    grupo1 = gCaso, 
    grupo2 = gControl)
  
  # Add DE results to resFinal
  resFinal[["DEResults"]] <- resComp
  
  # Return!
  return(resFinal)
  },
simplify = F, USE.NAMES = T)

saveRDS(todo, file = "20260409_infoNeededForNormalizationComparison.rds")


#..........................................................................####
#--- Aggregating data across 100 dataset####

todo <- readRDS(file = "20260409_infoNeededForNormalizationComparison.rds")

  ## MATRIX 1: DAP ####
  # Para a matriz de Proteínas x normalizatión con DA ou non DA:

    ### 1. % DA proteins by normalization ####

# Cálculos
percDAProts <- sapply(todo, "[[", 1)
df1 <- as.data.frame(t(percDAProts)) 
res1 <- apply(df1, 2, Biostatech::getMeanIC)
dfRes1 <- sapply(res1, "[[", 1)
rownames(dfRes1) <- colnames(res1$Log$lista)
dfRes1

# Boxplot graphical representation
dfLong <- df1 %>% 
  tibble::rownames_to_column("ID") %>%
  as.data.frame() %>%
  dplyr::select(ID, everything()) %>%
  tidyr::pivot_longer(
    cols = colnames(df1), 
    names_to = "Normalizations", 
    values_to = "percDAP")

Biostatech::plotBox( #- Typical boxplot
  base = dfLong, 
  varResumen = "percDAP", 
  varGrupo = "Normalizations")$grafico


ggplot(dfLong, aes(x = Normalizations, y = percDAP)) + #- Paired boxplot
  geom_boxplot(
    aes(fill = Normalizations),
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_line(
    aes(group = ID),
    color = "grey70",
    alpha = 0.5,
    linewidth = 0.4
  ) +
  geom_jitter(
    aes(color = Normalizations),
    width = 0.08,
    size = 2,
    alpha = 0.8
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
    axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black")
  ) + 
  ggplot2::ylim(c(0, 1)) +
  ylab("% DAP")

resFT <- friedman.test(y = dfLong$percDAP, groups = dfLong$Normalizations, blocks = dfLong$ID)
resFT
resPWWT <- pairwise.wilcox.test(
  x = dfLong$percDAP, 
  g = dfLong$Normalizations, 
  paired = TRUE, 
  p.adjust.method = "holm")
resPWWT # VSN vs Log

ggpubr::ggviolin( #- Violin boxplot
  data = dfLong, 
  x = "Normalizations", 
  y = "percDAP",
  fill = "Normalizations",
  palette = viridis::viridis(n = length(unique(dfLong$Normalizations)), alpha = 0.1),
  alpha = 0.3,
  add = c("jitter", "boxplot"), 
  ylab = "% DAP") + 
  ggplot2::ylim(c(0, 1))
summary(dfLong$percDAP)

# Mean + IC => forest plot
Biostatech::plotForest(etiquetas = colnames(dfRes1),
                       estPunt = dfRes1["Media",],
                       LI = dfRes1[5,],
                       LS = dfRes1[6,],
                       vertical = T,
                       tituloX = "Media (IC%) porcentaje proteínas significativas")




    ### 2. % DA proteins comunes / total DE ####
propDACommonProts <- sapply(todo, "[[", 2)

# Graphical representation
dfPlot <- data.frame(
  ID = names(propDACommonProts), 
  Prop = unname(propDACommonProts)
)
dfPlot <- dfPlot[order(dfPlot$Prop), ]
dfPlot$ID <- factor(dfPlot$ID, levels = dfPlot$ID)

Biostatech::getNumTable(x = dfPlot[, "Prop", drop = F])
Biostatech::plotBarUnivar(var = dfPlot$ID, freqVar = dfPlot$Prop, anguloX = 65,
                          freqRel = T, limY = c(0,1), interact = F)$grafico
Biostatech::plotBox(base = dfPlot, varResumen = "Prop", interact = F)$grafico



    ### 3. Kappa index by pairs of normalizations ####
df3 <- sapply(todo, "[[", 4, simplify = F) 

# Quickly way to compute Kappa mean
tmp = Reduce('+', df3)
result = tmp/length(df3)
result <- as.matrix(result)
cores <- Biostatech::colorPalette(gradiente = T, removeWhite = F)[100:200] #110:1

corrplot::corrplot.mixed(result, is.corr = F, lower = "ellipse", upper = "number", 
                         order = "hclust", hclust.method = "single", 
                         col.lim = c(0, 1), title="Índice de Kappa", 
                         mar=c(0,0,1,0), lower.col = cores,
                         upper.col = cores)

# Getting IC for forest plot
listaDF <- lapply(df3, function(dat){ # symmetric matrix to df format for each dataset
  dat <- as.matrix(dat)
  dfFinal <- as.data.frame(as.table(dat)) # all pairs
  # # Execute for selecting unique pairs of comparison
  # dfFinal <- subset(
  #   x = dfFinal,
  #   subset = match(Var1, colnames(dat)) > match(Var2, colnames(dat))
  # )
})

df <- data.frame(matrix(data = NA, nrow = nrow(listaDF[[1]]), # from list of df to a unique df
                        ncol = 2 + length(listaDF)))
colnames(df) <- c("Var1", "Var2", names(listaDF))
df$Var1 <- listaDF[[1]][,1]
df$Var2 <- listaDF[[1]][,2]
for (i in 1:length(listaDF)){ # Just fill in with info from each dataset
  df_i <- listaDF[[i]]
  if(any((df_i[,1] != df$Var1), (df_i[,2] != df$Var2))){
    warning("algo va mal!")
    break
  }
  df[names(listaDF)[i]] <- df_i$Freq
}

res3 <- apply(df[,3:ncol(df)], 1, Biostatech::getMeanIC) # Computing Mean and IC!!
res3 <- sapply(res3, "[[", 1) # only table, not tablaFormato
res3 <- as.data.frame(t(res3[c(1, 5:6),])) # subsetting mean, LI and LS
colnames(res3) <- c("Media", "LI", "LS")
res3$Metodo1 <- df$Var1
res3$Metodo2 <- df$Var2
res3 <- res3[, c(4:5, 1:3)]

res3 <- res3[which(res3$Metodo1 != res3$Metodo2),] # removing pairs of the same normalization
res3$Metodo2 <- factor( # kepping the order
  res3$Metodo2, 
  levels = c("Log", "Median", "Quantile", "RLR", "VSN", "CyclicLoess", "TI", "Mean"),
  ordered = T)
res3 <- res3 %>% dplyr::arrange(Metodo2) # sorting

Biostatech::plotForest(
  etiquetas = paste0(res3$Metodo1, " - ", res3$Metodo2),
  estPunt = res3$Media,
  LI = res3$LI, #sizePoint = 2.5,
  LS = res3$LS, #sizeLetra = 75, 
  interact = F,
  tituloX = "Mean [IC95%] of Kappa index")

ggpubr::ggboxplot( #- Violin boxplot
  data = res3, 
  x = "Metodo2", 
  y = "Media",
  fill = "Metodo2",
  palette = viridis::viridis(n = length(unique(res3$Metodo2)), alpha = 0.1),
  alpha = 0.7,
  add = c("jitter"), 
  ylab = "Mean of concordance") + 
  ggplot2::ylim(c(0, 1))



  ## MATRIX 2: LOGFC ####
    # Para a matriz de proteínas x normalización con logFc
### (Getting data ####

listResComp <- sapply(todo, "[[", 6, simplify = F)  

dataFCMatrix <- lapply(listResComp, function(listDfDAA){
  # listDfDAA <- df4$PXD055210
  
  #.--- Format data properly ####
  dfFC <- sapply(listDfDAA, "[[", 1)
  rownames(dfFC) <- rownames(listDfDAA$Log)
  dfFC <- as.data.frame(dfFC)
  
  #.--- Statistics about logFC of each protein ####
  dfResumen <- as.data.frame(matrix(data = NA, ncol = 6, 
                                    nrow = nrow(dfFC)))
  rownames(dfResumen) <- rownames(dfFC)
  colnames(dfResumen) <- c("Media", "SD", "CV", "SameSign", 
                           "TotalPos", "TotalNeg")
  
  dfResumen$SameSign <- as.vector(apply(dfFC, 1, function(i) any(all(i > 0), all(i<0))))
  dfResumen$TotalPos <- as.vector(apply(dfFC, 1, function(i) sum(i > 0)))
  dfResumen$TotalNeg <- as.vector(apply(dfFC, 1, function(i) sum(i < 0)))
  dfResumen$Media <- as.vector(apply(dfFC, 1, function(i) mean(i, na.rm = T)))
  dfResumen$SD <- as.vector(apply(dfFC, 1, function(i) sd(i, na.rm =T)))
  dfResumen$CV <- (dfResumen$SD/abs(dfResumen$Media))*100
  table(dfResumen$SameSign)
  dfResumen$SameSignInt <- dfResumen$SameSign
  dfResumen$SameSignInt <- factor(dfResumen$SameSignInt,
                                  levels = c(F, T),
                                  labels = c("Different sign", "Equal sign"))
  
  dfAll <- merge(dfFC, dfResumen, by = 0)
  colnames(dfAll)[1] <- "Proteins"
  
  #.--- OUTPUT metrics ####
  # (1) percentage of prots with equal sign logFC
  porcEqualProts <- unname(prop.table(table(dfResumen$SameSignInt))["Equal sign"])
  
  # (2) logFC of prots with different sign
  rangoError <- dfResumen %>% filter(SameSign == F) %>% pull(Media) %>% range
  
  # (3) correlation between de variability of logFC by each method (CV) and the 
  # mean of the logFC.
  dfResumen$mediaAbs <- abs(dfResumen$Media)
  tipoCor <- ifelse(nortest::ad.test(dfResumen$CV)$p.value < 0.05, "spearman", "pearson")
  corCV_medialogFC <- cor(x = dfResumen$mediaAbs, y = dfResumen$CV, 
                          method = tipoCor, use = "complete.obs")
  # Biostatech::tableDescr(SameSignInt ~ mediaAbs + CV, dfResumen)$tablaFormato
  # Future possibility: difference relative mean between diff sign and equal sign
  
  #### Return ####
  aggregatedMetrics <- list(
    percEqualSign = porcEqualProts, 
    minMeanFCDiffSign = rangoError[1], 
    maxrangeMeanFCDiffSign = rangoError[2], 
    corCV_meanFC = corCV_medialogFC, 
    dfDetailed = dfAll
  )
  return(aggregatedMetrics)
  
}
)
### ) ####


    ### 4. % proteínas con igual logFC entre normalizations ####
equalProts <- sapply(dataFCMatrix, "[[", 1)
Biostatech::getNumTable(equalProts)$tabla
Biostatech::getMeanIC(equalProts)$tablaFormato


    ### 5. Correlación da media FC entre normalizacións vs CV ####
corCV_media <- sapply(dataFCMatrix, "[[", 4)
Biostatech::getMeanIC(corCV_media)$tablaFormato
Biostatech::getNumTable(corCV_media)$tablaFormato
# Canto máis altos son of valores de logFC, menor variabilidade hai entre normalizacións


    ### 6. Promedio rango de diferencias de signo entre normalizacion ####
minRango <- sapply(dataFCMatrix, "[[", 2)
maxRango <- sapply(dataFCMatrix, "[[", 3)
minRes <- Biostatech::getMeanIC(minRango)
minRes$tablaFormato
maxRes <- Biostatech::getMeanIC(maxRango)
maxRes$tablaFormato

dfPlot <- data.frame(
  ID = c(names(minRango), names(maxRango)), 
  Type = c(rep("Minimum", length(minRango)), rep("Maximum", length(maxRango))),
  Values = c(unname(minRango), unname(maxRango))
)

ggpubr::ggboxplot( #- Violin boxplot
  data = dfPlot, 
  x = "Type", 
  y = "Values",
  fill = "Type",
  palette = viridis::viridis(n = length(unique(dfPlot$Type)), alpha = 0.1),
  alpha = 0.7,
  add = c("jitter"), 
  ylab = "Range values") + 
  ggplot2::ylim(c(-4, 4))



    ### FIGURA CONXUNTA ####
# Extract info
listPlot <- sapply(dataFCMatrix, "[[", 5, simplify = F)
listPlot <- sapply(listPlot, function(df) df[, c("Media", "SameSignInt")], 
                   simplify = F)
listPlot <- sapply(names(listPlot), function(nome) {
  df <- listPlot[[nome]]
  df$Dataset <-  nome
  df
  }, simplify = F)

# Table info
dfPlot <- do.call(rbind, listPlot)
colnames(dfPlot)
ordenDatasets <- dfPlot %>% 
  group_by(Dataset) %>%
  summarise(maximo = max(abs(Media))) %>%
  arrange(maximo) %>% 
  pull(Dataset)

dfPlot$Dataset <- factor(dfPlot$Dataset, levels = rev(ordenDatasets))
dfPlot$ID <- as.numeric(as.factor(dfPlot$Dataset))
table(dfPlot$ID)
dfPlot$SameSignInt <- factor(dfPlot$SameSignInt, 
                             levels = c("Equal sign", "Different sign"))

# Plot
sizeLetra <- 14
Biostatech::plotScatter(
  base = dfPlot, 
  varX = "Media", 
  varY = "ID",
  varGrupo = "SameSignInt",
  adjustLine = F, 
  color = Biostatech::colorPalette(n = 2), 
  interact = F, 
  sizeDots = 1, 
  tituloX = "Mean FC between normalizations", 
  tituloY = "Dataset", 
  tituloLeyenda = "Proteins with..."
  )$grafico + 
  ggplot2::annotate(
    "rect", xmin = minRes$lista[5], xmax = minRes$lista[6], 
    ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.2
  ) +
  ggplot2::annotate(
    "rect", xmin = maxRes$lista[5], xmax = maxRes$lista[6], 
    ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.2
  ) +
  ggplot2::geom_vline(
    xintercept = maxRes$lista[[1]], 
    linewidth = 0.3, 
    colour = "darkred") +
  ggplot2::geom_vline(
    xintercept = minRes$lista[[1]], 
    linewidth = 0.3, 
    colour = "darkred") +
  # guides(shape = guide_legend(title = "Proteins with...", override.aes = list(size = 5))) +
  theme(axis.text = element_text(size = sizeLetra),
        axis.title = element_text(size = sizeLetra+2),
        legend.text = element_text(size = sizeLetra-2), 
        legend.title = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.spacing.x = unit(4, 'cm'),
        legend.position = "top"
  ) + 
  xlim(-20, 20) +
  scale_y_continuous("Dataset",  
                     labels = as.character(unique(dfPlot$Dataset)), 
                     breaks = (1:100)) +
  scale_x_continuous(breaks = seq(-20, 20, 5))



#





