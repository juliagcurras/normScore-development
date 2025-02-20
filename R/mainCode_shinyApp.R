###############################################################################-

###################        NORMALIZATION ASSESSMENT   ########################

###############################################################################-

# 2023/09/05 - Julia García Currás
# Somos lo que hacemos cada día. De manera que la excelencia no es un acto, es un hábito. – Aristóteles

setwd("C:/Users/julia/Documents/GitHub/normScore/R")

# 0. Previous work ####

# Load/install libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("NormalyzerDE")
# library(NormalyzerDE)
library(dplyr)
library(Biostatech)
library(ggplot2)
# library(MSnSet.utils)
# library(EnhancedVolcano)

# BiocManager::install("vsn")
# require(vsn);require(preprocessCore);require(matrixStats);
# require(limma);require(stats);require(MASS)

# Load trial data
# archivo <- "Data/trialData.txt"
archivo <- "Data/trialDatashort.xlsx"
# archivo <- "Data/trialDatashort_PossibleErrors.xlsx"
# archivo <- "Data/trialData.xlsx"
# archivo <- "Data/pamelaDataLiver.xlsx"
# archivo <- "Data/raquelData.xlsx"
intensityMatrix <- as.data.frame(readxl::read_excel(path = archivo, sheet = 1, col_names = TRUE))
# intensityMatrix$`Sample Name` <- c(intensityMatrix[1:50,1], intensityMatrix[1:103,1])
# intensityMatrix <- read.table(archivo, header = TRUE, sep = "\t", 
#                               dec = ".",  check.names=FALSE)
# intensityMatrix <- as.data.frame(dfDatos)
intensityMatrix[,-1] <- apply(intensityMatrix[,-1], 2, as.numeric)
rownames(intensityMatrix) <- intensityMatrix$`Sample Name`
intensityMatrix <- intensityMatrix[,-1]
# intensityMatrix <- as.data.frame(apply(intensityMatrix, 2, as.numeric))
logIntensityMatrix <- as.matrix(log(intensityMatrix, base = 2))
logIntensityMatrix[is.infinite(logIntensityMatrix)] <- NA

archivoG <- "./Data/LabelsShort.xlsx"
# archivoG <- "./Data/pamelaGrupos.xlsx"

groupsData <- as.data.frame(readxl::read_excel(path = archivoG, sheet = 1, 
                                               col_names = T))[-3]
colnames(groupsData) <- c("Muestras", "Grupos")
grupos <- levels(as.factor(groupsData$Grupos))


# functions

paleta2 <- function (gradiente = FALSE, show = FALSE, n = NULL, removeWhite = F) {
  paletaDisc <- c("#003C72", "#005B9A", "#1786A3", "#2FB2AD", 
                  "#C3E5BC", "#BCD8E5", "#9389C7", "#D1BCE5", "#C67DD8", 
                  "#CBCBCB")
  gradientePal0 <- grDevices::colorRampPalette(colors = c("#003C72", 
                                                          "#005B9A", "#1786A3", "#2FB2AD", "#BCD8E5", "white", 
                                                          "white", "#D1BCE5", "#C67DD8", "#9389C7", "#544797", 
                                                          "#2E2753"))
  if (removeWhite){
    gradientePal0 <- grDevices::colorRampPalette(colors = c("#003C72", 
                                                            "#005B9A", "#1786A3", "#2FB2AD", "#BCD8E5", "#D1BCE5", "#C67DD8", 
                                                            "#9389C7", "#544797", "#2E2753"))
  }
  gradientePal <- gradientePal0(200)
  if (show) {
    par(mfrow = c(2, 1))
    plot(rep(1, length(paletaDisc)), col = paletaDisc, pch = 15, 
         cex = 10, ylab = "", yaxt = "n", frame.plot = FALSE, 
         xlab = "Paleta discreta", xaxt = "n")
    plot(rep(1, 200), col = gradientePal, pch = 15, cex = 10, 
         ylab = "", yaxt = "n", frame.plot = FALSE, xlab = "Paleta en gradiente", 
         xaxt = "n")
    par(mfrow = c(1, 1))
  }
  else if (!show) {
    if (is.null(n)) {
      if (!gradiente) {
        return(paletaDisc)
      }
      else if (gradiente) {
        return(gradientePal)
      }
    }
    else if (!is.null(n) & is.numeric(n)) {
      if (n == 1) {
        cores <- "#005B9A"
      }
      else if (n == 2) {
        cores <- paletaDisc[c(1, 4)]
      }
      else if (n == 3) {
        cores <- paletaDisc[c(1, 3, 4)]
      }
      else if (n == 4) {
        cores <- paletaDisc[c(1, 3, 4, 7)]
      }
      else if (n == 5 & n < length(paletaDisc)) {
        cores <- paletaDisc[c(1, 3, 4, 7, 9)]
      }
      else if (n <= length(paletaDisc)) {
        cores <- paletaDisc[1:n]
      }
      else if (n > length(paletaDisc)) {
        cores <- gradientePal0(n)
      }
      return(cores)
    }
    else {
      stop("Ha habido algún error.")
    }
  }
}




# 1. Initial assessment ####

#.---- 1.1 Total intensity graph ####
df <- as.data.frame(apply(intensityMatrix, 2, sum, na.rm = T))
colnames(df) <- "Intensity"
df$Samples <- factor(rownames(df), levels = rownames(df), ordered = T)
rownames(df) <- 1:nrow(df)
df$Color <- paleta2(n = nrow(df), removeWhite = TRUE)
rownames(df) <- df$Samples

TIplot <- (ggplot2::ggplot(df, ggplot2::aes(x = Samples, y = Intensity, 
                          text = paste0("Muestra: ", Samples, "\n",
                                        "Intensidad total: ", Intensity))) + 
             ggplot2::geom_bar(stat = "identity", fill = df$Color) + 
             ggplot2::theme_minimal() +
             # ggplot2::ylim(c(0, max(df$Intensity))) +
             ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
                            axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
                            axis.text.x = element_text(angle = 90)) +
             ylab("Intensidad total") + xlab("Muestras")) %>%
  plotly::ggplotly(., tooltip = "text") %>% 
  plotly::config(modeBarButtonsToRemove = c("autoScale2d", "lasso2d", 
                                            "select2d", "pan2d"), 
                                         displaylogo = FALSE)



#.---- 1.2 Boxplot intensity graph ####

IDplot <- Biostatech::plotBoxMultivar(base = as.data.frame(logIntensityMatrix), 
                                       varResumen = colnames(intensityMatrix), 
                                       tituloX = "Muestras", 
                                       tituloY = "Logaritmo de la intensidad total")



#.---- 1.3 Mean intensity graph ###



# 2. Normalization ####

meanNorm <- function(rawMatrix){
  
  # Do calculations
  colMeans <- colMeans(rawMatrix, na.rm = TRUE)
  avgColMean <- mean(colMeans, na.rm = TRUE)
  
  # Create empty matrix
  normMatrix <- matrix(nrow = nrow(rawMatrix), ncol = ncol(rawMatrix), 
                       byrow = TRUE)
  # Normalization 
  normFunc <- function(colIndex) {
    (rawMatrix[rowIndex, colIndex]/colMeans[colIndex]) * 
      avgColMean
  }
  for (rowIndex in seq_len(nrow(rawMatrix))) {
    normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), 
                                     normFunc, 0)
  }
  # Log transformation
  normLog2Matrix <- log2(normMatrix)
  colnames(normLog2Matrix) <- colnames(rawMatrix)
  
  # Output
  return(normLog2Matrix)
}
# meanData <- meanNorm(rawMatrix = intensityMatrix)

medianNorm <- function(rawMatrix) {
  
  rawMatrix <- as.matrix(rawMatrix)
  
  # Do calculations
  colMedians <- matrixStats::colMedians(rawMatrix, na.rm = TRUE)
  meanColMedian <- mean(colMedians, na.rm = TRUE)
  
  # Create empty matrix
  normMatrix <- matrix(nrow = nrow(rawMatrix), ncol = ncol(rawMatrix), 
                       byrow = TRUE)
  
  # Normalization 
  normFunc <- function(colIndex) {
    (rawMatrix[rowIndex, colIndex]/colMedians[colIndex]) * 
      meanColMedian
  }
  for (rowIndex in seq_len(nrow(rawMatrix))) {
    normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), 
                                     normFunc, 0)
  }
  # Log transformation
  normLog2Matrix <- log2(normMatrix)
  colnames(normLog2Matrix) <- colnames(rawMatrix)
  
  # Output
  normLog2Matrix
}
# medianData <- medianNorm(rawMatrix = intensityMatrix)


GINorm <- function(rawMatrix){
  
  # Do calculations
  colSums <- colSums(rawMatrix, na.rm = TRUE)
  colSumsMedian <- stats::median(colSums)
  
  # Create empty matrix
  normMatrix <- matrix(nrow = nrow(rawMatrix), ncol = ncol(rawMatrix), 
                       byrow = TRUE)
  
  # Normalization 
  normFunc <- function(colIndex) {
    (rawMatrix[rowIndex, colIndex]/colSums[colIndex]) * 
      colSumsMedian
  }
  for (rowIndex in seq_len(nrow(rawMatrix))) {
    normMatrix[rowIndex, ] <- vapply(seq_len(ncol(rawMatrix)), 
                                     normFunc, 0)
  }
  
  # Log transformation
  normLog2Matrix <- log2(normMatrix)
  colnames(normLog2Matrix) <- colnames(rawMatrix)
  
  # Output
  normLog2Matrix
}
# giData <- GINorm(rawMatrix = intensityMatrix)

quantileNorm <- function(log2Matrix){
  log2Matrix <- as.matrix(log2Matrix)
  normMatrix <- preprocessCore::normalize.quantiles(log2Matrix, 
                                                    copy = TRUE)
  colnames(normMatrix) <- colnames(log2Matrix)
  normMatrix
}
# quantileData <- quantileNorm(log2Matrix = logIntensityMatrix)


cyclicLoessNorm <- function (log2Matrix) {
  log2Matrix <- as.matrix(log2Matrix)
  normMatrix <- limma::normalizeCyclicLoess(log2Matrix, method = "fast")
  colnames(normMatrix) <- colnames(log2Matrix)
  normMatrix
}
# cycLoessData <- cyclicLoessNorm(log2Matrix = logIntensityMatrix)

VSNNorm <- function (rawMatrix) {
  rawMatrix <- as.matrix(rawMatrix)
  
  normMatrix <- suppressMessages(vsn::justvsn(rawMatrix))
  colnames(normMatrix) <- colnames(rawMatrix)
  normMatrix
}
# vsnData <- VSNNorm(rawMatrix = intensityMatrix)


RLRNorm <- function(log2Matrix){
  
  log2Matrix <- as.matrix(log2Matrix)
  
  log2Matrix[is.infinite(log2Matrix)] <- NA
  
  sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm = TRUE)
  
  calculateRLMForCol <- function(colIndex, sampleLog2Median, 
                                 log2Matrix) {
    lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex]) ~ 
                         sampleLog2Median, na.action = stats::na.exclude)
    coeffs <- lrFit$coefficients
    coefIntercept <- coeffs[1]
    coefSlope <- coeffs[2]
    globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept)/coefSlope
    globalFittedRLRCol
  }
  globalFittedRLR <- vapply(seq_len(ncol(log2Matrix)), calculateRLMForCol, 
                            rep(0, nrow(log2Matrix)), 
                            sampleLog2Median = sampleLog2Median, 
                            log2Matrix = log2Matrix)
  colnames(globalFittedRLR) <- colnames(log2Matrix)
  globalFittedRLR
}
# rlrData <- RLRNorm(log2Matrix = logIntensityMatrix)


MADNormalization <- function(log2Matrix){
  
  log2Matrix <- as.matrix(log2Matrix)
  
  sampleLog2Median <- matrixStats::colMedians(log2Matrix, 
                                              na.rm = TRUE)
  sampleMAD <- matrixStats::colMads(log2Matrix, na.rm = TRUE)
  madMatrix <- t(apply(log2Matrix, 1, function(row) ((row - 
                                                        sampleLog2Median)/sampleMAD)))
  madPlusMedianMatrix <- madMatrix + mean(sampleLog2Median)
  colnames(madPlusMedianMatrix) <- colnames(log2Matrix)
  madPlusMedianMatrix
}
# madData <-  MADNormalization(log2Matrix = logIntensityMatrix)

# Outros: Trimmed Mean of M Values (TMM), EigenMS (EIG), Locally Weighted Scatterplot Smoothing (LOW),



# mydata <- list(Mean = meanData, Median = medianData, GI = giData, 
#               Quantile = quantileData, CyclicLoess = cycLoessData, 
#               RLR = rlrData, MAD = madData)

mydata <- list(log = logIntensityMatrix,
               Mean = meanNorm(rawMatrix = intensityMatrix), 
               Median = medianNorm(rawMatrix = intensityMatrix), 
               GI = GINorm(rawMatrix = intensityMatrix),
               Quantile = quantileNorm(log2Matrix = logIntensityMatrix),
               VSN = VSNNorm(rawMatrix = intensityMatrix),
               CyclicLoess = cyclicLoessNorm(log2Matrix = logIntensityMatrix),
               RLR = RLRNorm(log2Matrix = logIntensityMatrix),
               MAD = MADNormalization(log2Matrix = logIntensityMatrix))




# 3. Assessment ####

#.---- 3.1 visual aids ####

#.--------- 3.1.1 Distribution graphs ####

### ..Boxplot ####

IDNormPlot <- lapply(1:length(mydata), function(i) {
  bbdd <- as.data.frame(mydata[[i]])
  plotBoxMultivar1(base = bbdd, 
                   varResumen = colnames(bbdd), 
                   tituloX = "Muestras", 
                   tituloY = paste0("Intensidad total normalizada por ", names(mydata)[i]))
})
names(IDNormPlot) <- names(IDNormPlot)
names(IDNormPlot) <- names(IDNormPlot)

### RLE plot

plotRLE <- function(dfDatos, normalizacion) {
  medianaProt <- apply(dfDatos, 1, stats::median)
  
  rleData <- as.data.frame(log(t(t(dfDatos) / medianaProt), base = 2))
  
  plotBoxMultivar1(base = rleData, 
                   varResumen = colnames(rleData), 
                   titulo = normalizacion,
                   tituloX = "Muestras", 
                   tituloY = "RLE")
  
}
rlePlots <- lapply(names(mydata), function(i) plotRLE(dfDatos = mydata[[i]], normalizacion = i))
names(rlePlots) <- names(rlePlots)


### meanSDplot

meanSDplot <- function(dfDatos, normalizacion, color){
  
  meanSamples <- apply(dfDatos, 2, mean, na.rm = T)
  sdSamples <- apply(dfDatos, 2, sd, na.rm = T)
  
  dfAux <- data.frame(Media = meanSamples,
                      DesvEst = sdSamples, 
                      Muestras = colnames(dfDatos))
  
  dfAux <- dfAux %>% arrange(Media)
  dfAux$Orden <- 1:nrow(dfAux)
  plotScatter(base = dfAux, varX = "Orden", varY = "DesvEst", 
              corType = "Pearson", titulo = normalizacion, 
              tituloX = "Muestras ordenadas según media de intensidad", 
              tituloY = "Desviación estándar (SD)", color = color)
  
}
cores <- Biostatech::colorPalette(n=length(mydata))
names(cores) <- names(mydata)
meanSDPlots <- lapply(names(mydata), function(i) meanSDplot(dfDatos = mydata[[i]], 
                                                            normalizacion = i, 
                                                            color = as.vector(cores[i])))
names(meanSDPlots) <- names(mydata)
meanSDPlots[[1]]
meanSDPlots[[9]]

#.--------- 3.1.2 Metrics ####

# Necesitamos grupos homogéneos... 


### Coefficient of Variation
# Pódese adaptar fácilmente a outras métricas como PEVC cambiando solo a función
# de cálculo
coefVariation <- function(x, na.rm = TRUE) {
  (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))*100
} 

cvGruposProt <- function(grupo, dfGrupos, dfDatos){
  df <- dfDatos[, colnames(dfDatos) %in% dfGrupos[dfGrupos$Grupos == grupo, "Muestras"]]
  apply(df, 1, coefVariation)
}

getPCV <- function(dfDatos, grupos, dfGrupos){
  
  # Calculate cv by groups for each protein
  dfGruposProt <- as.data.frame(sapply(grupos, cvGruposProt, 
                                       dfGrupos = dfGrupos, 
                                       dfDatos = dfDatos))
  # Calculate mean and mean IC of CV for each protein
  df1 <- as.vector(as.matrix(
    sapply(apply(dfGruposProt, 2, getMeanIC1), "[[", 1) # LI e LS
  ))
  names(df1) <- as.vector(sapply(grupos, paste0, c(": Media CV", ": LI", ": LS"), simplify = T))
  
  return(df1)
}

getPCV(dfDatos = mydata$Mean, grupos = grupos, dfGrupos = groupsData)
dfPCV <- data.frame(lapply(mydata, getPCV, grupo = grupos, dfGrupos = groupsData))
# rownames(dfPCV) <- grupos

plotForestGroup(etiquetas = rep(colnames(dfPCV), 2), 
                estPunt = as.vector(t(as.matrix(dfPCV[seq(1, nrow(dfPCV), 3),]))), 
                LI = as.vector(t(as.matrix(dfPCV[seq(2, nrow(dfPCV), 3),]))), 
                LS = as.vector(t(as.matrix(dfPCV[seq(3, nrow(dfPCV), 3),]))), 
                grupos = rep(grupos, each = ncol(dfPCV)), vertical = F,
                tituloX = "Media del coeficiente de variación por grupos - PVC (%)",
                referenceLine = F
                )



plotForestGroup <- function(etiquetas, estPunt, LI, LS, grupos = NULL, tipoEstPunt = "EP", titulo = "", tituloY = "",
                       tituloX = "Estimación puntual", vertical = T, referenceLine = T, color = NULL, interact = TRUE) {
  library(magrittr)
  
  if (tituloX == "Estimación puntual") {
    if (tipoEstPunt == "HR") {
      tituloX <- "Hazard Ratio (HR)"
    } else if (tipoEstPunt == "OR") {
      tituloX <- "Odds Ratio (OR)"
    }
  }
  
  # Datos
  if (is.data.frame(etiquetas)) {
    if (!all(colnames(etiquetas) == c("Variable", "Categoría", "Referencia"))) {
      stop("Input dataframe is not the required one from script_tablas.R")
    }
    # varPrev <- "" # se puede implementar para que no se repita el nombre de variable cuando son más de 2 categorías
    etiquetas <- apply(etiquetas, 1, function(fila) {
      if (fila[2] == "") {
        paste0("<b>", fila[1], "</b>")
      } else { # if (varPrev != fila[1]) {
        # varPrev <<- fila[1]
        paste0("<b>", fila[1], "</b>:\n <i>", fila[2], "</i> vs <i>", fila[3], "</i>")
        # } else if (varPrev == fila[1]){
        #   paste0("<i>", fila[2], "</i> vs <i>", fila[3], "</i>")
      }
    })
  }
  # remove(varPrev)
  
  if (!is.null(grupos) & length(grupos) != length(etiquetas)){
    stop("La cantidad de grupos introducidos no coincide con la cantidad de etiquetas!!")
  }
  
  datos <- data.frame(
    Etiquetas = etiquetas,
    OR = as.numeric(estPunt),
    li = as.numeric(LI),
    ls = as.numeric(LS),
    Grupos = (if (is.null(grupos)) {1} else {grupos})
  )
  datos$Grupos <- factor(datos$Grupos)
  nGrupos <- length(levels(datos$Grupos))
  
  pd <- ggplot2::position_dodge(0.6) # move them .2 to the left and right
  
  cores <- c("#1786A3", "#005B9A")
  if (is.null(grupos)){
    if (!is.null(color) & length(color) == 2) {
      cores <- color
    } else if (!is.null(color) & length(color) != 2) {
      stop("Número de colores introducido no válido. Se establecerán los colores por defecto.")
    }
  } else {
    cores <- Biostatech::colorPalette(n = nGrupos)
    if (!is.null(color) & length(color) == nGrupos) {
      cores <- color
    } else if (!is.null(color) & length(color) != nGrupos) {
      warning("Número de colores introducido no válido. Se establecerán los colores por defecto.")
    }
  }

  
  grafico <- ggplot2::ggplot(
    data = datos,
    ggplot2::aes(
      x = factor(Etiquetas, levels = rev(unique(Etiquetas))),
      y = OR,
      shape = Grupos, 
      color = Grupos, 
      text = paste0(
        "<b>", Etiquetas, "</b>",
        "<br>___________________",
        "<br>", tipoEstPunt, ": ", format(OR, nsmall = 3L),
        "<br>Límite inferior: ", format(li, nsmall = 3L),
        "<br>Límite superior: ", format(ls, nsmall = 3L)
      )
    )
  )
  
  
  grafico <- grafico +
    ggplot2::geom_point(position = pd, size = 1.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = li, ymax = ls),
                           width = 0.05, linewidth = 0.5,
                           position = pd) +
    ggplot2::ylim(0, max(datos$ls)) +
    ggplot2::ggtitle(titulo) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = tituloY, y = tituloX) +
    ggplot2::theme_minimal(base_family = "Calibri") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    ) 
  
  if (referenceLine){
    grafico <- grafico +ggplot2::geom_hline(
      yintercept = ifelse(tipoEstPunt %in% c("OR", "HR"), 1, 0),
      linewidth = 0.3, colour = "black"
    )
    
  }
  
  if (vertical){
    grafico <- grafico + ggplot2::coord_flip()
  }
  
  if (is.null(grupos)){
    grafico <- grafico +
      ggplot2::geom_point(position = pd, colour = cores[2], size = 1.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = li, ymax = ls),
                             width = 0.05, linewidth = 0.5,
                             position = pd, colour = cores[1]
      ) +
      ggplot2::theme(legend.position = "none") 
  } else {
    grafico <- grafico + scale_color_manual(values = cores)
  }
  
  
  # Interactive
  if (interact) {
    m <- list(
      l = 0,
      r = 150,
      b = 10,
      t = 50,
      pad = 4
    )
    grafico <- plotly::ggplotly(grafico, tooltip = "text")
    grafico <- grafico %>%
      # plotly::layout(legend = list(orientation = "h", x=0.1 ,y =-0.2)) %>%
      plotly::config(
        modeBarButtonsToRemove = c(
          "autoScale2d", "lasso2d", "select2d",
          "pan2d", "hoverCompareCartesian", "hoverClosestCartesian"
        ),
        displaylogo = FALSE
      ) %>%
      plotly::layout(autosize = T, margin = m)
  }
  
  return(list(grafico = grafico))
}

# etiquetas <- LETTERS[1:6]
# etiquetas <- data.frame(Variable = rep(c("Variable 1", "Variable 2", "Variable 2"), 2),
#                            Categoría = rep(c("Cat1", "Cat2.1", "Cat2.2"), 2),
#                            Referencia = rep(c("Ref1", "Ref2", "Ref2"),2))
# LI <- list(1.2, 1.1, 1.6, 1.4, 1.5, 1.5)
# estPunt <- list(1.5, 1.2, 1.8, 1.9, 1.75, 1.95)
# LS <- list(1.7, 1.3, 2.0, 2.3, 1.9, 2.6)
# grupos <- rep(c("Control", "Case"), each = 3)
# grupos <- NULL




PCVplot <- plotBoxMultivar1(base = dfPCV, 
                 varResumen = colnames(dfPCV), 
                 tituloX = "Normalizaciones", 
                 tituloY = "PVC (%)")
PCVplot





### Coefficient of correlationn (Pearson)

getCorrelationVector <- function(df, dfGrupos, metodo = "pearson"){
  df <- as.data.frame(df)
  allCorrs <- lapply(unique(dfGrupos$Grupos), function(i){
    # Select samples
    samplesByGroup <- dfGrupos %>% 
      dplyr::filter(Groups == i) %>%
      dplyr::pull(Samples)
    dfCor <- df %>% 
      dplyr::select(all_of(samplesByGroup))
    # Estimate correlation
    tabCor <- stats::cor(dfCor, use = "complete.obs", method = metodo)
    # Get unique pairs of correlation
    tabCor[lower.tri(tabCor)] <- NA # remove duplicated corr
    diag(tabCor) <- NA # remove variance diagonal
    tabCor <- stats::na.omit(reshape2::melt(tabCor))
  }
  )
  final <- as.data.frame(dplyr::bind_rows(allCorrs)) %>% dplyr::pull(value)
  return(final)
}

allVectorsCorr <- sapply(mydata, getCorrelationVector, 
              dfGrupos = groupsData, 
              metodo = "pearson")

dfPlot <- data.frame(sapply(allVectorsCorr, "length<-", max(lengths(allVectorsCorr))))

plotBoxMultivar(dfPlot, varResumen = colnames(dfPlot))



# 4. Comparisons ####

#.--- 4.1 Test t Student ####

doTestT <- function(df, dfGrupos, g1, g2) {
  
  # grupo 1 - control
  df <- as.data.frame(df)
  df[sapply(df, is.infinite)] <- NA
  
  samplesG1 <- dfGrupos[dfGrupos$Grupos == g1, "Muestras"]
  samplesG2 <- dfGrupos[dfGrupos$Grupos == g2, "Muestras"]
  df <- df[, c(samplesG1, samplesG2)]
  
  # log2 FC change group2 - group1 
  df$logFC  <- apply(df, 1, function(x) {
      mean(x[samplesG2], na.rm =T) - mean(x[samplesG1], na.rm = T)
    }
  )
  
  df$AveExpr  <- apply(df, 1, function(x) {
    (mean(x[samplesG2], na.rm =T) + mean(x[samplesG1], na.rm = T))/2
  }
  )
  
  #T-test with equal variance
  df[, c( "t", "P.Val")] <- t(apply(df, 1, function(x) {
      res <- t.test(x[samplesG2], x[samplesG1],
                    alternative = "two.sided", 
                    var.equal = TRUE)
      return(c(res$statistic, res$p.value))
    }
  ))
  
  #Benjamini-Hochberg correction for multiple testing
  df$adj.P.Val <- p.adjust(df$P.Val, method = "BH")
  df <- df[, c("logFC", "AveExpr", "t", "P.Val", "adj.P.Val")]
  return(df)
}

df <- as.data.frame(mydata[["Mean"]])
df[sapply(df, is.infinite)] <- NA
df$genes <- substr(rownames(intensityMatrix), start = 4, stop = 16)
rownames(df) <- df$genes
df1 <- doTestT(df = df, dfGrupos = groupsData, g1 = "Control", g2 = "Case")
resultTContrast <- lapply(mydata, doTestT, dfGrupos = groupsData, g1 = "Control", g2 = "Case")

df1$genes <- substr(rownames(intensityMatrix), start = 4, stop = 10)
df1 <- resultTContrast[[6]]

# plot volcano
EnhancedVolcano(df1,
                lab = rownames(df1),
                x = 'logFC',
                y = 'adj.P.Val',
                pointSize = 3.0,
                labSize = 3)


#.--- 4.2 Limma contrast ####
# Do limma contrast

doLimmaContrast <- function(df, dfGrupos , g1 = "Control", g2 = "Case"){
  
  library(limma)
  df <- as.data.frame(df)
  df[sapply(df, is.infinite)] <- NA
  
  #Define the design vector
  idRows <- sapply(colnames(df), function(i){
    which(dfGrupos$Muestras == i)
  })
  cond <- as.factor(dfGrupos[idRows, "Grupos"])
  design <- model.matrix(~0+cond)
  colnames(design) = gsub("cond", "", colnames(design))
  
  #Make contrasts
  contrast =  makeContrasts(contrasts=paste0(g2,"-", g1), levels=design)
  fit1 <- lmFit(df, design)
  fit2 <- contrasts.fit(fit1, contrasts = contrast)
  fit3 <- eBayes(fit2)
  
  # get results
  limmaTable <- limma::topTable(fit3, coef=1, number=Inf, sort.by="none")
  
  # Output
  return(limmaTable)
}


df <- as.data.frame(mydata[["Mean"]])
rownames(df) <- substr(rownames(intensityMatrix), start = 4, stop = 15)
df$genes <- NULL
# resultLimma <- doLimmaContrast(df = df, dfGrupo = groupsData, g1 = "Control", g2 = "Case")
resultLimma <- doLimmaContrast(df = df, dfGrupo = groupsData, g1 = "WT_Liver", g2 = "KO_Liver")
resultTContrast <- lapply(mydata, doLimmaContrast, dfGrupos = groupsData, g1 = "Control", g2 = "Case")

EnhancedVolcano(resultLimma,
                lab = rownames(resultLimma),
                x = 'logFC',
                y = 'adj.P.Val',
                pointSize = 3.0,
                labSize = 2)


# Volcano plot - handmade ####

# Get y axis values
rownames(df1) <- substr(rownames(intensityMatrix), start = 4, stop = 15)

data <- df1
maxAdjP <- 0.05
minFC <- 1
logFC = 'logFC'
adjPVal = 'adj.P.Val'
proteinas = NULL
g1='Control'
g2 = 'Case'
maxAdjP = 0.05
minFC = 1
testType = "T-test"

plotVolcano <- function(data, logFC = 'logFC', adjPVal = 'adj.P.Val',
                        proteinLabel = T, proteinas = NULL, 
                        g1='Control', g2 = 'Case', testType = "T-test", 
                        maxAdjP = 0.05, minFC = 1 #, interact = T
                        ){ 

  if (!any(logFC %in% colnames(data))){
    stop("El nombre de la variable que contiene el logFC no se encuentra en la base introducida.")
  } else if (!any(adjPVal %in% colnames(data))){
    stop("El nombre de la variable que contiene el valor p ajustado no se encuentra en la base introducida.")
  } 

  df <- data.frame(logFC = as.numeric(data[[logFC]]), 
                   adjPVal = as.numeric(data[[adjPVal]]),
                   Log10adjPval = (-1*log10(data[[adjPVal]]))
                   )
  
  if (proteinLabel){
    if (is.null(proteinas)){
      df$Proteinas <- rownames(data)
    } else if (is.character(proteinas)){
      if (!any(proteinas %in% colnames(data))){
        stop("El nombre de la variable que contiene las proteínas no se encuentra en la base introducida.")
      } else {
        df$Proteinas <- data[[proteinas]]
      }
    } else if (is.vector(proteinas) & length(proteinas) == nrow(df)){
      df$Proteinas <- proteinas
    } else {
      stop("Error con el argumento del parémetro proteinas.")
    }
  }
  

  
  # Clasificar según significación y logFC
  df$DiffAbund <- apply(
    df[, c("adjPVal", "logFC")], 1, function(x) {
      if (x[1] <= maxAdjP & x[2] >= minFC) {
        return( paste("Cantidad mayor en", g2))
      } else if (x[1] <= maxAdjP & x[2] <= (-1*minFC)) {
        return(paste("Cantidad mayor en", g1) )
      } else {
        return('No significativas')
      }
    }
  )
  
  df$DiffAbund <- factor(df$DiffAbund, 
                         levels = c("No significativas", paste("Cantidad mayor en", g1),
                                    paste("Cantidad mayor en", g2)))
  
  maxFC <- max(c(abs(min(df$logFC)), max(df$logFC)))*1.5
  maxPval <- max(df$Log10adjPval)*1.25
  
  
  grafico <- ggplot(
    df,
    aes(x = logFC, y = Log10adjPval, color = DiffAbund)
    ) +
    ggplot2::geom_point(shape = 19, size=3, alpha = 0.7)+
    ggplot2::geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65", 
                        show.legend = F) +
    # ggplot2::geom_vline(xintercept = 0, colour = "gray65") +
    ggplot2::geom_vline(xintercept = -1*minFC, colour = "gray65", 
                        show.legend = F) +
    ggplot2::geom_vline(xintercept = minFC, colour = "gray65", 
                        show.legend = F) +
    ggplot2::ggtitle(paste0("Comparación ", g2, " vs ", g1)) +
    xlim((0-maxFC),(0+maxFC)) +
    ylim(c(0, maxPval)) +
    scale_color_manual(values = c("#005B9A", "#DC143C", "#32CD32")) +
    ggplot2::theme_minimal() + # base_family = "Calibri"
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
      axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      legend.title = ggplot2::element_blank(), 
      legend.position = 'top',
      legend.text = ggplot2::element_text(size=14),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      plot.caption = element_text(size = 10, face = "italic"),
      plot.title = element_text(size = 18, hjust = 0.5)
      ) +
    # ggplot2::guides(color = guide_legend(override.aes = list(size = 5)))+
    ggplot2::labs(x = paste("Log2 FC", g2, "-", g1),
                  y = "-Log10 Adj. P-value",
                  caption = paste0("Contraste: ", testType, 
                                   "; Valores p ajustados <=", maxAdjP,
                                   "; -", minFC, "<= Log2 FC >= ", minFC)) 
  
  if (proteinLabel){
    grafico  <- grafico +
      ggplot2::geom_text(data = subset(df, logFC >= minFC | logFC <= -1*minFC),
                         aes(logFC, Log10adjPval, label = Proteinas),
                         alpha = 0.6, hjust = 0.5, vjust = -0.6, show.legend = F
      )
  }
  
  
  # if (interact){ # NO ENSEÑA LA INFORMACIÓN SOBRE LOS PUNTOS PORQUE LOS TOMA COMO CATEGORÍAS
  #   grafico <- suppressWarnings(plotly::ggplotly(grafico,
  #                                                tooltip = "text",
  #                                                dynamicTicks = TRUE,
  #                                                layerData = 1))
  #   grafico <- grafico %>% plotly::layout(showlegend = TRUE) %>%
  #     plotly::config(modeBarButtonsToRemove = c("autoScale2d", "lasso2d",
  #                                               "select2d", "pan2d"),
  #                    displaylogo = FALSE)
  # }
  
  return(list(grafico = grafico))
}

plotVolcano(data = df1, proteinLabel = T)



# 5. Assessment ####
df # datos mean
df1 # resultados contraste
# df[sapply(df, is.infinite)] <- NA



dfResults <- xlsx::read.xlsx(file = "Data/comparisonResults/comparisonResults.xlsx", sheetIndex = 1)
protSign <- dfResults[-1,1]

dfGrupos <- groupsData
library(factoextra)

df <- df[,-147] # req(dfQuant()) # non queremos a info dos genes, esta nos rownames

g1 <- "Control"
g2 <- "Case"
if ((!(g1 %in% dfGrupos$Grupos)) | (!(g2 %in% dfGrupos$Grupos))){
  return(NULL)
}

samplesG1 <- dfGrupos[dfGrupos$Grupos == g1, "Muestras"]
samplesG2 <- dfGrupos[dfGrupos$Grupos == g2, "Muestras"]
dfQuant <- df[, c(samplesG1, samplesG2)]

dfResults <- df1 # req(resultComp())
maxAdjP = 0.05
minFC = 1
dfGrupos<- groupsData

getPCA <- function(dfQuant, dfResults, dfGrupos, maxAdjP = 0.05, minFC = 1){
  
  # Get significant prots
  protSign <- dfResults %>% 
    filter(adj.P.Val < maxAdjP) %>% 
    filter(logFC <= -1*minFC | logFC >= minFC) %>% rownames
  
  # Select quantitative data from significative prots
  dfQuant$genes <- rownames(dfQuant)
  
  dfAux <- dfQuant %>% 
    dplyr::filter(genes %in% protSign) %>% 
    dplyr::select(-genes) 
  dfAux <- as.data.frame(t(dfAux))
  dfAux$IDs <- rownames(dfAux)
  
  # Identify groups for each patient
  df1 <- merge(dfAux, dfGrupos, by.x = "IDs", by.y = "Muestras") 
  df1 <- df1[,c(2:(ncol(df1)-1), 1, ncol(df1))]
  df1$Grupos <- factor(df1$Grupos)
  
  # Identify groups for each patient
  # dfGrupos <- as.data.frame(dfGrupos)
  # idRows <- sapply(rownames(dfAux), function(i){
  #   which(dfGrupos$Muestras == i)
  # }, simplify = T, USE.NAMES = F)
  # dfAux$Grupos <- as.factor(dfGrupos[idRows, "Grupos"])
  
  # Change infinite data for NA
  df1[sapply(df1, is.infinite)] <- NA
  
  if (sum(is.na(df1)) != 0){
    recuentoNARow <- apply(df1, 2, function(j) sum(is.na(j)))
    if (any(recuentoNARow == ncol(df1))){ # remove empty cols
      df1[, names(recuentoNARow[recuentoNARow == ncol(df1)])] <- NULL
    }
    df2 <- VIM::kNN(df1[,-c(ncol(df1)-1, ncol(df1))])[1:(ncol(df1)-2)] # imputación por k nearest neighbours
  } else { # Imputation was not necessary
    df2 <- df1[,-c(ncol(df1)-1, ncol(df1))]
  }
  
  # Imputation (if necessary)
  # if (sum(is.na(dfAux)) != 0){
  #   recuentoNARow <- apply(dfAux, 2, function(j) sum(is.na(j)))
  #   if (any(recuentoNARow == ncol(dfAux))){ # remove empty cols
  #     dfAux[, names(recuentoNARow[recuentoNARow == ncol(dfAux)])] <- NULL
  #   }
  #   df2 <- VIM::kNN(dfAux[,-c(ncol(dfAux)-1, ncol(dfAux))])[1:(ncol(dfAux)-2)] # imputación por k nearest neighbours
  # } else { # Imputation was not necessary
  #   df2 <- dfAux
  # }
  
  #.--- 5.1 PCA ####
  res.pca <- prcomp(df2, scale = FALSE)
  # PCA main figure
  pcaBiplot <- fviz_pca_biplot(res.pca, repel = TRUE, habillage = df1$Grupos,
                               # col.ind = dfAux$Grupos, palette = "jco",
                               palette = colorPalette(n = length(levels(df1$Grupos))),
                               addEllipses = TRUE, label = "var",
                               col.var = "black", title = "Biplot") + 
    theme_minimal() +
    theme(text = element_text(size = 16),
          title = element_text(size = 20),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))
  # fig1 <- fviz_eig(res.pca, addlabels=TRUE) + #barfill=cores[4], barcolor = cores[4], linecolor = cores[1])
  #   theme_minimal()
  # correlation components
  varCor <- get_pca_var(res.pca)
  # pcaCorComponets <- corrplot::corrplot(varCor$cos2, is.corr=FALSE)
  
  
  #.--- 5.2 Clustering ####
  df3 <- df1[,-c(ncol(df1)-1, ncol(df1))]
  rownames(df3) <- df1$IDs
  d <- dist(df3)
  hc1 <- hclust(d, method = "complete")
  # cophecor <- cor(d, cophenetic(hc1))
  
  # Expected and observed groups
  df1$NewGroups <- factor(cutree(hc1, k = 2), levels = 1:2, 
                          labels = c("Grupo 1", "Grupo 2"))
  
  dfClases <- as.data.frame(summary(arsenal::tableby(NewGroups~Grupos, df1, test = F, 
                           cat.stats = c("countpct", "N")), text = T))
  # Tree representation
  library(dendextend)
  df1$Cores <- factor(df1$Grupos, labels = colorPalette(n = 2))
  hc2 <- as.dendrogram(hc1)
  labels(hc2) <- df1[order.dendrogram(hc2), "IDs"]
  labels_colors(hc2) <- as.vector(df1[order.dendrogram(hc2), "Cores"])
  
  return(list(figura = pcaBiplot, infoVar = varCor, gruposC = dfClases, figuraC = hc2))
}

abc <- getPCA(dfQuant = df, dfResults = df1, dfGrupos = dfGrupos)

abc$figura
corrplot::corrplot(abc$infoVar$cos2, is.corr=FALSE)


#.--- 5.2 Clustering ####
# Clustering analysis 
df1 <- dfAux # req(dfQuant())
dfResults <- df1 # req(resultComp())
maxAdjP = 0.05
minFC = 1


# Clustering analysis 
df3 <- df1[,df1[,-c(ncol(df1)-1, ncol(df1))]]
d <- dist(df3)
hc1 <- hclust(d, method = "complete")
cophecor <- cor(d, cophenetic(hc1))
df1$NewGroups <- factor(cutree(hc1, k = 2), levels = 1:2, 
                        labels = c("Grupo 1", "Grupo 2"))

summary(arsenal::tableby(NewGroups~Grupos, df1, test = F, 
                         cat.stats = c("countpct", "N")), text = T)


library(dendextend)
df1$Cores <- factor(df1$Grupos, labels = colorPalette(n = 2))

hc2 <- as.dendrogram(hc1)
labels(hc2) <- df1[order.dendrogram(hc2), "IDs"]
labels_colors(hc2) <- as.vector(df1[order.dendrogram(hc2), "Cores"])
# set(hc2, "labels_cex", 0.8)
plot(hc2)
rect.dendrogram(hc2, k = 2, border = colorPalette()[9])












# ANEXO 1 | Errores ####

#.- Test does not work changing gruop names ####

setwd("C:/Users/julia/Documents/GitHub/Normalization")

intensityMatrix <- readxl::read_xlsx(path = "Data/moreData/AEGIS.xlsx", col_names = T)
intensityMatrix[,-1] <- apply(intensityMatrix[,-1], 2, as.numeric) 
intensityMatrix[,1] <- changeNames(protNames = intensityMatrix[,1])
dfRaw <- as.data.frame(intensityMatrix[,-1])

logIntensityMatrix <- as.matrix(log(dfRaw, base = 2))
logIntensityMatrix[is.infinite(logIntensityMatrix)] <- NA

normData <- doNormalization(listaNorm = c("Media", "Mediana", "MAD"), rawData = dfRaw,
                            logData = logIntensityMatrix)

normData$Log <- logIntensityMatrix
normData <- normData[c(length(normData), 1:(length(normData)-1))]

dfGrupos <- readxl::read_xlsx(path = "Data/moreData/AEGISgrupos.xlsx", col_names = T)
colnames(dfGrupos) <- c("Muestras", "Grupos")

df <- normData$MAD

resultComp <- doTestT(df = df, dfGrupos = as.data.frame(dfGrupos),
                      g1 = "G1",
                      g2 = "G9")
t(apply(df[121,], 1, function(x) { # NA error
  
  evalG1NA <- sum(is.na(x[samplesG1])) # When there's enought observations
  evalG2NA <- sum(is.na(x[samplesG2]))
  
  if (((length(x[samplesG1])-evalG1NA) < 3 ) | 
      ((length(x[samplesG2])-evalG2NA) < 3 )) {
    return(c(NA, NA))
  }
  
  res <- stats::t.test(x[samplesG2], x[samplesG1],
                       alternative = "two.sided",
                       var.equal = TRUE)
  return(c(res$statistic, res$p.value))
}))


dfGrupos$Grupos <- sapply(dfGrupos$Grupos, function(i){
  if (!is.na(suppressWarnings(as.numeric(i)))){
    return(paste0("G", i))
  } else {
    i
  }
})





# ANEXO 2 | CÓDIGO INTERESANTE ####

#.- Volcano plot con ggplot ####

# Get y axis values
# df$Log10adjPval <- -1*log10(df$adjPval)
#
# Add the categorical column for easier visualization
# df$Diff_Abund <- apply(
#   df, 1, function(x) {
#     if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] >= minFC) {
#       return( paste("Up in", gr2) )
#     } else if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] <= -1*minFC) {
#       return( paste("Up in", gr1) )
#     } else {
#       return('Non-significant')
#     }
#   }
# )
# 
# 
# ggplot(
#   df,
#   aes(x = Log2FC, y = Log10adjPval, colour = Diff_Abund )
# ) +
#   geom_point(shape=19, size=2, alpha = 0.6)+
#   geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65") +
#   geom_vline(xintercept = 0, colour = "gray65") +
#   geom_vline(xintercept = -1*minFC, colour = "gray65") +
#   geom_vline(xintercept = minFC, colour = "gray65") +
#   ggtitle(
#     paste(
#       "T-test ", g1, " vs ", g2,
#       " Adjusted P-value<=", maxAdjP, " Log2 FC>=", minFC,
#       sep=""
#     )
#   ) +
#   theme_classic() +
#   theme(
#     legend.title = element_blank(), legend.text = element_text(size=12),
#     plot.title = element_text(size=16)
#   ) +
#   labs(x = paste("Log2 FC", g2, "-", g1), y = "-Log10 Adj. P-value" ) +
#   geom_text(
#     data = subset(df, Log2FC >=0.9 | Log2FC <= -0.8),
#     aes(Log2FC, Log10adjPval, label = genes),
#     alpha = 0.6, hjust = 0.5, vjust = -0.6
#   )




#.- Contrastes con MSnSet.utils ####

# Old code 
# swath <- as.data.frame(mydata[["Mean"]])
# swath$SampleName <- rownames(intensityMatrix)
# swath <- swath[,c(ncol(swath), 1:(ncol(swath)-1))] # ID ten que ir de primeira col necesariamente
# dfGrupos <- groupsData
# 
# write.csv(x=as.matrix(swath), file="Data/swath_trial_2.csv", row.names = F)
# swathPath <- "Data/swath_trial_2.csv"
# 
# msn <- readMSnSet2(swathPath, ecol=2:ncol(swath), fnames=1, header=T)
# msn$group <- as.vector(sapply(colnames(swath)[-1], function(i) { # simplemente para asegurarse del orgen de las muestras, por si es diferente entre dfGrupos$muestra y colnames(swath)
#   dfGrupos[dfGrupos$Muestras == i, "Grupos"]
# }))
# msn$group = factor(msn$group, levels = c("Control", "Case"))
# 
# 
# t_res1 <- limma_a_b(eset = msn, model.str = "~ group", 
#                     coef.str = "group")
# t_res1$Protein = rownames(t_res1)
# 
# 
# proteinas <- substr(swath$SampleName, start = 11, stop = 15)
# 
# EnhancedVolcano(t_res1,
#                 lab = proteinas,
#                 x = 'logFC',
#                 y = 'P.Value',
#                 pointSize = 3.0,
#                 labSize = 6)


# Adaptado para la shiny app
#
# doContrast <- function(data, dfGrupos, controlG = "Control", 
#                        caseG = "Case"){
#   
#   write.csv(x=as.matrix(data), file="swath_trial.csv", row.names = F)
#   swathPath <- "swath_trial.csv"
#   
#   msn <- readMSnSet2(swathPath, ecol=2:ncol(data), fnames=1, header=T)
#   msn$group <- as.vector(sapply(colnames(data)[-1], function(i) { # simplemente para asegurarse del orgen de las muestras, por si es diferente entre dfGrupos$muestra y colnames(swath)
#     dfGrupos[dfGrupos$Muestras == i, "Grupos"]
#   }))
#   msn$group = factor(msn$group, levels = c(controlG, caseG))
#   
#   
#   t_res1 <- limma_a_b(eset = msn, model.str = "~ group", 
#                       coef.str = "group")
#   t_res1$Protein <- rownames(t_res1)
#   
#   
#   return(t_res1)
# }

# abc <- doContrast(data = swath, dfGrupos = dfGrupos)








