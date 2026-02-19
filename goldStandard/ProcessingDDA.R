###############################################################################-

############    Processing raw data from public datasets      #################-

###############################################################################-


# Julia G Curras - 2026/01/27
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")

#.............................................................................
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
inputDir <- "C:/Users/julia/Documents/GitHub/normScore/goldStandard/Datasets/"
outDir <- "C:/Users/julia/Documents/GitHub/normScore/goldStandard/ProcessedDatasets/"



#.........................................................................####
# PXD055210 ####
idDataset <- "PXD055210"
dfRaw <- read.table(file = "Datasets/PXD055210_proteinGroups.txt", 
                    header = T, sep = "\t")
# Ending in E: EDU | ending in D: DMSO

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs, 
                              starts_with("Intensity.18PC"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant) #, -Intensity)
  # Set 0 to NA #
df[df == 0] <- NA

  # Filters #
df <- df %>% 
    # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
    # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

  # Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL

colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.18", replacement = "")

  ## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = substr(x = colnames(df), start = 4, stop = 4)
)



## Processing ####
## Filtering ###
# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(data, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))






#...........................................................................####
# PXD058699 ####
idDataset <- "PXD058699"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs, 
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant) #, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL

# Group in white
colnames(df)[1] <- "C"
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.rep", replacement = "C.rep")
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)[1:4] <- paste0(colnames(df)[1:4], ".rep0")

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = substr(x = colnames(df), start = 1, stop = 1)
)



## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.1, 0.2, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.2)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))


## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD025933 ####
idDataset <- "PXD025933"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant) #, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins


# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL

# Removing intensity column 
df <- df %>% dplyr::select(-Intensity)


# Group in white
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = "Organoid.", replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = "Stem.Cell", replacement = "StemCell")

allElements <- strsplit(x = colnames(df), split = ".", fixed = T)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = sapply(strsplit(x = colnames(df), split = ".", fixed = T), "[[", 1)
)
dm


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD031710 ####
idDataset <- "PXD031710"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins


# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)


# Group in white
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "G")

allElements <- strsplit(x = colnames(df), split = ".", fixed = T)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = sapply(strsplit(x = colnames(df), split = ".", fixed = T), "[[", 1)
)
dm



## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD066495 ####
idDataset <- "PXD066495"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins


# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "Rep")


## Design matrix ###
tabDM <- readxl::read_xlsx(path = "Datasets/PXD066495_Sample_Legend.xlsx", 
                        sheet = 1, col_names = T)
allElements <- strsplit(x = tabDM$Description, split = " ", fixed = T)
dm <- data.frame(
  Samples = paste0("Rep", sapply(strsplit(tabDM$`File name`, split = "_", fixed = T), "[[", 3)),
  Groups = paste0(sapply(allElements, "[[", 1), "_", sapply(allElements, "[[", 2))
)
dm


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD065898 ####
idDataset <- "PXD065898"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins


# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "S")
colnames(df) <- gsub(x= colnames(df), pattern = "5uM_3h_", replacement = "")


## Design matrix ###
allElements <- strsplit(x = colnames(df), split = "_", fixed = T)
dm <- data.frame(
  Samples = sapply(allElements, "[[", 1),
  Groups = substr(colnames(df), start = 7, stop = 50)
)
dm
colnames(df) <- sapply(allElements, "[[", 1)

# quedámonos con Rabbit e ER que é o que indican na tabla supplementaria 
# (e son os que menos valores faltantes teñen)
dm <- dm %>% filter(Groups %in% c("Nutlin_ER", "Nutlin_IgG_Rabbit", 
                                  "DMSO_ER", "DMSO_IgG_Rabbit")) # Rabbit moitos valore faltantes...
dm <- dm %>% filter(Groups %in% c("Nutlin_ER","DMSO_ER"))
df <- df[, dm$Samples]



## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD064630 ####
idDataset <- "PXD064630"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins


# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")

## Design matrix ###
tabDM <- read.table(file = paste0("Datasets/", idDataset, "_experimentalDesign.txt"), 
                    header = T)
allElements <- strsplit(x = tabDM$PTM, split = "_", fixed = T)
dm <- data.frame(
  Samples = tabDM$PTM,
  Groups = paste0(sapply(allElements, "[[", 1), "_", sapply(allElements, "[[", 2))
)
colnames(df) %in% dm$Samples


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))


#...........................................................................####
# PXD062018 ####
idDataset <- "PXD062018"
dfRaw <- readxl::read_xlsx(path = "Datasets/PXD062018_proteinGroups.xlsx", col_names = T)
colnames(dfRaw) <- gsub(colnames(dfRaw), pattern = " ", replacement = ".")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA
df <- as.data.frame(df)
table(df$Potential.contaminant)
table(df$Only.identified.by.site)
table(df$Reverse)

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(is.na(Only.identified.by.site)) %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(is.na(Potential.contaminant)) %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins


# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = substr(colnames(df), start = 1, stop = 5)
)
dm
colnames(df) %in% dm$Samples



## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD033101 ####
idDataset <- "PXD033101"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")

## Design matrix ###
allElements <- strsplit(x = colnames(df), split = "_", fixed = T)
dm <- data.frame(
  Samples = colnames(df),
  Groups = sapply(allElements, "[[", 1)
)
colnames(df) %in% dm$Samples
table(dm$Groups)

## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.4, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
set.seed(9396)
listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))



#...........................................................................####
# PXD022561 ####
idDataset <- "PXD022561"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
# allElements <- strsplit(x = colnames(df), split = "_", fixed = T)

## Log transforation 
df0 <- df
df <- log(df,  base = 2)

## Design matrix ###
dfReplicates <- data.frame(
  Samples = colnames(df), 
  Replicates = sapply(strsplit(colnames(df), "_"), "[[", 1),
  Groups = rep(c("C", "C", "DM", "DM", "C", "DM", "DM", "C", "C", "DM"), each=4)
)
dm <- data.frame(
  Samples = unique(dfReplicates$Replicates),
  Groups = c("C", "C", "DM", "DM", "C", "DM", "DM", "C", "C", "DM")
)

# Dealing with technical replicates for the sample sample #
df <- Biomics::doJoinReplicates(dfQ = df, dm = dfReplicates, sample_col = "Replicates", 
                          rep_col = "Samples", max_missing_prop = 0.5)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
dataLog <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla
data <- 2^dataLog

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD043324 ####
idDataset <- "PXD043324"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_combined_protein (1).tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.ID,
                              Protein.Probability, 
                              Indistinguishable.Proteins,
                              ends_with(".Intensity"), 
                              -ends_with("Total.Intensity"), 
                              -ends_with("Unique.Intensity"),
                              -ends_with("MaxLFQ.Intensity"), 
                              -starts_with("ZIKV_pool"), 
                              -starts_with("ZIKV_HU"),
                              -starts_with("ZIKV_Pierce")
                              )
colnames(df)

# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Indistinguishable.Proteins == "") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Protein.Probability >0.99) %>% 
  dplyr::select(-Protein.Probability, -Indistinguishable.Proteins) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.ID
df$Protein.ID <- NULL
nrow(df)
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = ".Intensity", replacement = "")
# allElements <- strsplit(x = colnames(df), split = "_", fixed = T)

df0 <- df
df <- log(df, base = 2)


## Replicates ###
dfReplicates <- data.frame(
  Samples = paste0("S", sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 2)),
  Replicates = colnames(df)
)

# Dealing with technical replicates for the sample sample #
df <- Biomics::doJoinReplicates(dfQ = df, dm = dfReplicates, sample_col = "Samples", 
                                rep_col = "Replicates", max_missing_prop = 0.5)


## Design matrix ###
dm <- data.frame(
  Samples = paste0("S", c(1:10, 11:15, 17, 20:21, 23:29)), 
  Groups = c(rep("CTR", 10), rep("CZS+", 8), rep("CZS+", 7))
)

# Seleccionar muestras grupos en df principla
df <- df[, dm$Samples]


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
dataLog <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla
data <- 2^dataLog

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD058626 ####
idDataset <- "PXD058626"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_combined_protein.tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.ID,
                              Protein.Probability, 
                              Indistinguishable.Proteins,
                              ends_with(".Intensity"), 
                              -ends_with("Total.Intensity"), 
                              -ends_with("Unique.Intensity"),
                              -ends_with("MaxLFQ.Intensity") 
)
colnames(df)

# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Indistinguishable.Proteins == "") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Protein.Probability >0.99) %>% 
  dplyr::select(-Protein.Probability, -Indistinguishable.Proteins) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.ID
df$Protein.ID <- NULL
nrow(df)
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = ".Intensity", replacement = "")

allElements <- strsplit(colnames(df), "_", fixed = T)

dm <- data.frame(
  Samples = colnames(df), 
  Groups = gsub("_[0-9]+$", "", colnames(df))
)
dm$GroupsNew <- factor(dm$Groups, levels = unique(dm$Groups), labels = paste0("G", 1:length(unique(dm$Groups))))
dm$SamplesNew <- paste0(dm$GroupsNew, "_", as.integer(sub(".*_([0-9]+)$", "\\1", dm$Samples)))
df <- df[, dm$Samples]
colnames(df) <- dm$SamplesNew
dm <- dm %>% 
  dplyr::select(SamplesNew, GroupsNew) %>% 
  rename(Samples = SamplesNew, 
         Groups = GroupsNew)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
# dataLog <- Biomics::doImputation(df = as.matrix(dataLog))
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD029547 ####
idDataset <- "PXD029547"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_combined_protein.tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.ID,
                              Protein.Probability, 
                              Indistinguishable.Proteins,
                              ends_with(".Intensity"), 
                              -ends_with("Total.Intensity"), 
                              -ends_with("Unique.Intensity"),
                              -ends_with("MaxLFQ.Intensity") 
)
colnames(df)

# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Indistinguishable.Proteins == "") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Protein.Probability >0.99) %>% 
  dplyr::select(-Protein.Probability, -Indistinguishable.Proteins) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.ID
df$Protein.ID <- NULL
nrow(df)
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = ".Intensity", replacement = "")

### Design matrix ###
dm <- read.table(file = paste0("Datasets/", idDataset, "_dm.txt"), header = T)
dm <- dm %>% dplyr::select(sampleName, condition) %>% rename(Samples = sampleName,
                                                             Groups = condition)
dm$Samples <- paste0(toupper(substr(dm$Samples, 1, 1)), substr(dm$Samples, 2, nchar(dm$Samples)))
dm
colnames(df) %in% dm$Samples # Ok



## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
# dataLog <- Biomics::doImputation(df = as.matrix(dataLog))
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))



#...........................................................................####
# PXD038236 ####
idDataset <- "PXD038236"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_combined_protein (1).tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.ID,
                              Protein.Probability, 
                              Indistinguishable.Proteins,
                              ends_with(".Intensity"), 
                              -ends_with("Total.Intensity"), 
                              -ends_with("Unique.Intensity"),
                              -ends_with("MaxLFQ.Intensity") 
)
colnames(df)

# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Indistinguishable.Proteins == "") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Protein.Probability >0.99) %>% 
  dplyr::select(-Protein.Probability, -Indistinguishable.Proteins) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.ID
df$Protein.ID <- NULL
nrow(df)
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = ".Intensity", replacement = "")

### Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(colnames(df), start = 1, stop = 2)
)
dm



## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
# dataLog <- Biomics::doImputation(df = as.matrix(dataLog))
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))



#...........................................................................####
# PXD022614 ####
idDataset <- "PXD022614"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")

## Design matrix ###
allElements <- strsplit(x = colnames(df), split = ".", fixed = T)
dm <- data.frame(
  Samples = colnames(df),
  Groups = c(rep("Adult", 3), rep("AdultSperm", 3), rep("Young", 3))
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
# ## Imputation ###
# set.seed(9396)
# listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD019504 ####
idDataset <- "PXD019504"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = substr(colnames(df), start = nchar(colnames(df)), stop = nchar(colnames(df)))
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.2)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))


#...........................................................................####
# PXD054553 ####
idDataset <- "PXD054553"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
df <- df[,c(7:9, 16:18)] # Only WT samples


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD054636 ####
idDataset <- "PXD054636"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))


#...........................................................................####
# PXD054817 ####
idDataset <- "PXD054817"
dfRaw <- readxl::read_xlsx(
  path = paste0("Datasets/", idDataset,"_proteinGroups.xlsx"), 
  sheet = "MQ_UROtsa", col_names = T
)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(`Protein IDs`,
                              starts_with("Intensity "), #"LFQ.Intensity
                              Reverse, `Only identified by site`,
                              `Potential contaminant`, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(is.na(`Only identified by site`)) %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(is.na(Reverse)) %>% 
  dplyr::filter(is.na(`Potential contaminant`)) %>%
  dplyr::select(-`Only identified by site`, -Reverse, -`Potential contaminant`) # Keep´only interesting proteins
df <- as.data.frame(df)

# Adjusting rownames #
rownames(df) <- df$`Protein IDs`
df$`Protein IDs` <- NULL
rownames(df)[1:4]
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity ", replacement = "")


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")
data <- as.data.frame(data)
dataLog <- as.data.frame(dataLog)


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))



#...........................................................................####
# PXD054817 (II) ####
idDataset <- "PXD054817"
dfRaw <- readxl::read_xlsx(
  path = paste0("Datasets/", idDataset,"_proteinGroups.xlsx"), 
  sheet = "MQ_HBLAK", col_names = T
)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(`Protein IDs`,
                              starts_with("Intensity "), #"LFQ.Intensity
                              Reverse, `Only identified by site`,
                              `Potential contaminant`, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(is.na(`Only identified by site`)) %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(is.na(Reverse)) %>% 
  dplyr::filter(is.na(`Potential contaminant`)) %>%
  dplyr::select(-`Only identified by site`, -Reverse, -`Potential contaminant`) # Keep´only interesting proteins
df <- as.data.frame(df)

# Adjusting rownames #
rownames(df) <- df$`Protein IDs`
df$`Protein IDs` <- NULL
rownames(df)[1:4]
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity ", replacement = "")


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = rep(c("ctrl", "IP"), each = 3)
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")
data <- as.data.frame(data)
dataLog <- as.data.frame(dataLog)


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".2.rds"))



#...........................................................................####
# PXD047566 ####
idDataset <- "PXD047566"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)[1:24] <- paste0("M", colnames(df)[1:24])


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df),
  Groups = rep(c("M3", "M3T", "M6", "M6T", "M8", "M8T", "WT", "WTT"))
)
colnames(df) %in% dm$Samples
table(dm$Groups)


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD012431 ####
idDataset <- "PXD012431"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")


## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD012431_sdrf.xlsx", sheet = 1)
dm <- dm %>% 
  filter(`comment[technical replicate]` == 1) %>%
  mutate(Samples = gsub(x = `comment[data file]`, pattern = "_1.RAW", replacement = "")) %>%
  mutate(Groups = `Characteristics[organism part]`) %>%
  select(Samples, Groups) %>%
  as.data.frame
dm
colnames(df) %in% dm$Samples
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD018900 ####
idDataset <- "PXD018900"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)
grupos <- substr(x = colnames(df),
                 start = 1, 
                 stop = nchar(colnames(df))-7)
# colnames(df) <- substr(x = colnames(df),
#                        start = nchar(colnames(df))-5, 
#                        stop = nchar(colnames(df)))
table(table(colnames(df))>1)


## Design matrix ###
dm <- data.table::fread(file=  "Datasets/PXD018900_experimentalDesignTemplate.txt")
dm <- data.frame(
  Samples = colnames(df), 
  Groups = grupos
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un X% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
set.seed(9396)
listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD019296 ####
idDataset <- "PXD019296"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[", 1:2, simplify = F), paste0, collapse="_")
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD016670 ####
idDataset <- "PXD016670"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = c(rep("Ctl", 3), sapply(strsplit(x = colnames(df)[4:9], split = "_", fixed = T), "[", 2))
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD019139 ####
idDataset <- "PXD019139"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[", 1)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD019103 ####
idDataset <- "PXD019103"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 6, stop = 6)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD029776 ####
idDataset <- "PXD029776"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Only.identified.by.site, -Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df) <- paste0("M", colnames(df))
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD041237 ####
idDataset <- "PXD041237"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df) <- paste0("M", colnames(df))
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD039491 ####
idDataset <- "PXD039491"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD040288 ####
idDataset <- "PXD040288"
dfRaw1 <- data.table::fread(
  file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), check.names = F)
dfRaw2 <- data.table::fread(
  file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), skip = 2868, 
  check.names = F)
colnames(dfRaw2) <- colnames(dfRaw1)
dfRaw <- rbind(dfRaw1, dfRaw2)
colnames(dfRaw) <- gsub(pattern = " ", replacement = ".", x = colnames(dfRaw))


## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 1)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD022038 ####
idDataset <- "PXD022038"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")


## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = nchar(colnames(df))-2)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
set.seed(9396)
listaNorm <- lapply(listaNorm, Biomics::doImputation)


## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD007182 ####
idDataset <- "PXD007182"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")


## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = rep(c("Ctrl", "XAV939"), each = 3)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD042498 ####
idDataset <- "PXD042498"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")


## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = nchar(colnames(df))-3)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))







#...........................................................................####
# PXD038019 ####
idDataset <- "PXD038019"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = nchar(colnames(df))-3)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD069093 ####
idDataset <- "PXD069093"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD069093_KinomeMetadata.xlsx", sheet = 1)
dm <- dm %>% 
  as.data.frame() %>%
  select(File, Treatment) %>% 
  rename(Samples = File, 
         Groups = Treatment) %>%
  mutate(Samples = gsub(x = Samples, pattern = "200928-Pollok-", replacement = "")) %>%
  mutate(Samples = gsub(x = Samples, pattern = ".raw", replacement = "")) %>%
  mutate(Samples = gsub(x = Samples, pattern = "-", replacement = ".")) %>% 
  mutate(Groups = gsub(x = Groups, pattern = " + ", replacement = "_", fixed = T))
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD068908 ####
idDataset <- "PXD068908"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = nchar(colnames(df))-1)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7, 0.9))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

## Log and visualization ###
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
set.seed(9396)
listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD005025 ####
idDataset <- "PXD005025"
dfRaw <- data.table::fread(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"))

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(`Protein IDs`,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, `Only identified by site`,
                              `Potential contaminant`, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(`Only identified by site` != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(`Potential contaminant` != "+") %>%
  dplyr::select(-Reverse, -`Potential contaminant`, -`Only identified by site`) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$`Protein IDs`
df$`Protein IDs` <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.table::fread(file = "Datasets/PXD005025_experimentalDesignTemplate.txt")
dm <- dm %>%
  select(Experiment) %>%
  mutate(Samples = Experiment, 
         Groups = gsub(x = Experiment, pattern = "'", replacement = "", fixed = T)) %>%
  select(Samples, Groups)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

## Log and visualization ###
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
set.seed(9396)
listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD010785 ####
idDataset <- "PXD010785"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = F, sep = "\t")
dfRaw <- data.table::fread(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"))

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(`Protein IDs`,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, `Only identified by site`,
                              `Potential contaminant`, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(`Only identified by site` != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(`Potential contaminant` != "+") %>%
  dplyr::select(-Reverse, -`Potential contaminant`, -`Only identified by site`) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$`Protein IDs`
df$`Protein IDs` <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD010785_samples.xlsx")
dm <- dm %>%
  select(`MaxQuant search label`, treatment) %>%
  rename(Samples = `MaxQuant search label`, 
         Groups = treatment) %>% 
  arrange(Groups) %>%
  as.data.frame()
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD022701 ####
idDataset <- "PXD022701"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = rep(c("C", "E"), each = 3)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD030132 ####
idDataset <- "PXD030132"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")



## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD028841 ####
idDataset <- "PXD028841"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              Potential.contaminant, -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)
df <- df %>% select(-library)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = gsub("[^[:alpha:]]", "", colnames(df))
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

## Log and visualization ###
dataLog <- log(data, base =2)
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
set.seed(9396)
listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))



#...........................................................................####
# PXD015470 ####
idDataset <- "PXD015470"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse,Only.identified.by.site,
                              -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::select(-Reverse, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")




## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD020238 ####
idDataset <- "PXD020238"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, 
                              Potential.contaminant,
                              -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Potential.contaminant != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::select(-Reverse, -Potential.contaminant) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
colnames(df) <- paste0("S", colnames(df))
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.table::fread(file = "Datasets/PXD020238_metadata_NaBut.txt")
dm
dm <- dm %>%
  select(MQ_ID, Treatment) %>%
  rename(Samples = MQ_ID, 
         Groups = Treatment) %>% 
  mutate(Samples = paste0("S", Samples)) %>% as.data.frame
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")




## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))






#...........................................................................####
# PXD020245 ####
idDataset <- "PXD020245"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant,
                              -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Potential.contaminant != "+") %>% 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
colnames(df) <- paste0("S", colnames(df))
table(table(colnames(df))>1)

## Design matrix ###
dm <- data.table::fread(file = "Datasets/PXD020245_metadata_Mouse.txt")
dm <- dm %>%
  mutate(Groups = paste0(Treatment, "_", Tissue_type), 
         Samples = paste0("S", MQ_ID)) %>%
  as.data.frame %>%
  select(Samples, Groups) 
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))









#...........................................................................####
# PXD062946 ####
idDataset <- "PXD062946"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant,
                              -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Potential.contaminant != "+") %>% 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
# dm <- data.table::fread(file = "Datasets/PXD062946_experimentalDesignTemplate.txt")
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[", 1:2, simplify = F), paste0, collapse = "_")
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))









#...........................................................................####
# PXD011804 ####
idDataset <- "PXD011804"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant,
                              -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Potential.contaminant != "+") %>% 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)

## Design matrix ###
# dm <- data.table::fread(file = "Datasets/PXD062946_experimentalDesignTemplate.txt")
dm <- data.frame(
  Samples = colnames(df), 
  Groups = gsub("[^[:alpha:]]", "", colnames(df))
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
dataLog <- log(data, base =2)

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))









#...........................................................................####
# PXD005847 ####
idDataset <- "PXD005847"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_proteinGroups.txt"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.IDs,
                              starts_with("Intensity"), #"LFQ.Intensity
                              Reverse, Only.identified.by.site,
                              Potential.contaminant,
                              -Intensity)
# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Potential.contaminant != "+") %>% 
  dplyr::filter(Only.identified.by.site != "+") %>% 
  dplyr::filter(Reverse != "+") %>% 
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) # Keep´only interesting proteins

# Adjusting rownames #
df <- as.data.frame(df)
rownames(df) <- df$Protein.IDs
df$Protein.IDs <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "Intensity.", replacement = "")
colnames(df)
table(table(colnames(df))>1)


## Replicate matrix ####
df <- log(df, base =2)
dfReplicate <- data.frame(
  Replicates = colnames(df), 
  Samples = substr(x = colnames(df), start = 1, stop = nchar(colnames(df))-2)
)

df <- Biomics::doJoinReplicates(dfQ = df, dm = dfReplicate, 
                          sample_col = "Samples", rep_col = "Replicates")


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
dm
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df[,1:ncol(df)] <- sapply(colnames(df), function(x) as.numeric(df[,x]))
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))
data <- as.data.frame(data)

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
dataLog <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla
data <- 2^dataLog

## Log and visualization ###
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))




#...........................................................................####
# PXD030595 ####
idDataset <- "PXD030595"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_combined_protein.tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.ID,
                              Protein.Probability, 
                              Indistinguishable.Proteins,
                              ends_with("Total.Intensity"), 
                              -ends_with("Razor.Intensity"), 
                              -ends_with("Unique.Intensity"),
                              -ends_with("MaxLFQ.Intensity") 
)
colnames(df)

# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Indistinguishable.Proteins == "") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Protein.Probability >0.99) %>% 
  dplyr::select(-Protein.Probability, -Indistinguishable.Proteins) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.ID
df$Protein.ID <- NULL
nrow(df)
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = ".Total.Intensity", replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = "X2021.5.", replacement = "")
colnames(df) <- substr(x = colnames(df), start = 4, stop = nchar(colnames(df)))

## Design matriz ###
dm <- readxl::read_excel(path = "Datasets/PXD030595_Table_S1.xlsx")
dm <- dm %>%
  select(`Experimental ID...2`, `Class (verbose)`) %>%
  rename(Samples = `Experimental ID...2`, 
         Groups = `Class (verbose)`) %>% 
  as.data.frame
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
# dataLog <- Biomics::doImputation(df = as.matrix(dataLog))
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#...........................................................................####
# PXD056771 ####
idDataset <- "PXD056771"
dfRaw <- read.table(file = paste0("Datasets/", idDataset,"_combined_protein.tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% dplyr::select(Protein.ID,
                              Protein.Probability, 
                              Indistinguishable.Proteins,
                              ends_with(".Intensity"), 
                              -ends_with("Total.Intensity"), 
                              -ends_with("Unique.Intensity"),
                              -ends_with("MaxLFQ.Intensity") 
)
colnames(df)
table(table(df$Protein.ID)>1)
df <- df[!duplicated(df$Protein.ID),]

# Set 0 to NA #
df[df == 0] <- NA

# Filters #
df <- df %>% 
  # 1. Only  identified by site: Remove protein groups in which all peptides have a modified cysteine 
  dplyr::filter(Indistinguishable.Proteins == "") %>% 
  # 2. Reverse: remove protein groups in which some peptide had been found on the reverse library 
  dplyr::filter(Protein.Probability >0.99) %>% 
  dplyr::select(-Protein.Probability, -Indistinguishable.Proteins) # Keep´only interesting proteins

# Adjusting rownames #
rownames(df) <- df$Protein.ID
df$Protein.ID <- NULL
nrow(df)
colnames(df)
colnames(df) <- gsub(x= colnames(df), pattern = ".Intensity", replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = "X", replacement = "G")
colnames(df)


## Design matriz ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
table(dm$Groups)

# Checking
table(colnames(df) %in% dm$Samples)
table(dm$Samples %in% colnames(df))
dm <- dm %>% filter(Samples %in% colnames(df))
df <- df[, dm$Samples]
colnames(df)[which(!(colnames(df) %in% dm$Samples))]
dm$Samples[which(!(dm$Samples %in% colnames(df)))]
nrow(dm) == ncol(df)
table(dm$Groups)
paste0(unique(dm$Groups), collapse = ";")


## Processing ####
## Filtering ###

# Empty rows?
df$Miss <- rowSums(is.na(df))
table(df$Miss == (ncol(df)-1)) # hay filas vacias, fuera
data <- df[which(df$Miss < (ncol(df)-1)), -ncol(df)]
table(rowSums(is.na(data)) == (ncol(data)))

# Checking NA by sample and protein
res <- Biomics::plotBarNA(data = data, interact = F)
res$graficoMuestra
res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.2, 0.5, 0.7, 0.8))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- log(data, base =2)
# dataLog <- Biomics::doImputation(df = as.matrix(dataLog))
min(dataLog, na.rm = "always")
pheatmap::pheatmap(dataLog, show_rownames = F, scale = "row")


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))



