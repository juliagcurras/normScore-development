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
# PXD046983 ####
idDataset <- "PXD046983"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% 
  dplyr::select(
    Protein.Group, 
    D..Data.Andes.20220511_LIMMiceRetina.20220511_RawData_SWATH.20220511_255_OD_LIM_T1.wiff:D..Data.Andes.20220511_LIMMiceRetina.20220511_RawData_SWATH.20220511_336_OS_CTL_T2.wiff)
  
  # Set 0 to NA #
df[df == 0] <- NA

  # Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "D..Data.Andes.20220511_LIMMiceRetina.20220511_RawData_SWATH.20220511_",
  replacement = "")
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = ".wiff",
  replacement = "")
colnames(df) 

  ## Design matrix ###
dataDM <- read.table(file = "Datasets/PXD046983.sdrf.tsv", 
                 header = T, sep = "\t") # Just to check groups: Os_ctl is normal, OD_lim is miopy
table(dataDM$characteristics.disease., dataDM$characteristics.individual.)


## Design matrix ###
dfReplicates <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(colnames(df), "_"), "[[", 3),
  ID = sapply(strsplit(colnames(df), "_"), "[[", 1)
)
dfReplicates$Replicates <- paste0(dfReplicates$ID, "_", dfReplicates$Groups)
View(dfReplicates)

# Dealing with technical replicates for the sample sample #
df <- Biomics::doJoinReplicates(dfQ = df, dm = dfReplicates, sample_col = "Replicates", 
                                rep_col = "Samples", max_missing_prop = 0.5)

dm <- dfReplicates %>% 
  dplyr::select(Replicates, Groups) %>%
  rename(Samples = Replicates)
dm <- dm[duplicated(dm),]



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


#.........................................................................####
# PXD050996 ####
idDataset <- "PXD050996"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"), 
                    header = T, sep = "\t")

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw %>% 
  dplyr::select(
    Protein.Group,
    H..Proteomics.Results.2023.Nov.6.2023.biowires.2nd.expt.CA_1.raw:H..Proteomics.Results.2023.Nov.6.2023.biowires.2nd.expt.MP_21.raw)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "H..Proteomics.Results.2023.Nov.6.2023.biowires.2nd.expt.",
  replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = ".raw", replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = "CP11", replacement = "CP_11") # Error de separación
colnames(df) <- gsub(x= colnames(df), pattern = "CS_5", replacement = "CA_5") # Error de escritura de grupos
colnames(df) 

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(colnames(df), split = "_", fixed = T), "[[", 1)
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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla

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




#.........................................................................####
# PXD052720 ####
idDataset <- "PXD052720"
# dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
#                     header = T, sep = "\t")
dfRaw <- readxl::read_xlsx(path = "Datasets/PXD052720_report.xlsx", sheet = 1, col_names = T)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "D:\\Project\\ZZY\\XA05851DA\\XA05851DA_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".raw", replacement = "")
colnames(df) 

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(colnames(df), split = "_", fixed = T), "[[", 1)
)
colnames(df) %in% dm$Samples
table(dm$Groups) # sy = dengue patients

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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.1, 0.3, 0.5))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

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




#.........................................................................####
# PXD062678 ####
idDataset <- "PXD062678"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "G..Avinash_timsTOF.Saliva_timsTOF.",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".d", replacement = "")
colnames(df) <- gsub(x= colnames(df), pattern = "Saliva_Sample_rawFiles.", replacement = "", fixed = T)
colnames(df)

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD062678_Metadata.file.xlsx")
dm <- as.data.frame(dm)
table(dm$Group)
dm <- dm %>% 
  filter(Group %in% c("Control", "OSCC", "Pre-Malignant Lesions")) %>% 
  rename(Samples = `File name`)
dm$Samples <- gsub(x= dm$Samples, pattern = ".d", replacement = "")
dm$Samples <- gsub(x= dm$Samples, pattern = "-", replacement = ".", fixed = T)
df <- df[, dm$Samples]
dm$Samples <- sapply(strsplit(dm$Samples, split = "_", fixed = T), "[[", 1)
colnames(df) <- sapply(strsplit(colnames(df), split = "_", fixed = T), "[[", 1)


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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla

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


#.........................................................................####
# PXD063236 ####
idDataset <- "PXD063236"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 5:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(x= colnames(df), pattern = "S..Proteomics.CMB.CMB.1514.",
                     replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = "_EvoAurEl6_20SPDDIAPASEF.CMB.1514",
                     replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = "_EvoAurEl6_20SPDDIAPASEF.CMB.1314",
                     replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".d", replacement = "", fixed = T)
colnames(df) <- sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 3)
colnames(df) <- gsub(x= colnames(df), pattern = "S2.", replacement = "", fixed = T)
colnames(df)



## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 1)
)
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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.3)$tabla

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





#.........................................................................####
# PXD063383 ####
idDataset <- "PXD063383"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "E..FelixW.AGSchaefer_BV2_PU1.20230629_E1_EasyLC1_CollID_BN065_50cm_FW_Collab02_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".raw", replacement = "", fixed = T)
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
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



#.........................................................................####
# PXD063923 ####
idDataset <- "PXD063923"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "D..Data.Matthias.DIA.RAW.FIle.230526_MG_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".raw", replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = "RX_1710", replacement = "Rx_1710", fixed = T)
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
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




#.........................................................................####
# PXD058791 ####
idDataset <- "PXD058791"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "F..SJ.202409.24092507_SML_",
  replacement = "", fixed = T)
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "F..SJ.202409.24092508_SML_",
  replacement = "", fixed = T)
colnames(df) <- sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[", 1)
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = nchar(colnames(df))-1)
)
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




#.........................................................................####
# PXD058655 ####
idDataset <- "PXD058655"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "H..Tyler_Cooper.Kathleen_NatureCellBio_C48Protoemics.GPF.DIA_June2024.tc_Kathleen_GPFDIA_60min_",
  replacement = "", fixed = T)
df <- df %>% select(-starts_with("H"))
colnames(df) <- sapply(strsplit(x = colnames(df), split = ".mzML", fixed = T), "[", 1)
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1)
)
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



#.........................................................................####
# PXD051732 ####
idDataset <- "PXD051732"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "L..MSdata.timsTOF.Data.TimoRisch.2024.03.17.Clinical_isolates_CYS_mutants.",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".d", replacement = "", fixed = T)
colnames(df) <- sapply(sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[", 1:3, simplify = F), paste0, collapse = "_")
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[", 1:2, simplify = F), paste0, collapse = "_")
)
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




#.........................................................................####
# PXD055964 ####
idDataset <- "PXD055964"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
# dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "D..Mass_spectrometry.Raw_data.Joakim.Lund.",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".d", replacement = "", fixed = T)
colnames(df)

## Design matrix ###
dm <- readxl::read_excel(path = paste0("Datasets/", idDataset, "_Metadata.xlsx"), sheet = 1)
dm <- as.data.frame(dm)
dm$Groups <- paste0(dm$`i8-SRF`, "_", dm$`ANGii treatment`)
table(dm$Groups)
colnames(df) %in% dm$RawfileID
dm <- dm %>% select(RawfileID, Groups) %>% rename(Samples = RawfileID)
table(dm$Groups)
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



#.........................................................................####
# PXD060381 ####
idDataset <- "PXD060381"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "Z..2.ASTRAL.Files.2024.Puvanesarajah..Varun.24.060.Martinez_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = "_24.060.raw", replacement = "", fixed = T)
colnames(df)

## Design matrix ###

dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
dm <- as.data.frame(dm)
table(dm$Groups)
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




#.........................................................................####
# PXD056859 ####
idDataset <- "PXD056859"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "C..Data.Ref.962.combined.analysis.",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = ".raw", replacement = "", fixed = T)
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
dm <- as.data.frame(dm)
dm$Groups <- factor(dm$Groups, levels = c("WT", "He", "Ho"), 
                    labels = c("WT", "Het", "Hom"))
table(dm$Groups)
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




#.........................................................................####
# PXD047128 ####
idDataset <- "PXD047128"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "C..sergo.amber_hart.2212.crude.urine.221220_",
  replacement = "", fixed = T)
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "Orbi3_SK_SER_G60_T85_Sheffield_Hart_Sample_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x= colnames(df), pattern = "_DIA_1ug.raw", replacement = "", fixed = T)
colnames(df) <- sapply(strsplit(x = colnames(df), split = "_"), "[", 2)
colnames(df)


## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = substr(x = colnames(df), start = 1, stop = 2)
)
dm <- as.data.frame(dm)
table(dm$Groups)
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




#.........................................................................####
# PXD064948 ####
idDataset <- "PXD064948"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 5:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA


# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "Z..MS.backup.TIMS.ToF.HT.2024.10_October.20241003_SR_chcondocytes_",
  replacement = "", fixed = T)
colnames(df) <- paste0("S", sapply(strsplit(x = colnames(df), split = "_"), "[", 1))

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = rep(c("Unfixed", "Fixed"), each = 4) # Fig S2 suppl mat (in boxplots sample 1 to 4 constitute one group)
)
dm <- as.data.frame(dm)
table(dm$Groups)
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





#.........................................................................####
# PXD051201 ####
idDataset <- "PXD051201"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA


# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "I..gespadas.FATR.2021MK013_FATR_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x = colnames(df), pattern = "_01_2ug.raw", replacement = "")
colnames(df)

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD051201_DM.xlsx", sheet = 1, col_names = T)
dm <- as.data.frame(dm[, -3])
table(dm$Groups)
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
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

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



#.........................................................................####
# PXD068597 ####
idDataset <- "PXD068597"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 7:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA


# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "C..Users.whuang.Documents.P1959.P1959_1ug_",
  replacement = "", fixed = T)
colnames(df) <- gsub(x = colnames(df), pattern = "_MN_07.09.2025.raw", replacement = "")
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_"), "[[", 1)
)
table(dm$Groups)
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



#.........................................................................####
# PXD060921 ####
idDataset <- "PXD060921"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA


# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "Y..timsTOF.Flex.Data2024.Etienne.Toulouse.HCT.Sample.",
  replacement = "", fixed = T)
colnames(df) <- gsub(x = colnames(df), pattern = "diff_", replacement = "")
colnames(df) <- sapply(sapply(strsplit(x = colnames(df), split = "_"), "[", 1:2, simplify = F), paste0, collapse = "_")
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(x = colnames(df), split = "_"), "[[", 1)
)
table(dm$Groups)
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



#.........................................................................####
# PXD071471 ####
idDataset <- "PXD071471"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "Y..timsTOF.Flex.Data2024.Etienne.Toulouse.HCT.Sample.",
  replacement = "", fixed = T)

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD071471_SDRF_pregnancy.xlsx", sheet = 1)
dm <- dm %>% 
  dplyr::select(`comment[data file]`, `characteristics[disease]`) %>%
  rename(Samples = `comment[data file]`, Groups = `characteristics[disease]`)
table(dm$Groups)
dm$Groups <- factor(dm$Groups)
levels(dm$Groups) <- c("AR_3rdT", "AH_6MPP", "H_3rdT", "H_6MPP")
table(dm$Groups)
dm <- dm %>% filter(Samples %in% colnames(df))
dm <- as.data.frame(dm)
dm$Samples <- as.character(dm$Samples)
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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.7)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- as.matrix(log(data, base =2))

## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
listaNorm <- lapply(listaNorm, Biomics::doImputation)


## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))





#.........................................................................####
# PXD057069 ####
idDataset <- "PXD057069"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x= colnames(df), 
  pattern = "F..Hazal.HvP.R0941.hippocampus.raw...DIA.20230907_TIMS6_EVO4_PRI_X1050_P0146_R0941_44min_DIA_",
  replacement = "", fixed = T)
colnames(df) <- sapply(strsplit(x = colnames(df), split = "_"), "[[", 1)
colnames(df)

## Design matrix ###
dm <- read.table(file = "Datasets/PXD057069_HvP_Metadata.txt", header = T)
dm$Group <- paste0(dm$genotype, "_", dm$treatment)
table(dm$Group)
dm <- dm %>% select(sample_id, Group) %>%
  rename(Samples = sample_id, Groups = Group)
colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
dm$Samples <- paste0("S", dm$Samples)
colnames(df) <- paste0("S", colnames(df))


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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = .3)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- as.matrix(log(data, base =2))
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




#.........................................................................####
# PXD055158 ####
idDataset <- "PXD055158"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv.txt"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 5:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL

## Design matrix ###
dm <- read.table(file = "Datasets/PXD055158_experiment_annotation.txt", header = T)
dm <- dm %>% select(sample_name, condition) %>%
  rename(Samples = sample_name, Groups = condition)
dm
colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)


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
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = .5)$tabla

## Log and visualization ###
pheatmap::pheatmap(data, show_rownames = F)
dataLog <- as.matrix(log(data, base =2))
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




#.........................................................................####
# PXD045557 ####
idDataset <- "PXD045557"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "Y..HF_X_2022.Karakosta_Cassett_30446.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Replicates ###
dmReplicates <- data.frame(
  Samples = sapply(strsplit(colnames(df), "_"), "[[", 1), 
  Replicates = colnames(df)
)
dmReplicates
df <- Biomics::doJoinReplicates(dfQ = df, dm = dmReplicates, sample_col = "Samples",
                          rep_col = "Replicates", max_missing_prop = 0.5)


## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD045557_KEY_Phaco.xlsx", sheet = 1)
dm <- dm %>% select(PHACO, ...1) %>%
  rename(Samples = PHACO, Groups = ...1)
dm <- as.data.frame(dm)
dm
colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
dm <- dm %>% filter(Samples %in% colnames(df))
nrow(dm) == ncol(df)
table(dm$Groups)

colnames(df) <- paste0("S", colnames(df))
dm$Samples <- paste0("S", dm$Samples)


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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = .7)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))


## Normalization ###
listaNorm <- Biomics::doNormalization(rawData = data, logData = dataLog,
                                      listaNorm = c("Mean", "Median", "TI", "VSN", 
                                                    "Quantile", "CyclicLoess", "RLR"))
## Imputation ###
listaNorm <- lapply(listaNorm, Biomics::doImputation)

## Ouput ###
output <- list(
  data = data, 
  dataLog = dataLog, 
  listaNorm = listaNorm,
  dm = dm
)

saveRDS(object = output, file = paste0(outDir, idDataset, ".rds"))


#.........................................................................####
# PXD045554 ####
idDataset <- "PXD045554"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "Y..HF_X_2022.Karakosta_AH_30449.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Replicates ###
dmReplicates <- data.frame(
  Samples = sapply(strsplit(colnames(df), "_"), "[[", 1), 
  Replicates = colnames(df)
)
dmReplicates
df <- Biomics::doJoinReplicates(dfQ = df, dm = dmReplicates, sample_col = "Samples",
                          rep_col = "Replicates", max_missing_prop = 0.5)


## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD045554_KEY_AH.xlsx", sheet = 1)
dm <- dm %>% select(AH, ...1) %>%
  rename(Samples = AH, Groups = ...1)
dm <- as.data.frame(dm)
dm
colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
dm <- dm %>% filter(Samples %in% colnames(df))
nrow(dm) == ncol(df)
table(dm$Groups)

colnames(df) <- paste0("S", colnames(df))
dm$Samples <- paste0("S", dm$Samples)


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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))


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



#.........................................................................####
# PXD054270 ####
idDataset <- "PXD054270"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "D..Proteomics2023.ProteomicsUnit2023.KULeuven.MarcFransenLab.MCF001695.15Samples.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".wiff", replacement = "")
colnames(df)

## Design matrix ###
dm <- read.table(file = "Datasets/PXD054270_DIA_experimental_design_MCF001695.txt", header = T)
dm <- dm %>% select(-replicate) %>%
  rename(Samples = label, Groups = condition)
dm <- as.data.frame(dm)
dm
colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
dm <- dm %>% filter(Samples %in% colnames(df))
nrow(dm) == ncol(df)
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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))


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




#.........................................................................####
# PXD040451 ####
idDataset <- "PXD040451"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "D..Zoidakis_Filtered_30330.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Replicate matrix ###

dfReplicates <- data.frame(
  Samples = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1),
  Replicates = colnames(df)
)

df <- Biomics::doJoinReplicates(dfQ = df, dm = dfReplicates, sample_col = "Samples", rep_col = "Replicates")

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD040451_Key_30330.xlsx", sheet = 1, col_names = F)
colnames(dm) <- c("Samples", "Groups")
dm <- as.data.frame(dm)

colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
dm <- dm %>% filter(Samples %in% colnames(df))
nrow(dm) == ncol(df)
table(dm$Groups)

colnames(df) <- paste0("S", colnames(df))
dm$Samples <- paste0("S", dm$Samples)


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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0.5)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))


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




#.........................................................................####
# PXD051606 ####
idDataset <- "PXD051606"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "D..Dylan.Mitopep.plasma.20230731_Exploris480_DH_AD_Plasma_", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Design matrix ###
dm <- data.frame(
  Samples = colnames(df), 
  Groups = sapply(strsplit(colnames(df), split = "_", fixed = T), "[[", 1)
)
dm

colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))
pheatmap::pheatmap(dataLog)

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





#.........................................................................####
# PXD043635 ####
idDataset <- "PXD043635"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "E..Project_633_Zwart_Stefan.Seppak.20221110_Exploris_Evosep_EV1137_project_633_Seppak_Proteome_15spd_DIA_", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df) <- paste0("WZ", colnames(df))
colnames(df)

## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD043635_sample_info.xlsx", sheet = 1)
dm <- as.data.frame(dm)
dm$Samples <- gsub(pattern = "\tMM", x = dm$Samples, replacement = "")


colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))
pheatmap::pheatmap(dataLog)

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








#.........................................................................####
# PXD044220 ####
idDataset <- "PXD044220"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "Y..HF_X_2023.Dafou_40490.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Replicates matrix ###
dmReplicates <- data.frame(
  Samples = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1), 
  Replicates = colnames(df)
  
)

df <- Biomics::doJoinReplicates(
  dfQ = df, dm = dmReplicates, sample_col = "Samples", rep_col = "Replicates")


## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD044220_Dafou_KEY.xlsx", sheet = 1)
dm <- as.data.frame(dm)
dm <- dm %>% rename(Samples = `RawFile Number`, Groups = `Sample Description`)
dm$Groups <- rep(c("WT_D0", "WT_D3", "WT_D4", "KO_D0", "KO_D3", "KO_D4"), each = 3)
dm

colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
table(dm$Groups)

colnames(df) <- paste0("S", colnames(df))
dm$Samples <- paste0("S", dm$Samples)


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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))
pheatmap::pheatmap(dataLog)

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







#.........................................................................####
# PXD054374 ####
idDataset <- "PXD054374"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "D..Proteomics2024.KULeuven.MarcFrasenLab.PCF000078.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".wiff", replacement = "")
colnames(df)

## Design matrix ###
dm <- read.table(file = "Datasets/PXD054374_DIA_experimental_design_PCF000078.txt", header = T)
dm <- as.data.frame(dm)
dm <- dm %>% rename(Samples = label, Groups = condition) %>% select(-replicate)
dm

colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))
pheatmap::pheatmap(dataLog, scale = "row")

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







#.........................................................................####
# PXD047024 ####
idDataset <- "PXD047024"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "Y..HF_X_2023.Verginis_40485_CAFs.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Replicates matrix ####
dmReplicates <- data.frame(
  Samples = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1), 
  Replicate = colnames(df)
)
dmReplicates
df <- Biomics::doJoinReplicates(dfQ = df, dm = dmReplicates, )


## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD047024_40485_KEY.xlsx", sheet = 1, col_names = F)
dm <- as.data.frame(dm)
dm <- dm %>% rename(Samples = ...1, Groups = ...2) 
dm$Groups <- sapply(strsplit(x = dm$Groups, split = "_"), "[", 1)

colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
table(dm$Groups)

colnames(df) <- paste0("S", colnames(df))
dm$Samples <- paste0("S", dm$Samples)


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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))
pheatmap::pheatmap(dataLog, scale = "row")

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






#.........................................................................####
# PXD047024 ####
idDataset <- "PXD047024"
dfRaw <- read.table(file = paste0("Datasets/", idDataset, "_report.pg_matrix.tsv"),
                    # dfRaw <- read.table(file = "Datasets/report.pg_matrix.tsv",
                    header = T, sep = "\t")
colnames(dfRaw)

## Datasets ####
## Quantification matrix ###
# Selecting interesting cols #
df <- dfRaw[,c(1, 6:ncol(dfRaw))]
df <- as.data.frame(df)

# Set 0 to NA #
df[df == 0] <- NA

# Adjusting rownames #
rownames(df) <- df$Protein.Group
df$Protein.Group <- NULL
colnames(df) <- gsub(
  x = colnames(df), 
  pattern = "Y..HF_X_2023.Verginis_40485_CAFs.", 
  replacement = ""
)
colnames(df) <- gsub(x = colnames(df), pattern = ".raw", replacement = "")
colnames(df)

## Replicates matrix ####
dmReplicates <- data.frame(
  Samples = sapply(strsplit(x = colnames(df), split = "_", fixed = T), "[[", 1), 
  Replicate = colnames(df)
)
dmReplicates
df <- Biomics::doJoinReplicates(dfQ = df, dm = dmReplicates, )


## Design matrix ###
dm <- readxl::read_excel(path = "Datasets/PXD047024_40485_KEY.xlsx", sheet = 1, col_names = F)
dm <- as.data.frame(dm)
dm <- dm %>% rename(Samples = ...1, Groups = ...2) 
dm$Groups <- sapply(strsplit(x = dm$Groups, split = "_"), "[", 1)

colnames(df) %in% dm$Samples
dm$Samples %in% colnames(df)
nrow(dm) == ncol(df)
table(dm$Groups)

colnames(df) <- paste0("S", colnames(df))
dm$Samples <- paste0("S", dm$Samples)


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
# res$graficoProteina # proteinas con moitos valores faltantes

# Applying different thresholds for filtering
resFilt <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = c(0, 0.3, 0.5, 0.7))
resFilt$tablaFormato

# Finally: permitimos hasta un 30% de valores faltantes por grupo
data <- Biomics::filterMissing(df = data, dfGrupos = dm, threshold = 0)$tabla

# Imputation & log-transformation
dataLog <- as.matrix(log(data, base =2))
pheatmap::pheatmap(dataLog, scale = "row")

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
