###############################################################################-

############                Assessment by items              #################-

###############################################################################-


# Julia G Curras - 2026/01/27
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")
outDir <-  "C:/Users/julia/Documents/GitHub/normScore/goldStandard/ProcessedDatasets/"

#...........................................................................####
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
library(boot)

source(file = "../R/simulationFunction.R", encoding = "UTF-8")
source(file = "../R/scoreFunction.R", encoding = "UTF-8")

pairwiseAccuracy <- function(manual_rank, score_value) {
  stopifnot(length(manual_rank) == length(score_value))
  n <- length(manual_rank)
  ok <- 0
  tot <- 0
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (is.na(manual_rank[i]) || is.na(manual_rank[j]) ||
          is.na(score_value[i]) || is.na(score_value[j])) next
      
      if (manual_rank[i] == manual_rank[j]) next  # empate manual => no cuenta
      
      tot <- tot + 1
      manual_i_better <- manual_rank[i] < manual_rank[j]
      score_i_better  <- score_value[i] < score_value[j]  # menor = mejor
      
      if (manual_i_better == score_i_better) ok <- ok + 1
    }
  }
  
  if (tot == 0) return(NA_real_)  # todo empatado o faltan datos
  ok / tot
}


#...........................................................................####
# Initial trail ####
## Simulate some data without error ###
set.seed(9396)
data <- simulate_proteomics_clean()

listNorm <- Biomics::doNormalization(
  rawData = data$rawData, 
  logData = data$logData)
listNorm <- listNorm[-9]

## Use normScore to assess it ####
outScore <- normScore(
  normMatrixList = listNorm, 
  designMatrix = data$metadata,
  dfRaw = data$rawData,
  refGroup = "G1", 
  altGroup = "G2", 
  onlyFinalRank = F)
outScore$detailScore
outScore$detailRanking
outScore$finalRanking
outScore$bootstrapScore
outScore$graphic



#...........................................................................####
# Gold Standard ####
dataGS <- readxl::read_excel(path = "_Assessment.xlsx", sheet = "Main", col_names = T)[-1,]

## Apply normScore by dataset (1-25) #
allFiles <- dataGS[1:50, "ID", drop = T]
allFiles <- allFiles[!(allFiles == "PXD005025")]

finalList <- sapply(allFiles, function(i){
  # cat("\t* ", i , "\n")
  output <- readRDS(file = paste0(outDir, i, ".rds"))
  return(normScore(
    normMatrixList = output$listaNorm, 
    designMatrix = output$dm,
    dfRaw = output$data,
    onlyFinalRank = F, 
    onlyDetailRanking = T
    )$detailRanking)
},simplify = F, USE.NAMES = T)

resNormScore <- purrr::imap_dfr(finalList, function(df, dataset_id) {
  df %>%
    dplyr::select(-Total) %>%
    tibble::rownames_to_column("Normalization") %>%   # from rownames to column => normalizations
    tidyr::pivot_longer( # Reorder and make longer
      cols = -Normalization,
      names_to = "Item",
      values_to = "Ranking"
    ) %>%
    dplyr::mutate( # Adding ID dataset
      ID = dataset_id,
      Ranking = as.numeric(Ranking)
    ) %>%
    dplyr::select(ID, Item, Normalization, Ranking) %>%
    dplyr::arrange(ID, Item, Ranking)
})
# resNormScore$ID <- substr(resNormScore$ID, start = 1, stop = nchar(resNormScore$ID)-4)


## Retrieve manual ranking GS ####
items <- paste0("Item", 1:6) # Sheet names = items from score
resGS <- purrr::map_dfr(items, function(sh) {
  # load
  df <- readxl::read_excel(path = "_Assessment.xlsx", sheet = sh) %>%
    dplyr::select(-Number, -Final) %>%
    slice(1:50) %>%
    filter(ID != "PXD005025")# Number col

  # change format
  df %>%
    tidyr::pivot_longer(
      cols = -ID,
      names_to = "Normalization",
      values_to = "RankGS"
    ) %>%
    dplyr::mutate(
      Item = sh,
      RankGS = as.numeric(RankGS)
    ) %>%
    dplyr::select(ID, Item, Normalization, RankGS) %>%
    dplyr::arrange(ID, Item, RankGS)
})



## Compare ####

# Join #
res <- resNormScore %>%
  dplyr::inner_join(resGS,
             by = c("ID", "Item", "Normalization"), )

# Kendall Cor by item and dataset #
kendallCor <- res %>%
  slice(1:1200) %>%
  group_by(ID, Item) %>%
  summarise(
    tau = cor(RankGS,
              Ranking,
              method = "kendall",
              use = "complete.obs"),
    spearman = cor(RankGS,
              Ranking,
              method = "spearman",
              use = "complete.obs"),
    pairAccuracy = pairwiseAccuracy(manual_rank = RankGS, 
                                    score_value = Ranking),
    .groups = "drop"
  ) %>% 
  as.data.frame()

kendallCor %>%
  group_by(Item) %>%
  summarise(
    n_total = n(),
    n_valid = sum(!is.na(tau)),
    prop_valid = n_valid / n_total,

    # ---- Kendall ----
    tau_median = median(tau, na.rm = TRUE),
    tau_mean   = mean(tau, na.rm = TRUE),
    tau_sd     = sd(tau, na.rm = TRUE),
    prop_tau_gt_0_6 = mean(tau > 0.6, na.rm = TRUE),

    # ---- Spearman ----
    spearman_median = median(spearman, na.rm = TRUE),
    spearman_mean   = mean(spearman, na.rm = TRUE),
    spearman_sd     = sd(spearman, na.rm = TRUE),
    prop_spear_gt_0_7 = mean(spearman > 0.7, na.rm = TRUE),

    # ---- Pairwise accuracy ----
    acc_median = median(pairAccuracy, na.rm = TRUE),
    acc_mean   = mean(pairAccuracy, na.rm = TRUE),
    acc_sd     = sd(pairAccuracy, na.rm = TRUE),
    prop_acc_gt_0_8 = mean(pairAccuracy > 0.8, na.rm = TRUE),
    .groups = "drop"
  ) %>% View


ggplot(kendallCor, aes(x = tau)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ Item) +
  theme_minimal()




















