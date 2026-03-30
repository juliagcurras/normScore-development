###############################################################################-

############                Assessment by items              #################-

###############################################################################-


# Julia G Curras - 2026/01/27
rm(list=ls())
graphics.off()
outDir <-  "C:/Users/julia/Documents/GitHub/normScore/goldStandard/ProcessedDatasets/"
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")

#...........................................................................####
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
library(boot)

source(file = "../R/simulationFunction.R", encoding = "UTF-8")
source(file = "../R/scoreFunction.R", encoding = "UTF-8")
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")

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

rankMethods <- function(score_named_vec) {
  # score_named_vec: numeric, con names = métodos
  stopifnot(!is.null(names(score_named_vec)))
  names(sort(score_named_vec, decreasing = FALSE, na.last = NA))
}

hitAtK<- function(
    rankMethods, # normScore result
    targetSet, # normScore result
    k = 1, 
    type = c("top", "bottom")) {
  
  if (length(rankMethods) == 0 || length(targetSet) == 0)
    return(NA_integer_)
  
  if (type == "top") {
    selected <- head(rankMethods, k)
  } else if (type == "bottom") {
    selected <- tail(rankMethods, k)
  } else {
    stop("Error! Revisa el código")
  }
  
  return(as.integer(any(selected %in% targetSet)))
}


mrrAtK <- function(
    rankedMethods,  # normScore result
    targetSet, # normScore result
    type = c("top", "bottom")) {
  
  type <- match.arg(type)# fuerza a que solo sea "top" o "bottom"
  
  if (length(rankedMethods) == 0 || length(targetSet) == 0) {
    return(NA_real_)
  }
  
  n <- length(rankedMethods)
  pos <- match(targetSet, rankedMethods)   # posiciones en el ranking (NA si no está)
  pos <- pos[!is.na(pos)]
  
  if (length(pos) == 0) {
    return(0) # ninguno de los targetSet aparece (debería ser raro si todo está bien)
  }
  
  if (type == "top") {
    1 / min(pos)
  } else if (type == "bottom") {
    # posición desde abajo: 1 = último, 2 = penúltimo, ...
    posFromBottom <- n - pos + 1
    1 / min(posFromBottom)
  } else {
    stop("Error! Revisa el código")
  }
}


bottomCoverage <- function(rankedMethods, worstSet, k = NULL, prop = FALSE) {
  
  if (length(rankedMethods) == 0 || length(worstSet) == 0) {
    return(NA_real_)
  }
  
  # Si no se especifica k, usar tamaño del worstSet
  if (is.null(k)) {
    k <- length(worstSet)
  }
  
  # Evitar k mayor que número de métodos
  k <- min(k, length(rankedMethods))
  
  bottomK <- tail(rankedMethods, k)
  
  nCovered <- sum(bottomK %in% worstSet)
  
  if (prop) {
    return(nCovered / k)
  } else {
    return(nCovered)
  }
}

summaryHit <- function(res){
  df <- res %>%
    as.data.frame %>%
    rename(k1 = V1, k2 = V2, k3 = V3)
  # tabla <- df %>%
  #     summarise(
  #     n = n(),
  #     nHitAt1 = sum(k1, na.rm = TRUE),
  #     propHitAt1 = mean(k1, na.rm = TRUE),
  #     nHitAt2 = sum(k2, na.rm = TRUE),
  #     propHitAt2 = mean(k2, na.rm = TRUE),
  #     nHitAt3 = sum(k3, na.rm = TRUE),
  #     propHitAt3 = mean(k3, na.rm = TRUE)
  #   )
  return(Biostatech::getCatTable(df)$tablaFormato)
  # tabla <- Biostatech::getNumTable(df)$tablaFormato
  
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

## Apply normScore by dataset (1-25) ####
allFiles <- dataGS[1:25, "ID", drop = T]
# allFiles <- allFiles[!(allFiles == "PXD005025")]

finalList <- sapply(allFiles, function(i){
  cat("\t* ", i , "\n")
  output <- readRDS(file = paste0(outDir, i, ".rds"))
  return(normScore(
    normMatrixList = output$listaNorm, 
    designMatrix = output$dm,
    dfRaw = output$data,
    onlyFinalRank = F, 
    onlyDetailRanking = T
    )$detailRanking)
},simplify = F, USE.NAMES = T)
# saveRDS(object = finalList, file = "AssessmentFiles/normScore_dataGS_1_To_25.rds")
# finalList <- readRDS(file = "AssessmentFiles/normScore_dataGS_1_To_25.rds")
finalList <- readRDS(file = "AssessmentFiles/normScore_CORRECTED_dataGS_1_To_25.rds")


### Info for item assessment ####
# change wid format to longer format for comparison by item
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

### Info for global assessment ####
resGlobalNormScore <- sapply(finalList, rownames,
                             simplify = F, USE.NAMES = T)
  # df %>% rownames() # quitar normalization + puntuación score
    # tibble::rownames_to_column("Normalization") %>%
    # dplyr::select(Normalization, Total) %>% 
    # tibble::deframe()
  # return(finalVect)




## Retrieve manual ranking GS ####
### BY ITEM  ####
items <- paste0("Item", 1:6) # Sheet names = items from score
resGS <- purrr::map_dfr(items, function(sh) {
  # load
  df <- readxl::read_excel(path = "_Assessment.xlsx", sheet = sh) %>%
    dplyr::select(-Number, -Final) #%>%
    # slice(1:50) %>%
    # filter(ID != "PXD005025")# Number col

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


### GLOBAL  ####
normalizationNames <- c("Log", "Mean", "Median","TI", "VSN",
                        "Quantile", "CyclicLoess", "RLR", "MAD")
dfGS_Main <- readxl::read_excel(path = "_Assessment.xlsx", sheet = "Main")[-1,] 

# Best UNIQUE normalization
dfBest <- dfGS_Main %>%
    dplyr::select(ID, normScore) %>%
    tibble::deframe()
resGlobalGS_BEST <- sapply(dfBest, strsplit, split = ";", fixed = T,
                    simplify = T, USE.NAMES = T)
resGlobalGS_BEST <- sapply(resGlobalGS_BEST, "[[", 1, simplify = F)
if (!all(sapply(resGlobalGS_BEST, function(i) any(i %in% normalizationNames)))){
    stop("Algún nombre de normalization no coincide!!")
}

# Best normalizations
dfBestS <- dfGS_Main %>%
    dplyr::select(ID, normScore) %>%
    tibble::deframe()
resGlobalGS_BESTS <- sapply(dfBest, strsplit, split = ";", fixed = T,
                    simplify = T, USE.NAMES = T)
if (!all(sapply(resGlobalGS_BESTS, function(i) any(i %in% normalizationNames)))){
    stop("Algún nombre de normalization no coincide!!")
}

# Excluded normalizations
dfWorst <- dfGS_Main %>%
  dplyr::select(ID, Excluded) %>%
  filter(!is.na(Excluded)) %>%
  tibble::deframe()
resGlobalGS_WORST <- sapply(dfWorst, strsplit, split = ";", fixed = T,
                    simplify = T, USE.NAMES = T)
if (!all(sapply(resGlobalGS_WORST, function(i) any(i %in% normalizationNames)))){
    stop("Algún nombre de normalization no coincide!!")
}


allResultsGS <- list(
  byItem = resGS, 
  globalBest = resGlobalGS_BEST, 
  globalBests = resGlobalGS_BESTS, 
  globalWorst = resGlobalGS_WORST
)

saveRDS(allResultsGS, file = "AssessmentFiles/results_Gold_Standard.rds")



#...........................................................................####
# Comparison by item ####

# Join #
res <- resNormScore %>%
  dplyr::inner_join(resGS,
             by = c("ID", "Item", "Normalization"), )

# Kendall Cor by item and dataset #
kendallCor <- res %>%
  slice(1:1200) %>%
  group_by(ID, Item) %>%
  summarise(
    ## ---- Kendall ----
    tau = cor(RankGS,
              Ranking,
              method = "kendall",
              use = "complete.obs"),
    ## ---- Spearman ----
    spearman = cor(RankGS,
              Ranking,
              method = "spearman",
              use = "complete.obs"),
    ## ---- Pairwise accuracy ----
    pairAccuracy = pairwiseAccuracy(manual_rank = RankGS, 
                                    score_value = Ranking),
    .groups = "drop"
  ) %>% 
  as.data.frame()

# Summarize results #
    # Table #
tableByItem <- kendallCor %>%
  group_by(Item) %>%
  summarise(
    n_total = n(),
    n_valid = sum(!is.na(tau)),
    prop_valid = n_valid / n_total,

    # ---- Kendall ---
    tau_median = median(tau, na.rm = TRUE),
    tau_mean   = mean(tau, na.rm = TRUE),
    tau_sd     = sd(tau, na.rm = TRUE),
    prop_tau_gt_0_6 = mean(tau > 0.6, na.rm = TRUE),

    # ---- Spearman ---
    spearman_median = median(spearman, na.rm = TRUE),
    spearman_mean   = mean(spearman, na.rm = TRUE),
    spearman_sd     = sd(spearman, na.rm = TRUE),
    prop_spear_gt_0_7 = mean(spearman > 0.7, na.rm = TRUE),

    # ---- Pairwise accuracy ---
    acc_median = median(pairAccuracy, na.rm = TRUE),
    acc_mean   = mean(pairAccuracy, na.rm = TRUE),
    acc_sd     = sd(pairAccuracy, na.rm = TRUE),
    prop_acc_gt_0_8 = mean(pairAccuracy > 0.8, na.rm = TRUE),
    .groups = "drop"
  )
      # Graphical representation #
ggplot(kendallCor, aes(x = tau)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ Item) +
  theme_minimal()




#...........................................................................####
# Global Comparison ####

## Hit top ####
  # Best (only one) => Does first position from GS and normScore match with 1, 2, 3 position?
bestSet <- resGlobalGS_BEST[names(resGlobalNormScore)]
resHitBestAll <- sapply(1:3, function(j) {sapply(names(resGlobalNormScore), function(i) 
  hitAtK(
    rankMethods = resGlobalNormScore[[i]], 
    targetSet = bestSet[[i]],   
    k = j, 
    type = "top"))
  })
resHitBestAll
summaryHit(res = resHitBestAll)

  # Best (set of normalizations) => Do some of the valid normalizations from GS and the selected by the normScore ( for the 1, 2, 3 position) match?
bestSet <- resGlobalGS_BESTS[names(resGlobalNormScore)]
resHitBestsAll <- sapply(1:3, function(j) {sapply(names(resGlobalNormScore), function(i) 
  hitAtK(
    rankMethods = resGlobalNormScore[[i]], 
    targetSet = bestSet[[i]],   
    k = j, 
    type = "top"))
  })
resHitBestsAll
summaryHit(res = resHitBestsAll)

  # worst - control negativo, debería ser 0 todo
worstToAsses <- intersect(names(resGlobalGS_WORST), names(resGlobalNormScore))
worstSetGS <- resGlobalGS_WORST[worstToAsses]
worstSetNS <- resGlobalNormScore[worstToAsses]
resHitWorstAll <- sapply(1:3, function(j) {sapply(names(worstSetGS), function(i) 
  hitAtK(
    rankMethods = worstSetNS[[i]], 
    targetSet = worstSetGS[[i]],   
    k = j, 
    type = "top"))
})
resHitWorstAll
summaryHit(res = resHitWorstAll)



## Excluded normalizations ####
worstToAsses <- intersect(names(resGlobalGS_WORST), names(resGlobalNormScore))
worstSetGS <- resGlobalGS_WORST[worstToAsses]
worstSetNS <- resGlobalNormScore[worstToAsses]
resCovWorstAll <- sapply(names(worstSetGS), function(i) 
  bottomCoverage(
    rankedMethods = worstSetNS[[i]], 
    worstSet = worstSetGS[[i]], 
    # k = 5,
    prop = T))
table(resCovWorstAll)



## MRR ####
# Best (only one) => Does first position from GS and normScore match with 1, 2, 3 position?
bestSet <- resGlobalGS_BEST[names(resGlobalNormScore)]
resMrrBestAll <- sapply(names(resGlobalNormScore), function(i) 
  mrrAtK(
    rankedMethods = resGlobalNormScore[[i]], 
    targetSet = bestSet[[i]],   
    type = "top"))
summary(resMrrBestAll)

# Best (set of normalizations) => Do some of the valid normalizations from GS and the selected by the normScore ( for the 1, 2, 3 position) match?
bestSet <- resGlobalGS_BESTS[names(resGlobalNormScore)]
resMrrBestsAll <- sapply(names(resGlobalNormScore), function(i) 
  mrrAtK(
    rankedMethods = resGlobalNormScore[[i]], 
    targetSet = bestSet[[i]],   
    type = "top"))
table(resMrrBestsAll)

# worst
worstToAsses <- intersect(names(resGlobalGS_WORST), names(resGlobalNormScore))
worstSetNS <- resGlobalGS_WORST[worstToAsses]
worstSetGS <- resGlobalNormScore[worstToAsses]
resMrrWorstAll <- sapply(names(worstSetGS), function(i) 
  mrrAtK(
    rankedMethods = worstSetNS[[i]], 
    targetSet = worstSetGS[[i]],   
    type = "bottom")
  )
table(resMrrWorstAll)


## All
output <- list(
  tableByItem, 
  resHitBestAll, 
  resHitBestsAll, 
  resHitWorstAll,
  resCovWorstAll,
  resMrrBestAll,
  resMrrBestsAll,
  resMrrWorstAll
)

# saveRDS(object = output, file = "AssessmentFiles/normScore_dataGS_1_To_25_OUTPUT.rds")
# saveRDS(object = output, file = "AssessmentFiles/normScore_dataGS_CORRECTED_1_To_25_OUTPUT.rds")
output <- readRDS(file = "AssessmentFiles/normScore_dataGS_1_To_25_OUTPUT.rds")




#...........................................................................####
# Function - Format normScore results ####

formatNormScoreResults <- function(listByItem, finalList){
  
  ### Info for item assessment 
  # change wid format to longer format for comparison by item
  resNormScore <- purrr::imap_dfr(listByItem, function(df, dataset_id) {
    df %>%
      # dplyr::select(-Total) %>%
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
  
  ### Info for global assessment 
  resGlobalNormScore <- sapply(finalList, rownames,
                               simplify = F, USE.NAMES = T)
  
  
  ### Return
  resultado <- list(
    resByItem = resNormScore, 
    resGlobal = resGlobalNormScore
  )
  
  return(resultado)
}

resAllGS <- readRDS(file = "AssessmentFiles/results_Gold_Standard.rds")
resNSCorrected <- readRDS(file = "AssessmentFiles/normScore_CORRECTED_dataGS_All.rds")
resNS<- readRDS(file = "AssessmentFiles/normScore_dataGS_All.rds")

resAllNS_Corrected <- formatNormScoreResults(finalList = resNSCorrected)
resAllNS <- formatNormScoreResults(finalList = resNS)


# Function - Comparison by item ####

resultsByItem <- function(resNormScore, resGS){
  res <- resNormScore %>%
    dplyr::inner_join(resGS,
                      by = c("ID", "Item", "Normalization"))
  
  # Kendall Cor by item and dataset #
  kendallCor <- res %>%
    # slice(1:1200) %>%
    group_by(ID, Item) %>%
    summarise(
      ## ---- Kendall ---
      tau = cor(RankGS,
                Ranking,
                method = "kendall",
                use = "complete.obs"),
      ## ---- Spearman ---
      spearman = cor(RankGS,
                     Ranking,
                     method = "spearman",
                     use = "complete.obs"),
      ## ---- Pairwise accuracy ---
      pairAccuracy = pairwiseAccuracy(manual_rank = RankGS, 
                                      score_value = Ranking),
      .groups = "drop"
    ) %>% 
    as.data.frame()
  
  # Summarize results #
  # Table #
  tableByItem <- kendallCor %>%
    group_by(Item) %>%
    summarise(
      n_total = n(),
      n_valid = sum(!is.na(tau)),
      prop_valid = n_valid / n_total,
      
      # ---- Kendall ---
      tau_median = median(tau, na.rm = TRUE),
      tau_mean   = mean(tau, na.rm = TRUE),
      tau_sd     = sd(tau, na.rm = TRUE),
      prop_tau_gt_0_6 = mean(tau > 0.6, na.rm = TRUE),
      
      # ---- Spearman ---
      spearman_median = median(spearman, na.rm = TRUE),
      spearman_mean   = mean(spearman, na.rm = TRUE),
      spearman_sd     = sd(spearman, na.rm = TRUE),
      prop_spear_gt_0_7 = mean(spearman > 0.7, na.rm = TRUE),
      
      # ---- Pairwise accuracy ---
      acc_median = median(pairAccuracy, na.rm = TRUE),
      acc_mean   = mean(pairAccuracy, na.rm = TRUE),
      acc_sd     = sd(pairAccuracy, na.rm = TRUE),
      prop_acc_gt_0_8 = mean(pairAccuracy > 0.8, na.rm = TRUE),
      .groups = "drop"
    )
  
  resultados <- list(
    results = kendallCor,
    tabla = tableByItem
  )
  
  return(resultados)
}

summaryByItemCorrected <- resultsByItem(resNormScore = resAllNS_Corrected$resByItem, 
                                        resGS = resAllGS$byItem)
summaryByItemCorrected$tabla %>% View
summaryByItem <- resultsByItem(resNormScore = resAllNS$resByItem,
                               resGS = resAllGS$byItem)
summaryByItem$tabla %>% View
rbind(summaryByItemCorrected$tabla, summaryByItem$tabla) %>% View

# Da o mesmos resultado porque a versión correxida cos pesos só afecta ao global
summary(arsenal::comparedf(x = summaryByItem$tabla, y = summaryByItemCorrected$tabla))




# Function - Global comparison ####
resultsGlobal <- function(
    resGlobalNS, 
    resGlobalGS_BEST = NULL,
    resGlobalGS_BESTS = NULL, 
    resGlobalGS_WORST = NULL
    ){
  
  # Inicializando 
  resultados <- list()
  
  if (!is.null(resGlobalGS_BEST)){
    ## Hit top ###
    # Best (only one) => Does first position from GS and normScore match with 1, 2, 3 position?
    bestSet <- resGlobalGS_BEST[names(resGlobalNS)]
    resHitBestAll <- sapply(1:3, function(j) {sapply(names(resGlobalNS), function(i) 
      hitAtK(
        rankMethods = resGlobalNS[[i]], 
        targetSet = bestSet[[i]],   
        k = j, 
        type = "top"))
    })
    summaryHit(res = resHitBestAll)
    
    ## MRR ###
    # Best (only one) => Does first position from GS and normScore match with 1, 2, 3 position?
    bestSet <- resGlobalGS_BEST[names(resGlobalNS)]
    resMrrBestAll <- sapply(names(resGlobalNS), function(i) 
      mrrAtK(
        rankedMethods = resGlobalNS[[i]], 
        targetSet = bestSet[[i]],   
        type = "top"))
    
    resBest <- list(
      hitRes = resHitBestAll, 
      hitTable = summaryHit(res = resHitBestAll), 
      mrrRes = resMrrBestAll)
    
    resultados[["Best"]] <- resBest
  }
  
  if (!is.null(resGlobalGS_BESTS)){
    
    # HIT 
    # Best (set of normalizations) => Do some of the valid normalizations from GS and the selected by the normScore ( for the 1, 2, 3 position) match?
    bestSet <- resGlobalGS_BESTS[names(resGlobalNS)]
    resHitBestsAll <- sapply(1:3, function(j) {sapply(names(resGlobalNS), function(i) 
      hitAtK(
        rankMethods = resGlobalNS[[i]], 
        targetSet = bestSet[[i]],   
        k = j, 
        type = "top"))
    })
    
    # MRR
    # Best (set of normalizations) => Do some of the valid normalizations from GS and the selected by the normScore ( for the 1, 2, 3 position) match?
    bestSet <- resGlobalGS_BESTS[names(resGlobalNS)]
    resMrrBestsAll <- sapply(names(resGlobalNS), function(i) 
      mrrAtK(
        rankedMethods = resGlobalNS[[i]], 
        targetSet = bestSet[[i]],   
        type = "top"))
    table(resMrrBestsAll)
    
    resBests <- list(
      hitRes = resHitBestsAll, 
      hitTable = summaryHit(res = resHitBestsAll), 
      mrrRes = resMrrBestsAll)
    
    resultados[["Bests"]] <- resBests
  }
  
  
  if (!is.null(resGlobalGS_WORST)){
    
    ## Excluded normalizations ###
    worstToAsses <- intersect(names(resGlobalGS_WORST), names(resGlobalNS))
    worstSetGS <- resGlobalGS_WORST[worstToAsses]
    worstSetNS <- resGlobalNS[worstToAsses]
    resCovWorstAll <- sapply(names(worstSetGS), function(i) 
      bottomCoverage(
        rankedMethods = worstSetNS[[i]], 
        worstSet = worstSetGS[[i]], 
        # k = 5,
        prop = T))
    # table(resCovWorstAll)
    
    
    # # worst - control negativo, debería ser 0 todo
    # worstToAsses <- intersect(names(resGlobalGS_WORST), names(resGlobalNS))
    # worstSetGS <- resGlobalGS_WORST[worstToAsses]
    # worstSetNS <- resGlobalNS[worstToAsses]
    # resHitWorstAll <- sapply(1:3, function(j) {sapply(names(worstSetGS), function(i) 
    #   hitAtK(
    #     rankMethods = worstSetNS[[i]], 
    #     targetSet = worstSetGS[[i]],   
    #     k = j, 
    #     type = "top"))
    # })
    # resHitWorstAll
    # summaryHit(res = resHitWorstAll)
    
    # worst
    worstToAsses <- intersect(names(resGlobalGS_WORST), names(resGlobalNS))
    worstSetNS <- resGlobalGS_WORST[worstToAsses]
    worstSetGS <- resGlobalNS[worstToAsses]
    resMrrWorstAll <- sapply(names(worstSetGS), function(i) 
      mrrAtK(
        rankedMethods = worstSetNS[[i]], 
        targetSet = worstSetGS[[i]],   
        type = "bottom")
    )
    # table(resMrrWorstAll)
    
    
    resWorst <- list(
      coverageRes = resCovWorstAll, 
      mrrRes = resMrrWorstAll)
    
    resultados[["Worst"]] <- resWorst
  }
  
  return(resultados)
}

summaryGlobal <- resultsGlobal(resGlobalNS = resAllNS$resGlobal, 
                  resGlobalGS_BEST = resAllGS$globalBest,
                  resGlobalGS_BESTS = resAllGS$globalBests,
                  resGlobalGS_WORST = resAllGS$globalWorst)
summaryGlobal$Best$hitTable
summaryGlobal$Bests$hitTable
table(summaryGlobal$Worst$coverageRes)
table(summaryGlobal$Worst$mrrRes)

summaryGlobalCorrected <- resultsGlobal(resGlobalNS = resAllNS_Corrected$resGlobal, 
                  resGlobalGS_BEST = resAllGS$globalBest,
                  resGlobalGS_BESTS = resAllGS$globalBests,
                  resGlobalGS_WORST = resAllGS$globalWorst)
summaryGlobalCorrected$Best$hitTable
summaryGlobalCorrected$Bests$hitTable
table(summaryGlobalCorrected$Worst$coverageRes)
table(summaryGlobalCorrected$Worst$mrrRes)




