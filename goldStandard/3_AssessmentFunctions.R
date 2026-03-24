###############################################################################-

############                Assessment by items              #################-

###############################################################################-


# Julia G Curras - 2026/01/27

#...........................................................................####
# Set up ####
# Libraries
library(dplyr)
library(ggplot2)
library(boot)



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
# Function - Format normScore results ####


formatNormScoreResults <- function(finalList, item0 = T){
  
  ### Info for item assessment 
  # change wid format to longer format for comparison by item
  resNormScore <- purrr::imap_dfr(finalList, function(df, dataset_id) {
    df %>%
      dplyr::select(-Total, -TotalCorrected) %>%
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
  if (!item0){ # item0 = F si queremos saber os resultados sin corrección de item0
    finalList <- sapply(finalList, function(df){
      df <- df %>% dplyr::arrange(Total)
      df
    }, simplify = F, USE.NAMES = T)
  }
  # item0 = T e xa recuperamos os resultados ordenados directamente no normScore por TotalCorrected
  resGlobalNormScore <- sapply(finalList, rownames,
                               simplify = F, USE.NAMES = T)
  
  
  ### Return
  resultado <- list(
    resByItem = resNormScore, 
    resGlobal = resGlobalNormScore
  )
  
  return(resultado)
}


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