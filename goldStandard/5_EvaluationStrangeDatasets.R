###############################################################################-

############    Evaluation of strange datasets NS vs GS       #################-

###############################################################################-


# Julia G Curras - 2026/03/27
rm(list=ls())
graphics.off()
inputDir <- "C:/Users/julia/Documents/GitHub/normScore/goldStandard/AssessmentFiles/"
outDir <- "C:/Users/julia/Documents/GitHub/normScore/goldStandard/ProcessedDatasets/"
setwd("C:/Users/julia/Documents/GitHub/normScore/goldStandard")

source(file = "3_AssessmentFunctions.R")


#...........................................................................####
# Datasets ####
## Gold Standard results ####
resAllGS <- readRDS(file = paste0(inputDir, "results_Gold_Standard.rds")) # Manual ranking

## NormScore final results ####
# Item0x3 and Item2 scaledx0.1
resNS <- readRDS(file = paste0(inputDir, "normScore_dataGS_All_Item0_x3_Item2_01.rds")) # Originales
resAllNS <- formatNormScoreResults(finalList = resNS, item0 = T)



#...........................................................................####
# Results ####
## By item ####
summaryByItem <- resultsByItem(resNormScore = resAllNS$resByItem, 
                               resGS = resAllGS$byItem)
summaryByItem$results %>% View

# Selecting datasets with pairwise = 0 and cor = -1
summaryByItem$results %>% 
  filter(pairAccuracy == 0.2 | spearman < -0.5 | tau < -0.5) %>%
  arrange(Item)
# Revisamos estes datasets. 
# ID  Item       tau   spearman pairAccuracy
# 1   PXD041237 Item2 -1.000000 -1.0000000         0.00
# 2   PXD048564 Item2 -1.000000 -1.0000000         0.00
# 3 PXD054817.2 Item2 -1.000000 -1.0000000         0.00
# 4   PXD064948 Item2 -1.000000 -1.0000000         0.00
# 5   PXD065898 Item2 -1.000000 -1.0000000         0.00
# 6   PXD055964 Item3 -0.591608 -0.7457995         0.15
# 7   PXD014311 Item5 -0.439155 -0.5070926         0.20
# 8   PXD068597 Item5 -0.500000 -0.5773503         0.00




### Global ####
# Non corrected
summaryGlobal <- resultsGlobal(resGlobalNS = resAllNS$resGlobal, 
                               resGlobalGS_BEST = resAllGS$globalBest,
                               resGlobalGS_BESTS = resAllGS$globalBests,
                               resGlobalGS_WORST = resAllGS$globalWorst)
summaryGlobal$Best$hitTable
summaryGlobal$Bests$hitTable
table(summaryGlobal$Worst$coverageRes)
table(summaryGlobal$Worst$mrrRes)

# Corrected
summaryGlobalCorrected <- resultsGlobal(resGlobalNS = resAllNS_Corrected$resGlobal, 
                                        resGlobalGS_BEST = resAllGS$globalBest,
                                        resGlobalGS_BESTS = resAllGS$globalBests,
                                        resGlobalGS_WORST = resAllGS$globalWorst)
summaryGlobalCorrected$Best$hitTable
summaryGlobalCorrected$Bests$hitTable
table(summaryGlobalCorrected$Worst$coverageRes)
table(summaryGlobalCorrected$Worst$mrrRes)


summaryGlobal$Best$hitTable
summaryGlobalCorrected$Best$hitTable
