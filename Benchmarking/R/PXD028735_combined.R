###############################################################################\
##################             All together             ###################\
###############################################################################\

# Set directory
setwd("C:/Users/julia/Documents/GitHub/normScore-development/Benchmarking/R")
rm(list=ls())
graphics.off()

# Libraries ####
library(dplyr)
library(ggplot2)


# Load data ####

ddaRES <- readRDS(
  file = "../Data/PXD028735/PXD028735_DDA_Complete_Analysis.rds"
)
dfResDDA <- ddaRES$allMetrics

diaRES <- readRDS(
  file = "../Data/PXD028735/PXD028735_DIA_Complete_Analysis.rds"
)
dfResDIA <- diaRES$allMetrics



# Combining results ####

dfAll <- bind_rows(
  dfResDIA |> mutate(Dataset = "DIA"),
  dfResDDA |> mutate(Dataset = "DDA")
)

dfRanks <- dfAll |>
  select(Dataset, Normalizations, starts_with("Rank"))  |>
  group_by(Dataset) |> 
  mutate(
    RankExternal = rowMeans(cbind(RankMAE, RankAUC))) |>
  ungroup()

## Global concordance DIA-DDA ####
cor(
  dfRanks$RankNS,
  dfRanks$RankExternal,
  method = "kendall"
)


## Concordance by criteria & external ####
dfCorAll <- cor(dfRanks[, 3:9], method = "kendall")
corrplot::corrplot.mixed(dfCorAll, upper = "ellipse")



# Average position by method ####
meanRanks <- dfRanks |>
  group_by(Normalizations) |>
  summarise(
    NormScore = mean(RankNS),
    External = mean(RankExternal),
    .groups = "drop"
  ) |>
  arrange(External)

meanRanks

cor(
  meanRanks$NormScore,
  meanRanks$External,
  method = "kendall"
)


# SAVING ####
lista <- list(
  allMetrics = dfAll, 
  dfCor = dfCorAll,
  meanDDADIA = meanRanks)


saveRDS(
  lista, 
  file = "../Data/PXD028735/PXD028735_combined.rds"
)

