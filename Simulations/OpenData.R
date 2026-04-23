

# Comparing simulations - all results show be identically

# Julia G Curras - 2026/01/28
rm(list=ls())
graphics.off()
setwd("C:/Users/julia/Documents/GitHub/normScore-development/Simulations")
library(dplyr)


# TRIAL DATASETS ####
# Open
archivos <- list.files(path = "TrialDatasets/", pattern = "item2_result_n")

resList <- list()
for (i in archivos){
  data <- readRDS(file = paste0("TrialDatasets/", i))
  resList[[i]] <- data$item2
}

df <- do.call(cbind, resList)



# Compare: 
resComp <- apply(df, 1, function(x) length(unique(x)) == 1)
table(resComp) # equal elements by row => OKKKKK



# FINAL DATASETS ####
# Open
archivos <- list.files(path = "FinalDatasets/", pattern = "20Rep.rds")

for (i in archivos){
  cat("Archivo: ", i, "\n")
  data <- readRDS(file = paste0("FinalDatasets/", i))
  data <- as.data.frame(data)
  cat("Filas: ", nrow(data), "\n")
  if (i == "item0156_result_n32_20Rep.rds"){
    df <- data
  } else{
    df <- merge(df, data, by =c("n_proteins", "n_per_group", "especificos", "semilla"), all.x= T)
  }
}


df <- df %>% 
  dplyr::select(n_proteins:n_per_group, semilla, especificos, Item0, Item1, 
                item2, item3, item4, Item5, Item6, normScore) %>%
  rename(
    Item2 = item2,
    Item3 = item3,
    Item4 = item4
  ) %>%
  dplyr::arrange(across(all_of(c("n_proteins", "n_per_group", "semilla"))), especificos)
saveRDS(df, file = "allSimulations.rds")






