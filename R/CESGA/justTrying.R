# trying

library(dplyr)

cat("Hola mundo!")
iris <- iris %>% dplyr::select(Sepal.Length:Petal.Length)

saveRDS(iris, file = "trial.rds")