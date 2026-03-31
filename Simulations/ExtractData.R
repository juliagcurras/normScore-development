##############################################################################-

#####################      CESGA MQ Computations    ###########################

##############################################################################-

# Julia G Curras - 25/10/2024


# Initial setup #
setwd("C:/Users/julia/Documents/GitHub/normScore/Simulations/CESGA")
library(dplyr)


# Global variables ####
carpetas <- "OutputFiles"
fios <- c(1, 4, 8, 16, 24, 32, 48, 64)

# New
data <- readxl::read_xlsx(path = "Parametros.xlsx", sheet = 1, col_names = T)
dfFinal <- data %>% filter(Proteins == "All")
dfFinal <- as.data.frame(dfFinal)
dfFinal$TotalThreads <- dfFinal$Cores*dfFinal$Threads*dfFinal$Nodes

# New variables ####
## Speed up ####
dfFinal$SpeedUp <- dfFinal[which(dfFinal$TotalThreads == 1), "Time (s)"]/dfFinal[, "Time (s)"]

## Efficiency ####
dfFinal$Efficiency <- dfFinal$SpeedUp/dfFinal$TotalThreads

## Memory per thread ####
dfFinal$`Memory Utilized (GB)` <- ifelse(dfFinal$Units...11 == "MB", dfFinal$`Memory Utilized`/1024, dfFinal$`Memory Utilized`)
dfFinal$`Memory By Thread (GB)` <- dfFinal$`Memory Utilized (GB)`/dfFinal$TotalThreads
dfFinal$`Memory By Thread (GB)`


# Split by type ####
dfPrin <- dfFinal[1:9, ]
dfSec <- dfFinal[c(1, 10:12, 1, 13:15), ]
dfSec$ThreadGroup <- rep(c(16, 24), each = 4)


# Saving df ####
output <- list(dfPrin, dfSec)
saveRDS(object = output, file = "./OutputFiles/data.rds")


# Graphical representations ####
## Time ####
Biostatech::plotScatter(base = dfPrin, 
                        varY = "Time (min)", 
                        tituloY = "Time (minutes)",
                        varX = "TotalThreads",
                        adjustLine = T, 
                        sizeDots = 4,
                        adjustType = "auto")$grafico
Biostatech::plotScatter(base = dfSec, 
                        varY = "Time (min)", 
                        tituloY = "Time (minutes)",
                        varX = "Cores",
                        varGrupo = "Group",
                        adjustLine = T, 
                        sizeDots = 4,
                        adjustType = "auto")$grafico

## Memory by thread ####
Biostatech::plotScatter(base = dfPrin, 
                        varY = "Memory By Thread (GB)", 
                        varX = "TotalThreads", 
                        adjustLine = F)$grafico


## SpeedUp ####
Biostatech::plotScatter(base = dfPrin, 
                        varY = "SpeedUp", 
                        varX = "TotalThreads", 
                        adjustLine = T, 
                        sizeDots = 4,
                        adjustType = "auto")$grafico
Biostatech::plotScatter(base = dfSec, 
                        varY = "SpeedUp", 
                        varX = "Cores", 
                        varGrupo = "Group",
                        adjustLine = T, 
                        sizeDots = 4,
                        adjustType = "auto")$grafico


## Efficiency ####
Biostatech::plotScatter(base = dfSec, 
                        varY = "Efficiency", 
                        varX = "Cores", 
                        varGrupo = "Group", 
                        adjustLine = T, 
                        sizeDots = 4,
                        adjustType = "auto")$grafico

Biostatech::plotScatter(base = dfPrin, 
                        varY = "Efficiency", 
                        varX = "TotalThreads", 
                        adjustLine = F)$grafico






