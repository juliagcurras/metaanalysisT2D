######################################################-

############   Kappa index - peer review    ###########

######################################################-

# Julia G Currás - 26/03/2025

rm(list=ls())
graphics.off()
library(dplyr)


# Loading data ####
dfRPL <- readxl::read_xlsx(path = "peerReview.xlsx", 
                           sheet = "PerezLois", col_names = T)
dfJGC <- readxl::read_xlsx(path = "peerReview.xlsx", 
                           sheet = "GarciaCurras", col_names = T)
df <- merge(x = dfRPL, y = dfJGC, by = c("ID1", "ID2"), all = T)

# Select data ####
df <- df %>% dplyr::select(ID1, ID2, Decision_RPL, Decision_JGC, Title.x:DOI.x)
df$Decision_RPL <- ifelse(df$Decision_RPL %in% c("in", "In"), "In", 
                          ifelse(df$Decision_RPL %in% c("out", "Out"), "Out", "Error"))

# Kappa estimation ####
vcd::Kappa(table(df$Decision_RPL, df$Decision_JGC))

