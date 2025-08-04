##############################################################################-

######   SELECTING INTERESTING PROTEINS FOR TRADITIONAL METANALYSIS    #######-

##############################################################################-


# Julia G Currás - 2025/05/15
rm(list=ls())
graphics.off()
library(dplyr)



# Initial data ####
## Selecting list of proteins from each study that provided raw data ####
data <- readxl::read_xlsx(path = "../4_Processing/allProcessedProteinData.xlsx", 
                          sheet = "Quantitative (23 studies)", col_names = T)
rawData <- c("16_128", "284_268", "377_100", "425_58", "432_287", "446_99", 
             "440_211", "511_252", "944_401", "231_253", "499_339", "1046_653") # raw data studies
df <- data %>%
  filter(ID %in% rawData)
table(df$ID)
df$logFC <- as.numeric(df$logFC)
df$pvalAdjFinal <- as.numeric(df$pvalAdjFinal)

## Joining basal information to protein list ####
basal <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", 
                           sheet = "1. Basal data", col_names = T)
# data <- basal %>% 
#   dplyr::select(ID, Reference, `Type of sample (grouped)`) %>%
#   filter(ID %in% rawData)
sampleInfo <- basal %>% dplyr::select(ID, Reference, `Type of sample (grouped)`)


# Loading interesting prots form discusion, tables, validation, etc from different articles ####
intProts <- readxl::read_excel(path = "3_RelevantBiomarkers.xlsx", 
                               sheet = 1, col_names = T)
table(table(intProts$ProteinID))
intProts <- merge(intProts, sampleInfo, by = "ID", all.x = T)





# Finding common proteins ####
## ALL SAMPLES ####
## For each proteins, the number of studies in which they appeared
allFrec <- table(df$ProteinID)
table(allFrec) # 26 shared between 8 studies; 8 shared between 9; 1 protein was identified in 10 studies


### Selection (I): identified in 10 studies ####
## P02763 was identified in 10 studies out o 12 with raw data
interestingFrec <- allFrec[which(allFrec > 9)] # 1
interestingFrec # P02763 
df %>% filter(ProteinID == names(interestingFrec)) %>% View
datos10All<- df %>% 
  filter(ProteinID %in% names(interestingFrec)) %>%
  as.data.frame()
xlsx::write.xlsx(x = datos10All, file = "3_proteinsForMetanalisis.xlsx",
                 col.names = T, row.names = F, showNA = F, 
                 sheetName = "10_All")


### Selection (II): identified in >7 studies & interesting protein ####
# From proteins identified in >7 studies, those that appear in the list of interesting
# proteins were selected
prot8_9 <- df %>% filter(ProteinID %in% intProts$ProteinID) %>% 
  group_by(ProteinID) %>% 
  arrange(length(ProteinID)) %>%
  summarise(n = length(ProteinID)) %>% 
  filter(n> 7) %>%
  ungroup() %>% pull(ProteinID)
df %>% filter(ProteinID %in% prot8_9) %>% View

## Saving 
datos8_9_All <- df %>% 
  filter(ProteinID %in% prot8_9) %>%
  arrange(ProteinID) %>% 
  as.data.frame()
xlsx::write.xlsx(x = datos8_9_All, file = "3_proteinsForMetanalisis.xlsx",
                 col.names = T, row.names = F, showNA = F, append = T,  
                 sheetName = "8_9_All")



## BY SAMPLE ####
### EYE ####
eyeData <- basal %>% filter(`Type of sample (grouped)` == "Eye") %>% pull(ID)
dfEye <- df %>% dplyr::filter(ID %in% eyeData)

## For each protein, number of studies 
allFrec <- table(dfEye$ProteinID)
table(allFrec)
allFrec[which(allFrec > 1)] # Proteins identified in more than one study

## Interesting proteins 
prots <- intProts %>% filter(`Type of sample (grouped)` == "Eye" ) %>% pull(ProteinID)
dfEye %>% filter(ProteinID %in% prots) %>% View 
# Lysozyme C P61626 (n=4) / P02649 APOE (3) / P25311 AZGP1 (3) 

## Saving 
datosEye <- dfEye %>% 
  filter(ProteinID %in% c("P61626", "P02649", "P25311")) %>%
  arrange(ProteinID) %>% as.data.frame()
xlsx::write.xlsx(x = datosEye, file = "3_proteinsForMetanalisis.xlsx",
                 col.names = T, row.names = F, showNA = F, append = T,  
                 sheetName = "4_3_Eye")



### BLOOD ####
bloodData <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Serum", "Plasma", "Extracellular vesicles (blood)")) %>%
  pull(ID)
dfBlood <- df %>% dplyr::filter(ID %in% bloodData)

## For each protein, number of studies 
allFrec <- table(dfBlood$ProteinID)
table(allFrec) 
allFrec[which(allFrec > 1)] # Proteins identified in more than one study

## Interesting proteins
prots <- intProts %>% filter(ID %in% unique(dfBlood$ID)) %>% pull(ProteinID)
dfBlood %>% filter(ProteinID %in% prots) %>% View 
protsBlood3 <- dfBlood %>% 
  filter(ProteinID %in% prots) %>% 
  group_by(ProteinID) %>% 
  arrange(length(ProteinID)) %>%
  summarise(n = length(ProteinID)) %>% 
  filter(n > 2) %>%
  ungroup() %>% pull(ProteinID)

## Saving 
datosBlood <- dfBlood %>% 
  filter(ProteinID %in% protsBlood3) %>%
  arrange(ProteinID) %>% as.data.frame()
xlsx::write.xlsx(x = datosBlood, file = "3_proteinsForMetanalisis.xlsx",
                 col.names = T, row.names = F, showNA = F, append = T, 
                 sheetName = "3_Blood")





## URINE ####
urineData <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Urine")) %>%
  pull(ID)
dfUrine <- df %>% dplyr::filter(ID %in% urineData)

## For each protein, number of studies 
allFrec <- table(dfUrine$ProteinID)
table(allFrec)
allFrec[which(allFrec > 1)] # # Proteins identified in more than one study

## Interesting proteins
prots <- intProts %>% filter(ID %in% unique(dfUrine$ID)) %>% pull(ProteinID)
dfUrine %>% filter(ProteinID %in% prots) %>% View 

## Saving 
datosUrine <- dfUrine %>% 
  filter(ProteinID %in% c("P01011", "P02768")) %>%
  arrange(ProteinID) %>% as.data.frame()
xlsx::write.xlsx(x = datosUrine, file = "3_proteinsForMetanalisis.xlsx",
                 col.names = T, row.names = F, showNA = F, append = T, 
                 sheetName = "2_Urine")



