##############################################################################-

################    METANALYSIS WITH AMANIDE - BY SAMPLE  #######################-

##############################################################################-


# Julia G Currás - 2025/05/23
rm(list=ls())
graphics.off()
library(dplyr)
library(amanida)
# Basal data
basal <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", 
                           sheet = "1. Basal data", col_names = T)
basal <- basal[!(basal$infoPvalue == "No" & basal$infoAdjustedPvalue == "No"), ]
basal <- basal %>% filter(ID != "1098_977") 
table(basal$`Type of sample (grouped)`)
# Quantitative data (list of proteins by study)
dataOriginal <- readxl::read_xlsx(path = "../4_Processing/allProcessedProteinData.xlsx", 
                          sheet = "Quantitative (protein>1 study)", col_names = T)
# Diccionary 
diccionary <- dataOriginal %>%
  group_by(ProteinID) %>%
  reframe(
    Alias = {
      valores <- unique(na.omit(ProteinName))
      if (length(valores) == 0) as.character(ProteinID[1]) else if (length(valores) == 1) valores else paste(valores, collapse = "; ")
    },
    .groups = "drop"
  ) %>%
  dplyr::select(ProteinID, Alias)
diccionary <- as.data.frame(diccionary)
rownames(diccionary) <- diccionary$ProteinID
diccionary[which(diccionary$Alias == "-"), "Alias"] <- diccionary[which(diccionary$Alias == "-"), "ProteinID"]
sum(grepl(pattern = "NA;", x = diccionary$Alias, fixed = T)) # aliminar NA;
diccionary$Alias <- gsub(diccionary$Alias, pattern = "NA;", replacement = "", fixed = T)




# 1) BLOOD ####
## Processing data ####
# Selecting dataframes of interest #
basalQ <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Serum", "Plasma", "Saliva;Plasma")) 
df <- dataOriginal %>% filter(ID %in% c(basalQ$ID, "1027_15_Pla"))
table(df$ID)

# Removing proteins without logFc info #
sum(is.na(df$logFC))
df <- df %>% filter(!is.na(logFC))

# Finding duplicated proteins #
allFrec <- table(df$ProteinID)
table(allFrec) # 230 proteins únicas, 104 en al menos 2 estudios, 17 en al menos 3 estudios, 3 en 4 estudios
interestingFrec <- allFrec[which(allFrec > 1)]
df <- df %>% filter(ProteinID %in% names(interestingFrec)) # Keeping proteins that appears in 2 or more studies #

# Adding N info to protein data for omics-based metaanalysis #
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(Reference, TotalN)
data <- base::merge(df, basalQ_2, by = "Reference", all.x = T, all.y = T)

# Selecting relevant data for Amanida metanalysis #
sum(is.na(data$logFC))
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)

# Saving #
write.csv(x = df, file = "df_amanida_quantitativo_BLOOD.csv", 
          sep = ",", row.names = F, col.names = T, quote = F)


## AMANIDA META-ANALYSIS ####
# Open with Amanida #
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
datafile <- amanida_read(
  file = "df_amanida_quantitativo_BLOOD.csv", #_AllProts
  mode = "quan", coln, separator=",")

# Analyzing... #
amanida_result <- compute_amanida(datafile, comp.inf = F)

# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Displaying results #
volcano_plot(amanida_result, cutoff = c(0.05, 1.5))
vote_plot(amanida_result, counts = 2)
explore_plot(datafile, type = "sub", counts = 1)


## Saving Blood results ####
xlsx::write.xlsx(x = as.data.frame(amanida_result@stat), 
                 col.names = T, row.names = F, append = F,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "BLOOD - metaanalysis")
xlsx::write.xlsx(x = as.data.frame(amanida_result@vote), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "BLOOD - vote counting")


# EYE  ####
## Processing data ####
# Selecting dataframes of interest #
basalQ <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Eye")) 
df <- dataOriginal %>% filter(Reference %in% basalQ$Reference)
table(df$ID)

# Removing proteins without logFc info #
sum(is.na(df$logFC))
df <- df %>% filter(!is.na(logFC))

# Finding duplicated proteins #
allFrec <- table(df$ProteinID)
table(allFrec) # 1529 proteins únicas, 625 en al menos 2 estudios, 163 en al menos 3 estudios, 38 en 4 estudios
interestingFrec <- allFrec[which(allFrec > 1)]
df <- df %>% filter(ProteinID %in% names(interestingFrec))# Keeping proteins that appears in 2 or more studies #

# Adding N info #
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(Reference, TotalN)
data <- base::merge(df, basalQ_2, by = "Reference", all.x = T, all.y = T)

# Selecting relevant data for Amanida metanalysis #
sum(is.na(data$logFC))
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)

# Saving# 
write.csv(x = df, file = "df_amanida_quantitative_Samples.csv", 
          sep = ",", row.names = F, col.names = T, quote = F)

## AMANIDA META-ANALYSIS ####
# Open with Amanida #
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
datafile <- amanida_read(
  file = "df_amanida_quantitative_Samples.csv", #_AllProts
  mode = "quan", coln, separator=",")

# Analyzing... #
amanida_result <- compute_amanida(datafile, comp.inf = F)

# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Displaying results #
volcano_plot(amanida_result, cutoff = c(0.05, 2))
vote_plot(amanida_result, counts = 2)
explore_plot(datafile, type = "sub", counts = 1)


## Saving EYE results ####
xlsx::write.xlsx(x = as.data.frame(amanida_result@stat), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "EYE - metaanalysis")
xlsx::write.xlsx(x = as.data.frame(amanida_result@vote), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "EYE - vote counting")



# URINE ####
## Processing data ####
# Selecting dataframes of interest #
basalQ <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Urine")) 
df <- dataOriginal %>% filter(Reference %in% basalQ$Reference)
table(df$ID)

# Removing proteins without logFc info #
sum(is.na(df$logFC))
df <- df %>% filter(!is.na(logFC))

# Duplicated IDs across datasets #
allFrec <- table(df$ProteinID)
table(allFrec) # 1529 proteins únicas, 625 en al menos 2 estudios, 163 en al menos 3 estudios, 38 en 4 estudios
interestingFrec <- allFrec[which(allFrec > 1)]

# Adding N info #
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(Reference, TotalN)
data <- base::merge(df, basalQ_2, by = "Reference", all.x = T, all.y = T)

# Selecting relevant data for Amanida metanalysis #
sum(is.na(data$logFC))
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)

# Saving
write.csv(x = df, file = "df_amanida_quantitative_Samples.csv", 
          sep = ",", row.names = F, col.names = T, quote = F)


## AMANIDA META-ANALYSIS ####
# Open with Amanida #
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
datafile <- amanida_read(
  file = "df_amanida_quantitative_Samples.csv", #_AllProts
  mode = "quan", coln, separator=",")

# Analyzing... #
amanida_result <- compute_amanida(datafile, comp.inf = F)


# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Displaying results #
volcano_plot(amanida_result, cutoff = c(0.05, 4))
vote_plot(amanida_result, counts = 2)
explore_plot(datafile, type = "sub", counts = 1)

## Saving URINE results ####
xlsx::write.xlsx(x = as.data.frame(amanida_result@stat), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "URINE - metaanalysis")
xlsx::write.xlsx(x = as.data.frame(amanida_result@vote), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "URINE - vote counting")


# SALIVA ####
## Processing data ####
# Selecting dataframes of interest #
basalQ <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Saliva"))
df <- dataOriginal %>% filter(ID %in% c(basalQ$ID, "1027_15_Sal"))
table(df$ID)

# Removing proteins without logFc info #
sum(is.na(df$logFC))
df <- df %>% filter(!is.na(logFC))

# Finding duplicated proteins #
allFrec <- table(df$ProteinID)
table(allFrec) # 230 proteins únicas, 104 en al menos 2 estudios, 17 en al menos 3 estudios, 3 en 4 estudios

# Adding N info #
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(Reference, TotalN)
data <- base::merge(df, basalQ_2, by = "Reference", all.x = T, all.y = T)

# Selecting relevant data for Amanida metanalysis #
sum(is.na(data$logFC))
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)

# Saving
write.csv(x = df, file = "df_amanida_quantitative_Samples.csv", 
          sep = ",", row.names = F, col.names = T, quote = F)


## AMANIDA META-ANALYSIS ####
# Open with Amanida 
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
datafile <- amanida_read(
  file = "df_amanida_quantitative_Samples.csv", #_AllProts
  mode = "quan", coln, separator=",")

# Analyzing... #
amanida_result <- compute_amanida(datafile, comp.inf = F)

# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Displaying results #
volcano_plot(amanida_result, cutoff = c(0.05, 3))
vote_plot(amanida_result, counts = 2)
explore_plot(datafile, type = "sub", counts = 1)


## Saving SALIVA results ####
xlsx::write.xlsx(x = as.data.frame(amanida_result@stat), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "SALIVA - metaanalysis")
xlsx::write.xlsx(x = as.data.frame(amanida_result@vote), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "SALIVA - vote counting")





# EV ####
## Processing data ####
# Selecting dataframes of interest #
basalQ <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Extracellular vesicles (blood)")) 
df <- dataOriginal %>% filter(Reference %in% basalQ$Reference)
table(df$ID)

# Removing proteins without logFc info #
sum(is.na(df$logFC))
df <- df %>% filter(!is.na(logFC))

# Duplicated IDs across datasets #
allFrec <- table(df$ProteinID)
table(allFrec) # 1529 proteins únicas, 625 en al menos 2 estudios, 163 en al menos 3 estudios, 38 en 4 estudios

# Adding N info #
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(Reference, TotalN)
data <- base::merge(df, basalQ_2, by = "Reference", all.x = T, all.y = T)

# Selecting relevant data for Amanida metanalysis #
sum(is.na(data$logFC))
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)

# Saving
write.csv(x = df, file = "df_amanida_quantitative_Samples.csv", 
          sep = ",", row.names = F, col.names = T, quote = F)


## AMANIDA META-ANALYSIS ####
# Open with Amanida #
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
datafile <- amanida_read(
  file = "df_amanida_quantitative_Samples.csv", #_AllProts
  mode = "quan", coln, separator=",")

# Analyzing... #
amanida_result <- compute_amanida(datafile, comp.inf = F)

# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Displaying results #
volcano_plot(amanida_result, cutoff = c(0.05, 6))
vote_plot(amanida_result, counts = 2)
explore_plot(datafile, type = "sub", counts = 1)

## Saving EV results ####
xlsx::write.xlsx(x = as.data.frame(amanida_result@stat), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "EV - metaanalysis")
xlsx::write.xlsx(x = as.data.frame(amanida_result@vote), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "EV - vote counting")




# VAT ####
## Processing data ####
# Selecting dataframes of interest #
basalQ <- basal %>% 
  filter(`Type of sample (grouped)` %in% c("Visceral Adipose Tissue"))
df <- dataOriginal %>% filter(Reference %in% basalQ$Reference)
table(df$ID)

# Removing proteins without logFc info #
sum(is.na(df$logFC))
df <- df %>% filter(!is.na(logFC))

# Duplicated IDs across datasets #
allFrec <- table(df$ProteinID)
table(allFrec) # 284 proteins únicas, 6 en 2 estudios

# Adding N info #
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(Reference, TotalN)
data <- base::merge(df, basalQ_2, by = "Reference", all.x = T, all.y = T)

# Selecting relevant data for Amanida metanalysis #
sum(is.na(data$logFC))
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)

# Saving
write.csv(x = df, file = "df_amanida_quantitative_Samples.csv", 
          sep = ",", row.names = F, col.names = T, quote = F)


## AMANIDA META-ANALYSIS ####
# Open with Amanida #
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
datafile <- amanida_read(
  file = "df_amanida_quantitative_Samples.csv", 
  mode = "quan", coln, separator=",")

# Analyzing... #
amanida_result <- compute_amanida(datafile, comp.inf = F)

# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Displaying results #
volcano_plot(amanida_result, cutoff = c(0.05, 2))
vote_plot(amanida_result, counts = 2)
explore_plot(datafile, type = "sub", counts = 1)

## Saving quantitative Amanida results ####
xlsx::write.xlsx(x = as.data.frame(amanida_result@stat), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "VAT - metaanalysis")
xlsx::write.xlsx(x = as.data.frame(amanida_result@vote), 
                 col.names = T, row.names = F, append = T,
                 file = "resultsBySample_Quantitative.xlsx", 
                 sheetName = "VAT - vote counting")




