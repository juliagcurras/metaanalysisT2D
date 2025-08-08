##############################################################################-

#################       METANALYSIS WITH AMANIDE     ########################-

##############################################################################-


# Julia G Currás - 2025/04/23
rm(list=ls())
graphics.off()
# setwd() # set working directory to the current directory
library(dplyr)
library(amanida)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(9396)



# Data ####
basal <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", 
                                sheet = "1. Basal data", col_names = T)


# Qualitative Analysis ####
data <- readxl::read_xlsx(path = "../4_Processing/allProcessedProteinData.xlsx", 
                        sheet = "Qualitative (protein>1 study)", col_names = T)

## Dictionary ####
diccionary <- data %>%
  group_by(ProteinID) %>%
  reframe(
    Alias = {
      valores <- unique(na.omit(ProteinName))
      if (length(valores) == 1) valores else paste(valores, collapse = "; ")
    },
    .groups = "drop"
  ) %>%
  dplyr::select(ProteinID, Alias)
table(duplicated(diccionary$ProteinID))
diccionary <- as.data.frame(diccionary)
rownames(diccionary) <- diccionary$ProteinID

## Preparing data ####
df <- data %>%
  filter(!is.na(logFC)) %>%
  dplyr::select(ProteinID, logFC, ID)
df$Behaviour <- ifelse(is.na(df$logFC), NA, 
                       ifelse(df$logFC < 0, "Down", 
                              ifelse(df$logFC > 0, "Up", NA)))
# Adding 1098_977 info (we only have the direction of the logFC but not the value)
df1098 <- data[which(data$ID == "1098_977"), c("ProteinID", "logFC", "ID")]
df1098$Behaviour <- "Up"
df1098[which(df1098$ProteinID == "P20700"), "Behaviour"] <- "Down"
table(df1098$Behaviour)
df <- rbind(df, df1098)
# Selecting rows
df <- df %>% 
  dplyr::select(ProteinID, Behaviour, ID) 
df[which(df$ID %in% c("1027_15_Pla", "1027_15_Sal")), "ID"] <- "1027_15"

# Change reference name
dfAux <- basal %>% dplyr::select(ID, Reference)
df <- merge(df, dfAux, by = "ID", all = T)[,-1]
colnames(df) <- c("Protein", "Change", "Reference")
# Gardalo para abrilo coa súa función
write.csv(x = df, file = "df_amanida_qualitativo.csv", 
          sep = ";", row.names = F, col.names = T, quote = F)

## Open with Amanida data ####
coln =  c("Protein", "Change", "Reference")
datafile <- amanida_read(
  file = "df_amanida_qualitativo.csv",
  mode = "qual", 
  coln = coln, 
  separator=",")

## Analyzing... ####
amanida_result <- amanida_vote(datafile)
amanida_result@vote %>% View
vote_plot(amanida_result, counts = 5)
# datafile$id_2 <- diccionary[datafile$id, "Alias"] # Non funciona moi ben, moitas proteínas sen nome
explore_plot(datafile, type = "sub", counts = 3)

## Saving qualitative Amanida results ####
output <- list(result = amanida_result, data = datafile)

### Vote counting results ####
# Adapted from Amanida #
counts = 3
cuts = counts
type = "sub"
trend = NULL
trend_l = NULL
N = NULL
vc = NULL
. = NULL
cont = NULL
lab = NULL
dataVote <- datafile
dt <- filter(mutate(summarise(dplyr::group_by(mutate(dplyr::group_by(
  mutate(dataVote, trend_l = case_when(trend == -1 ~ "Down-regulated", T ~ "Up-regulated")),id), 
  vc = sum(trend)), id, trend_l), cont = n(), total_N = sum(N), vc = unique(vc), lab = c("Vote-counting")), 
  cont = case_when(trend_l == "Down-regulated" ~ cont * -1, T ~ cont * 1)), vc >  cuts | vc < -1 * cuts)
# Structuration
data_wide <- dt %>%
  pivot_wider(names_from = trend_l, values_from = cont) %>%
  dplyr::select(-total_N, -lab, -vc)
dfFinal <- merge(data_wide, amanida_result@vote, all.X = T, by = "id")
dfFinal <- merge(dfFinal, diccionary, by.x = "id", by.y = "ProteinID", all.x = T)
dfFinal <- dfFinal %>% 
  dplyr::select(id, Alias, everything()) %>%
  dplyr::rename(
    `Protein ID` = id, 
    `Gene name` = Alias, 
    Votes = votes, 
    `N articles` = articles, 
    `Vote counting` = vote_counting
  ) %>%
  arrange(desc(Votes))
dfFinal[is.na(dfFinal)] <- 0
xlsx::write.xlsx(x = dfFinal, 
                 col.names = T, row.names = F, append = T,
                 file = "../5_1_QualitativeAssessment/results_Qualitative.xlsx", 
                 sheetName = "Vote counting")





#.-----------------------------------------------------------------------.####
# Quantitative ####
data <- readxl::read_xlsx(path = "../4_Processing/allProcessedProteinData.xlsx", 
                          sheet = "Quantitative (protein>1 study)", col_names = T)
basalQ <- basal[!(basal$infoPvalue == "No" & basal$infoAdjustedPvalue == "No"), ]
basalQ <- basalQ %>% filter(ID != "1098_977") 

## Dictionary ####
diccionary <- data %>%
  group_by(ProteinID) %>%
  reframe(
    Alias = {
      valores <- unique(na.omit(ProteinName))
      if (length(valores) == 0) as.character(ProteinID[1]) else if (length(valores) == 1) valores else paste(valores, collapse = "; ")
    },
    .groups = "drop"
  ) %>%
  dplyr::select(ProteinID, Alias)
table(duplicated(diccionary$ProteinID))
diccionary <- as.data.frame(diccionary)
rownames(diccionary) <- diccionary$ProteinID
diccionary[which(diccionary$Alias == "-"), "Alias"] <- diccionary[which(diccionary$Alias == "-"), "ProteinID"]
sum(grepl(pattern = "NA;", x = diccionary$Alias, fixed = T)) # aliminar NA;
diccionary$Alias <- gsub(diccionary$Alias, pattern = "NA;", replacement = "", fixed = T)

## Preparing data ####
data$ID <- sapply(data$ID, function(ids){
  ifelse(nchar(ids) > 8, stringr::str_sub(ids, end = -5) , ids)
})

# Adding N info
basalQ$TotalN <- basalQ$SampleSizeControls + basalQ$SampleSizeDiabetics
basalQ_2 <- basalQ %>% dplyr::select(ID, TotalN)
data <- base::merge(data, basalQ_2, by = "ID", all.x = T, all.y = T)

# Selecting relevant data
sum(is.na(data$logFC))
data %>% filter(is.na(logFC)) %>% View
df <- data %>%
  filter(!is.na(logFC)) %>% 
  dplyr::select(ProteinID, pvalAdjFinal, logFC, TotalN, Reference)
df$logFC <- as.numeric(df$logFC)
df$logFC <- 2^df$logFC
coln = c("Protein", "AdjPval", "logFC", "N total", "References")
colnames(df) <- coln
df <- df %>% filter(Protein != "Q6ZMP0") # aggregated logFC= 30, OUTLIERS, REMOVING IT SO 2734 PROTEINS

df$References <- gsub(df$References, pattern = ", ", replacement = " ", fixed = T)
# Saving
write.csv(x = df, file = "df_amanida_quantitativo.csv", #_AllProts (df para todas as proteinas sin quedarnos solo coas de al menos en 2 estudios)
          sep = ",", row.names = F, quote = F)

## Open with Amanida ####
datafile <- amanida_read(
  file = "df_amanida_quantitativo.csv", #_AllProts
  mode = "quan", coln, separator=",")

## Analyzing... ####
amanida_result <- compute_amanida(datafile, comp.inf = F)

## Saving quantitative Amanida results ####
### Aggregated values ####
res <- merge(amanida_result@stat , diccionary, by.x = "id", by.y = "ProteinID", all = T)
res$trend <- factor(res$trend, levels = c(-1, 1), labels = c("Downregulated", "Upregulated"))
res$N <- sapply(strsplit(res$reference, ";"), length)
res <- res %>% 
  arrange(desc(N), pval, abs(fc)) %>%
  dplyr::select(id, Alias, everything())
colnames(res) <- c("Protein ID", "Gene name", "Trend", 
                   "Aggregated adjusted p-value", "Fold change (FC)", 
                   "No. subjects", "References", "No. references")
xlsx::write.xlsx(x = as.data.frame(res), col.names = T, row.names = F, append = F,
                file = "results_Quantitative.xlsx", sheetName = "Metaanalysis")



## Displaying plots ####
amanida_result_Ori <- amanida_result
# Changing protein ids by protein names
amanida_result@stat$id <- diccionary[amanida_result@stat$id, "Alias"]
amanida_result@vote$id <- diccionary[amanida_result@vote$id, "Alias"]

# Plots
volcano_plot(amanida_result, cutoff = c(0.05, 4))

vote_plot(amanida_result, counts = 5)

# explore_plot(datafile, type = "sub", counts = 4)



## Assessment #####
amanida_result <- amanida_result_Ori
### Contribution analysis ####
dfRes <- amanida_result@stat 
tabAbs <- table(unlist(strsplit(dfRes$reference, split = "; ", fixed = T)))
tabRel <- tabAbs/length(unique(data$ProteinID))

tabOut <- data.frame(Reference = names(tabRel), 
                     N = as.vector(unname(tabAbs)), 
                     Percentage = as.vector(unname(tabRel)))
#### Saving quantitative Amanida results ####
xlsx::write.xlsx(x = tabOut, col.names = T, row.names = F, append = T,
                 file = "results_Quantitative.xlsx", sheetName = "Contribution analysis")


### Influence analysis ####
pathToFile <- paste0(getwd(), "/df_amanida_quantitativo.csv")

influenceAnalysis <- function(pathToFile){
  dfOriginal <- utils::read.table(file = pathToFile, header = T, sep = ",", 
                                  row.names = NULL, check.names = F)
  coln <- colnames(dfOriginal)
  datafile <- amanida_read(file = pathToFile, coln = coln, separator=",", mode = "quan")
  dfRes <- compute_amanida(datafile, comp.inf = F)@stat 
  
  vMAPE <- list()
  for (i in unique(dfOriginal$References)){
    # i <- unique(df$References)[1]
    # Exclude article by article
    df <- dfOriginal %>% dplyr::filter(References != i) 
    write.csv(x = df, file = "df_amanida_quantitativo_temporal.csv", #_AllProts (df para todas as proteinas sin quedarnos solo coas de al menos en 2 estudios)
              sep = ";", row.names = F, quote = F)
    datafile <- amanida_read(
      file = "df_amanida_quantitativo_temporal.csv", #_AllProts
      mode = "quan", coln, separator=",")
    ## Analyzing... 
    amanida_result_temp <- compute_amanida(datafile, comp.inf = F)
    dfResTemp <- amanida_result_temp@stat
    colnames(dfResTemp)[-1] <- paste0(colnames(dfResTemp)[-1], "_TEMP") 
    
    dfMerge <- merge(dfRes, dfResTemp, by = "id", all = T) 
    dfMerge <- dfMerge %>% dplyr::select(id, fc, fc_TEMP) %>% na.omit
    
    mapeTemp <- mean(abs((dfMerge$fc - dfMerge$fc_TEMP) / dfMerge$fc)) * 100
    vMAPE[[i]] <- mapeTemp
  }
  
  protsXarticle <- table(dfOriginal$References)
  dfMAPE <- data.frame(
    Datasets = names(vMAPE),
    MAPE = unlist(vMAPE, use.names = FALSE), 
    NProteins = as.vector(unname(protsXarticle[names(vMAPE)]))
  )
  
  # Plot 
  # Sort by MAPE magnitude
  dfMAPE$Datasets <- factor(dfMAPE$Datasets,
                          levels = dfMAPE$Datasets[order(dfMAPE$MAPE, decreasing = T)])
  
  pMAPE <- Biostatech::plotBarUnivar(var = dfMAPE$Datasets, 
                                     freqVar = dfMAPE$MAPE, 
                                     tituloY = "MAPE %", 
                                     sizeLetra = 20, 
                                     anguloX = 65)$grafico
  
  return(list(table = dfMAPE, graphMAPE = pMAPE))
}
salida <- influenceAnalysis(pathToFile = pathToFile)

#### Saving quantitative Amanida results ####
xlsx::write.xlsx(x = salida$table, col.names = T, row.names = F, append = T,
                 file = "results_Quantitative.xlsx", sheetName = "Influence analysis")





## Functional analysis ####
### GSEA Analysis ####
dfRes <- amanida_result@stat 
table(table(dfRes$id)>1)
dfRes$logFC <- log(dfRes$fc, base = 2)

#### Nominal p-value (Pval+logFC): Prepare list for analyzed ####
geneList <- dfRes %>%
  filter(!is.na(id), !is.na(logFC)) %>%
  mutate(ranking_metric = -log10(pval)*sign(logFC)) %>%
  group_by(id) %>%
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%
  arrange(-ranking_metric) %>% # sort descending (important!)
  tibble::deframe() # convert to named vector
geneList <- geneList[!is.infinite(geneList)]

##### GO ####
cgsea_res <- gseGO(geneList = geneList, 
                   ont = "ALL",
                   keyType = "UNIPROT",
                   OrgDb = "org.Hs.eg.db", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   nPermSimple = 1000, 
                   pvalueCutoff = 0.5,
                   pAdjustMethod = "BH",
                   seed = TRUE)
cgsea_res@result %>% View
cgsea_res_original <- cgsea_res

cgeCC <- cgsea_res
cgeCC@result <- cgeCC@result %>% filter(ONTOLOGY == "CC")
ridgeplot(cgeCC, fill = "pvalue",) + labs(x = "enrichment distribution")
gseaplot(cgeCC, by = "all", title = cgeCC$Description[1], geneSetID = 1)

cgeMF <- cgsea_res
cgeMF@result <- cgeMF@result %>% filter(ONTOLOGY == "MF")
ridgeplot(cgeMF, fill = "pvalue") + labs(x = "enrichment distribution")
gseaplot(cgeMF, by = "all", title = cgeMF$Description[1], geneSetID = 1)

cgeBP <- cgsea_res
cgeBP@result <- cgeBP@result %>% filter(ONTOLOGY == "BP")
ridgeplot(cgeBP, fill = "pvalue", showCategory = 15) + labs(x = "enrichment distribution")
gseaplot(cgeBP, by = "all", title = cgeBP$Description[2], geneSetID = 2)

xlsx::write.xlsx(x = as.data.frame(cgsea_res@result), 
                 col.names = T, row.names = F, append = T,
                 file = "results_Quantitative.xlsx", 
                 sheetName = "GSEA GO (nominal p-value)")


#### logFC: Prepare list for analyzed ####
geneList <- dfRes %>%
  filter(!is.na(id), !is.na(logFC)) %>%
  mutate(ranking_metric = logFC) %>%
  group_by(id) %>%
  summarise(ranking_metric = mean(logFC, na.rm = TRUE)) %>%
  arrange(-ranking_metric) %>% # sort descending (important!)
  tibble::deframe() # convert to named vector
geneList <- geneList[!is.infinite(geneList)]

##### GO ####
cgsea_res <- gseGO(geneList = geneList, 
                   ont = "ALL",
                   keyType = "UNIPROT",
                   OrgDb = "org.Hs.eg.db", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   nPermSimple = 1000, 
                   pvalueCutoff = 0.5,
                   pAdjustMethod = "BH",
                   seed = TRUE)
cgsea_res@result %>% View
cgsea_res_original <- cgsea_res

cgeCC <- cgsea_res
cgeCC@result <- cgeCC@result %>% filter(ONTOLOGY == "CC")
ridgeplot(cgeCC, fill = "pvalue", showCategory = 15) + labs(x = "enrichment distribution")
gseaplot(cgeCC, by = "all", title = cgeCC$Description[1], geneSetID = 1)

cgeMF <- cgsea_res
cgeMF@result <- cgeMF@result %>% filter(ONTOLOGY == "MF")
ridgeplot(cgeMF, fill = "pvalue", showCategory = 10) + labs(x = "enrichment distribution")
gseaplot(cgeMF, by = "all", title = cgeMF$Description[1], geneSetID = 1)

cgeBP <- cgsea_res
cgeBP@result <- cgeBP@result %>% filter(ONTOLOGY == "BP")
ridgeplot(cgeBP, fill = "pvalue", showCategory = 15) + labs(x = "enrichment distribution")
gseaplot(cgeBP, by = "all", title = cgeBP$Description[2], geneSetID = 2)

xlsx::write.xlsx(x = as.data.frame(cgsea_res@result), 
                 col.names = T, row.names = F, append = T,
                 file = "results_Quantitative.xlsx", 
                 sheetName = "GSEA GO (logFC)")



## Exploring important proteins ####
dfRes <- amanida_result@stat 
table(table(dfRes$id)>1)
dfRes$logFC <- log(dfRes$fc, base = 2)
dfRes$NTotal <- sapply(strsplit(dfRes$reference, split = ";", fixed = T), length)
dfRes %>% filter(pval < 0.05 & abs(logFC) > 3) %>% View
dfRes %>% filter(pval < 0.05 & NTotal >= 3 & abs(logFC) > 1) %>% View

### I) Important proteins #####
prots <- dfRes %>% filter(pval < 0.05 & abs(logFC) > 1.5) %>% pull(id) # Original

#### Interaction analysis - STRING ####
# Map #
targetID <- rbioapi::rba_string_map_ids(prots, species = 9606)
# Net # 
rbioapi::rba_string_network_image(ids = targetID$stringId, 
                                  required_score = 400,
                                  image_format = "image",
                                  node_labels_font_size = 20, 
                                  save_image = file.path("STRING_quantitative_filterProts_1.png"),
                                  species = 9606) # only query proteins

int_net <- rbioapi::rba_string_interactions_network(ids = targetID$stringId, # todas
                                           species = 9606,
                                           required_score = 400
                                           ) # medium confident
int_test <- rbioapi::rba_string_enrichment_ppi(
  ids = targetID$stringId, species = 9606, 
  required_score = 400,
  )

prots <- dfRes %>% filter(pval < 0.05 & abs(logFC) > 1.5) %>% 
  dplyr::select(-trend, -fc, -N_total)
dfInfo <- merge(prots, targetID, by.x = "id", by.y = "queryItem", all = T)
dfInfo <- dfInfo %>% 
  dplyr::select(id, preferredName, pval, logFC, annotation, NTotal, reference)
colnames(dfInfo) <- c("UniProt ID", "Protein name", "p-value", 
                      "logFC", "Description", "Nº references", "References")

# output <- list(tableNet = int_net, tableTest = int_test, tableInfo = dfInfo)
xlsx::write.xlsx(x = dfInfo, 
                 col.names = T, row.names = F, append = T,
                 file = "results_Quantitative.xlsx", 
                 sheetName = "Filtered proteins (I)")



### II) Important proteins #####
prots <- dfRes %>% filter(pval < 0.05 & NTotal >= 3 & abs(logFC) > 1) %>% pull(id) # Alternative

#### Interaction analysis - STRING ####
# Map #
targetID <- rbioapi::rba_string_map_ids(prots, species = 9606)
# Net # 
rbioapi::rba_string_network_image(ids = targetID$stringId, 
                                  required_score = 500,
                                  image_format = "image",
                                  node_labels_font_size = 20, 
                                  save_image = file.path("STRING_quantitative_filterProts_2.png"),
                                  species = 9606) # only query proteins

int_net <- rbioapi::rba_string_interactions_network(ids = targetID$stringId, # todas
                                                    species = 9606,
                                                    required_score = 500
) # medium confident
int_test <- rbioapi::rba_string_enrichment_ppi(
  ids = targetID$stringId, species = 9606, 
  required_score = 500 
)

prots <- dfRes %>% filter(pval < 0.05 & NTotal >= 3 & abs(logFC) > 1) %>% 
  dplyr::select(-trend, -fc, -N_total)
dfInfo <- merge(prots, targetID, by.x = "id", by.y = "queryItem", all = T)
dfInfo <- dfInfo %>% 
  dplyr::select(id, preferredName, pval, logFC, annotation, NTotal, reference) %>%
  arrange(desc(NTotal))
colnames(dfInfo) <- c("UniProt ID", "Protein name", "p-value", 
                      "logFC", "Description", "Nº references", "References")

xlsx::write.xlsx(x = dfInfo, 
                 col.names = T, row.names = F, append = T,
                 file = "results_Quantitative.xlsx", 
                 sheetName = "Filtered proteins (II)")





