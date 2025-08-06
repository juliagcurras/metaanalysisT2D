##############################################################################-

#################   PREPARING DATA FOR METANALYSIS    ########################-

##############################################################################-


# Julia G Currás - 2025/05/15
rm(list=ls())
graphics.off()
library(dplyr)

# Loading basal and protein data bases ####
data <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", sheet = "1. Basal data")
dfProts <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", sheet = "2. Protein data")


# Processing of Proteomic Data ####
## Get Uniprot ID ####
# Select Protein names from proteins without UniprotID
genes <- dfProts %>% filter(is.na(ProteinID)) %>% pull(ProteinName)
dfGenes <- as.data.frame(genes) 
write.table(dfGenes, file = "UniProtMapping/UNIPROT_GenesToMap_2025_05_15.txt", quote = F, row.names = F)
#.----------------------    MAPPING AND SEARCH INTO UNIPROT    ---------------.#
proteinID <- readxl::read_xlsx(path = "UniProtMapping/UNIPROT_idmapping_GenesToMap_2025_05_15.xlsx", 
                               sheet = 1, col_names = T)
proteinID <- as.data.frame(proteinID)
# Only reviewed entries
proteinID <- proteinID %>% filter(Reviewed == "reviewed")
# Dealing with >1 Uniprot ID for the same protein/gene name
dfGenes$Protein <- unname(unlist(sapply(dfGenes$genes, function(i){
  if (i %in% proteinID$From){
    resultado <- proteinID[which(proteinID$From == i), "Entry"]
    if (length(resultado) > 1){
      return(paste0(resultado, collapse = ";"))
    } else{
      return(resultado)
    }
  } else {
    return(NA)
  }
}, simplify = F)))
table(dfProts[which(is.na(dfProts$ProteinID)), "ProteinName"] == genes) # genes match lack-of-protein ID rows
dfProts[which(is.na(dfProts$ProteinID)), "ProteinID"] <- dfGenes$Protein
table(is.na(dfProts$ProteinID)) # 581 proteins without Uniprot ID => REMOVING
dfProts <- dfProts %>% filter(!is.na(ProteinID))


## Get Gene/protein name ####
proteins <- dfProts %>% filter(is.na(ProteinName)) %>% pull(ProteinID)
proteinsDF <- as.data.frame(proteins) 
write.table(proteinsDF, file = "UniProtMapping/UNIPROT_ProteinsToMap_2025_05_15.txt", quote = F, row.names = F)
#.----------------------    MAPPING AND SEARCH INTO UNIPROT    ---------------.#
genesID <- readxl::read_xlsx(path = "UniProtMapping/UNIPROT_idmapping_ProteinsToMap_2025_05_15.xlsx", 
                               sheet = 1, col_names = T)
genesID <- as.data.frame(genesID)
colnames(genesID) <- c("Uniprot", "GeneName")
# Dealing with >1 gene/protein name for the same UniProt ID
proteinsDF$Genes <- unname(unlist(sapply(proteinsDF$proteins, function(i){
  if (i %in% genesID$Uniprot){
    resultado <- genesID[which(genesID$Uniprot == i), "GeneName"]
    if (length(resultado) > 1){
      return(paste0(resultado, collapse = ";"))
    } else{
      return(resultado)
    }
  } else {
    return(NA)
  }
}, simplify = F)))
table(dfProts[which(is.na(dfProts$ProteinName)), "ProteinID"] == proteins) 
dfProts[which(is.na(dfProts$ProteinName)), "ProteinName"] <- proteinsDF$Genes
table(is.na(dfProts$ProteinName)) # 4053 proteins without gene/protein name



## Any duplicated UniProt ID within a specific dataset? ####
for (i in unique(dfProts$ID)){ 
  cat(i, ": \t")
  cat(any(duplicated(dfProts[which(dfProts$ID == i), "ProteinID"]))) # 1097_685, 1276_2
  cat(" with total of ")
  cat(sum(duplicated(dfProts[which(dfProts$ID == i), "ProteinID"]))) 
  cat("\n")
}

### 1097_685: Aggregation ####
df1097 <- dfProts %>% filter(ID == "1097_685")
df1097$ProteinID[duplicated(df1097$ProteinID)] # Duplicated proteins
dfNew1097 <- df1097 %>%
  group_by(ProteinID, ProteinName, ProteinDescription, Reference, ID) %>% # Aggregate by protein ID
  summarise( # veage of logFC and p-value
    logFC = mean(logFC),
    pvalAdj = mean(pvalAdj),
    .groups = "drop"  
  )
dfNew1097$FC <- NA
dfNew1097$pval <- NA
dfNew1097 <- dfNew1097 %>% dplyr::select(ID, ProteinID:ProteinDescription, logFC, FC, pval, pvalAdj, Reference)
# Join to main database
dfProts <- dfProts %>% filter(ID != "1097_685")
dfProts <- rbind(dfProts, dfNew1097)


### 1276_2: Aggregation ####
dfAux <- dfProts %>% filter(ID == "1276_2")
proteinsToChange <- dfAux$ProteinID[duplicated(dfAux$ProteinID)] # Duplicated proteins
dfAux %>% filter(ProteinID %in% proteinsToChange) %>% View
# Name homogenization for FAM21
dfAux[which(dfAux$ProteinName %in% c("FAM21A", "FAM21B")), "ProteinName"] <- "WASHC2A" 
dfAux[which(dfAux$ProteinName == "WASHC2A"), "ProteinDescription"] <- "family with sequence similarity 21 (	WASHC2A, FAM21A, FAM21B)"
# Average for duplicated ids
dfNew1276 <- dfAux %>%
  group_by(ProteinID, ProteinName, ProteinDescription, Reference, ID) %>%
  summarise(
    logFC = mean(logFC),
    pval = mean(pval),
    .groups = "drop"  
  )
dfNew1276$FC <- NA
dfNew1276$pvalAdj <- NA
dfNew1276 <- dfNew1276 %>% dplyr::select(ID, ProteinID:ProteinDescription, logFC, FC, pval, pvalAdj, Reference)
# Joint everything
dfProts <- dfProts %>% filter(ID != "1276_2")
dfProts <- rbind(dfProts, dfNew1276)
nrow(dfProts)
# Chagne from tibble to dataframe
dfProts <- as.data.frame(dfProts)



## Homogenization ####
### Pvalues ####
adjPvals <- c()
for (i in unique(dfProts$ID)){
  # Adjustment for slightly differences in study IDs between basal data (data) and protein data (dfProts)
  if (nchar(i) > 8) {
    i_data <- substr(i, 1, nchar(i) - 4) 
  } else {
    i_data <- i
  }
  cat("\n ", i_data, "\n")
  # Checking p-value info 
  infoPval <- unname(unlist(as.vector(data[data$ID == i_data, c("infoPvalue", "infoAdjustedPvalue", "TotalProteins")])))
  if (all(infoPval[1] == "No", infoPval[2] == "Yes")){ # Only adjusted pvalues 
    adjPvals <- c(adjPvals, dfProts[dfProts$ID == i, "pvalAdj"])
  } else if (all(infoPval[1] == "No", infoPval[2] == "No")){ # Any information
    adjPvals <- c(adjPvals, rep(NA, nrow(dfProts[dfProts$ID == i, ])))
  } else if (infoPval[1] == "Yes"){ # P-values info, then adjusted p-values are estimated
    if (is.na(infoPval[3])){
      adjPvals <- c(adjPvals, rep(NA, nrow(dfProts[dfProts$ID == i, ])))
    } else {
      pvalores <- as.vector(as.numeric(dfProts[dfProts$ID == i, "pval"]))
      ni <- as.numeric(infoPval[3])
      ajustados <- stats::p.adjust(
        p = pvalores, 
        method = "BH",
        n = ni # Total proteins
        )
      adjPvals <- c(adjPvals, ajustados)
    }
  } else {
    cat("error")
  }
}

# Using new adjusted p-values when no adjusted p-value is provided, but the 
# original adjusted p-values are used always instead the new ones if this 
# data was provided
dfProts$pvalAdjNew <- adjPvals
dfProts$ProbaPVAL <- ifelse(is.na(dfProts$pval) & is.na(dfProts$pvalAdj), NA,
                            ifelse(!is.na(dfProts$pval) & is.na(dfProts$pvalAdj), dfProts$pvalAdjNew,
                                   ifelse(is.na(dfProts$pval) & !is.na(dfProts$pvalAdj), dfProts$pvalAdj,
                                          ifelse(!is.na(dfProts$pval) & !is.na(dfProts$pvalAdj), dfProts$pvalAdj, "Error"
                                                 )
                                          )
                                   )
                            )
# Checking that the ifelse have worked properly
dfProts %>% dplyr::select(pval, pvalAdj, pvalAdjNew, ProbaPVAL, ID) %>% View
# Checking the differences between new adjusted p-values and provided adjusted pvalues
completePvalInfo <- data %>% filter(infoPvalue == "Yes" & infoAdjustedPvalue == "Yes") %>% pull(ID)
dfProts %>% 
  filter(ID %in% completePvalInfo) %>%
  dplyr::select(pval, pvalAdj, pvalAdjNew, ID) %>% View

# Ok, keeping info of adjustment
dfProts$pvalAdjFinal <- dfProts$ProbaPVAL
dfProts$ProbaPVAL<- NULL
dfProts$pvalAdjNew <- NULL


### logFC ####
dfProts$FC <- as.numeric(dfProts$FC)
dfProts$logfcNew <- log(dfProts$FC, base = 2)
dfProts$ProbalogFC <- ifelse(is.na(dfProts$FC) & is.na(dfProts$logFC), NA,
                            ifelse(!is.na(dfProts$FC) & is.na(dfProts$logFC), dfProts$logfcNew,
                                   ifelse(is.na(dfProts$FC) & !is.na(dfProts$logFC), dfProts$logFC,
                                          ifelse(!is.na(dfProts$FC) & !is.na(dfProts$logFC), dfProts$logFC, "Error"
                                                 )
                                          )
                                   )
                            )
# Checking differences between new and original logFC, as well as log10FC just 
# in case it was provided like that
dfProts$log10 <- log(dfProts$FC, base = 10)
dfProts %>% dplyr::select(FC, logFC, logfcNew, ProbalogFC, log10) %>% View
dfError <- dfProts %>% filter(logFC != logfcNew) 
table(dfError$ID) # Not real differences, only a matter of decimals
# 499_339 with extremly values of logFC. This is the only study in which new logFC
# will replace original logFC. For the rest, if logFC was provided in the original
# study, this value will be used. 
dfProts$log10 <- NULL
dfProts$ProbalogFC <- NULL
dfProts$logFC <- ifelse(is.na(dfProts$FC) & is.na(dfProts$logFC), NA,
                             ifelse(!is.na(dfProts$FC) & is.na(dfProts$logFC), dfProts$logfcNew,
                                    ifelse(is.na(dfProts$FC) & !is.na(dfProts$logFC), dfProts$logFC,
                                           ifelse(!is.na(dfProts$FC) & !is.na(dfProts$logFC), dfProts$logfcNew, "Error"
                                           )
                                    )
                             )
)
dfProts$logfcNew <- NULL

#.------- SAVING ------.####
xlsx::write.xlsx(x = dfProts,
                 file = "allProcessedProteinData.xlsx",
                 col.names = T, sheetName = "Qualitative (27 studies)",
                 showNA = F, row.names = F, append = F
                 )




## Final dataframe for qualitative analysis ####
df <- dfProts

### Duplicated protein IDs across datasets ####
length(unique(df$ID)) # 28 resultados de 27 estudios
allFrec <- table(df$ProteinID)
table(allFrec) # 13936 unique proteins; 1481 - 2 studies; 573 - 3 studies; 301 - 4 studies...
interestingFrec <- allFrec[which(allFrec > 1)] # Proteins that appear in, at least, 2 studies
interestingFrec[interestingFrec > 8] # Most frequent proteins
# Keeping proteins that appear in, at least, 2 studies
dfFinal <- df %>% filter(ProteinID %in% names(interestingFrec)) # Filter database (1)
length(unique(dfFinal$ID))

#.------- SAVING ------.####
xlsx::write.xlsx(x = dfFinal,
                 file = "allProcessedProteinData.xlsx",
                 col.names = T, sheetName = "Qualitative (proteins> 1 study)",
                 showNA = F, row.names = F, append = T
                 )


##  Final dataframe for quantitative information ####
# No overlapping proteins between plasma and saliva datasets (both from 1027 reference)
# No aggregation/exclusion is needed, independent results due to lack of overlapping
df1 <- df %>% filter(ID == "1027_15_Sal")
df2 <- df %>% filter(ID == "1027_15_Pla")
df1$ProteinID %in% df2$ProteinID
df2$ProteinID %in% df1$ProteinID

# Removing datasets without p-values or logFC
dfAmanida <- df %>% filter(!(ID %in% c("462_346", "1233_675", "57_73", "1098_977")))
length(unique(dfAmanida$ID)) # 24 comparisons from 23 studies (1097 with saliva and plasma)
# Checking logFC data
sum(is.na(dfAmanida$logFC))
# Cheking duplicated proteinID within each study
for (i in unique(dfAmanida$ID)){
  cat(i, ": \t")
  cat(any(duplicated(dfAmanida[which(dfAmanida$ID == i), "ProteinID"])))
  cat("\n")
}

### Duplicated protein IDs across datasets ##
allFrec <- table(dfAmanida$ProteinID)
table(allFrec) # 13923 unique proteins; 1477 - 2 studies; 574 - 3 studies; 296 - 4 studies...
interestingFrec <- allFrec[which(allFrec > 1)] # Proteins that appear in, at least, 2 studiess
interestingFrec[interestingFrec > 6]  # Most frequent proteins

### Filtering & output ####
dfAmanidaFiltered <- dfAmanida %>% filter(ProteinID %in% names(interestingFrec))


#.------- SAVING ------.####
xlsx::write.xlsx(x = dfAmanida,
                 file = "allProcessedProteinData.xlsx",
                 col.names = T, sheetName = "Quantitative (23 studies)",
                 showNA = F, row.names = F, append = T
)
xlsx::write.xlsx(x = dfAmanidaFiltered,
                 file = "allProcessedProteinData.xlsx",
                 col.names = T, sheetName = "Quantitative (protein>1 study)",
                 showNA = F, row.names = F, append = T
)




