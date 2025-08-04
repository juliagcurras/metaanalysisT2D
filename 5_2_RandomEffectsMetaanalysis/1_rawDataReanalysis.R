###############################################################################-
###################                                       #####################-
###################           RAW DATA ANALYSIS           #####################-
###################                                       #####################-
###############################################################################-


# Julia G Currás  - 2025/05/13

# Setup ####
library(dplyr)
library(ggplot2)
set.seed(9396)

doTestT <- function(df, dfGrupos, g1, g2) {
  
  # g1: control
  df <- as.data.frame(df)
  dfGrupos <- as.data.frame(dfGrupos)
  
  # Selecting columns from dataframe of each group
  samplesG1 <- as.character(dfGrupos[dfGrupos$Groups == g1, "Samples"])
  samplesG2 <- as.character(dfGrupos[dfGrupos$Groups == g2, "Samples"])
  df <- df[, c(samplesG1, samplesG2)]
  
  # log2FC: group2 (case) - group1 (control)
  df$FC  <- apply(df, 1, function(x) {
    mean(x[samplesG2], na.rm =T)/mean(x[samplesG1], na.rm = T)
  }
  )
  df$logFC <- log(df$FC, base = 2)
  
  # Fifference of means estimation
  df$DiffExpr  <- apply(df, 1, function(x) {
    (mean(x[samplesG2], na.rm =T) - mean(x[samplesG1], na.rm = T))
  }
  )
  
  # T-test with equal variance
  df[, c( "t", "P.Val")] <- t(apply(df, 1, function(x) {
    # When there isn't enough observations, returning NA
    evalG1NA <- sum(is.na(x[samplesG1])) #/length(x[samplesG1]) 
    evalG2NA <- sum(is.na(x[samplesG2]))
    if (((length(x[samplesG1])-evalG1NA) < 3) | 
        ((length(x[samplesG2])-evalG2NA) < 3 ) #|
    ) {
      return(c(NA, NA))
    }
    # Otherwise, test
    testDone <- try(stats::t.test(x[samplesG2], x[samplesG1],
                                  alternative = "two.sided",
                                  var.equal = TRUE), T)
    if (!(class(testDone) == "try-error")){ # if some error happend, then NA
      res <- testDone
      return(c(res$statistic, res$p.value))
    } else {
      return(c(NA, NA)) 
    }
  }
  ))
  
  #Benjamini-Hochberg correction for multiple testing
  df$adj.P.Val <- stats::p.adjust(df$P.Val, method = "BH")
  df <- df[, c("logFC", "FC", "P.Val", "adj.P.Val", "DiffExpr")]
  
  # Returning final vector
  return(df)
}


removeDiscordancesReplicates <- function(dfQ, nReps = 2){
  posFinal <- nReps
  dfResult <- NULL
  for (i in seq(1, ncol(dfQ)-1, nReps)){ # for each sample
    if (i != 1){
      posFinal <- posFinal + nReps
    } 
    dfAux <- dfQ[, i:posFinal] # its respective columns are selected
    
    resultProt <- apply(dfAux, 1, function(fila){ # and then only one value is established
      totalNA <- sum(is.na(fila)) 
      if (totalNA >=3){ # NA if there are 3 or 4 missing values
        return(NA)
      } else {
        return(mean(fila, na.rm = T)) # mean(replicates) if there are 0, 1 or 2 missing values
      }
    })
    
    # Finally, info from different samples is bound
    if (is.null(dfResult)){
      dfResult <- data.frame(resultProt)
    } else {
      dfResult <- cbind(dfResult, resultProt)
    }
  }
  colnames(dfResult) <- colnames(dfQ)[(seq(nReps, ncol(dfQ), nReps))]
  
  # Final object!
  return(dfResult)
}




#.------------------------------------------------------------------------ ####
# MAIN ANALYSIS OF RAW DATA #########
#.------------------------------------------------------------------------ ####


# 284_268 ####
## logFC estimations ####
# Group 2 = T2DM and Group 1 = control - inferred from Figure 5c
data <- xlsx::read.xlsx(file = "rawData/284_268_PubMed.xlsx", sheetIndex = 2)
colnames(data)[1:50] <- paste0(data[2, 1:50], "_S", 1:50)
colnames(data)[1] <- "1_S1"
data <- data[-c(1,2),]

# Selecting groups and good proteins
df <- data %>% 
  select(starts_with("1_"), starts_with("2_"), Potential.contaminant, 
         Student.s.T.test.Significant.2_1, 
         Student.s.T.test.Significant.2_1_, 
         Student.s.T.test.p.value.2_1:Student.s.T.test.Test.statistic.2_1, 
         Protein.IDs:Gene.names)
table(df$Potential.contaminant)
df <- df %>% filter(is.na(Potential.contaminant))

# P-values from two test-t, the first one was selected
table(df$Student.s.T.test.Significant.2_1)
table(df$Student.s.T.test.Significant.2_1_)

# Log2FC estimation #
df$FC <- apply(df[, c(1:16)], 1, function(i){
  g2 <- as.numeric(i[8:16]) # T2DM 
  g1 <- as.numeric(i[1:7]) # Control
  return(mean(g2, na.rm = T)/mean(g1, na.rm = T))
})
df$logFC <- log(df$FC, 2)

# P-value adjustment #
df$pvalAdj <- p.adjust(p = df$Student.s.T.test.p.value.2_1, 
                       method = "BH")

# Selecting outputs & preparing final database to export #
df %>% filter(pvalAdj < 0.05) %>% nrow
dfFinal <- df %>% 
  filter(Student.s.T.test.p.value.2_1 != "NaN") %>%
  select(Majority.protein.IDs, Gene.names, Protein.names,
         logFC, FC, Student.s.T.test.p.value.2_1, pvalAdj, 
         `1_S1`:`2_S30`)
dfFinal[, 5:ncol(dfFinal)] <- apply(dfFinal[, 5:ncol(dfFinal)], 2, as.numeric)
colnames(dfFinal) <- gsub(pattern = "1_", replacement = "Control_", x = colnames(dfFinal))
colnames(dfFinal) <- gsub(pattern = "2_", replacement = "T2DM_", x = colnames(dfFinal))
# Save
xlsx::write.xlsx(dfFinal, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = F, 
                 sheetName = "284_268")




# 377_100 ####
## Reanalysis ####
# Loading data from TS5 #
tab5 <- xlsx::read.xlsx(file = "rawData/377_100_PubMed_table5.XLSX", sheetIndex = 1)
colnames(tab5) <- tab5[1,] 
tab5 <- tab5[-1,]

# Selecting study groups #
data <- tab5 %>% 
  dplyr::select(Accession_id, Description, 
                starts_with("Abundance_DM"), 
                starts_with("Abundance_NC"))

# Keeping only abundance information #
dfQuant <- data[, -c(1,2)]
rownames(dfQuant) <- data$Accession_id
colnames(dfQuant) <- gsub(colnames(dfQuant), pattern = "Abundance_", replacement = "")
dfQuant[dfQuant == 0] <- NA

# Design matrix #
dfGrupos <- data.frame(
  Samples = colnames(dfQuant),
  Groups = gsub(colnames(dfQuant), pattern = "[0-9]", replacement = "")
)

# Removing proteins with more than 50% missing data in, at least, one group #
g1 <- "NC" # Control
g2 <- "DM" # T2DM
naLim <- 0.5
samplesG1 <- as.character(dfGrupos[dfGrupos$Groups == g1, "Samples"])
samplesG2 <- as.character(dfGrupos[dfGrupos$Groups == g2, "Samples"])
# dfQuant <- dfQuant[, c(as.character(samplesG1), as.character(samplesG2))]
indOkProts <- apply(dfQuant, 1, function(x) {
  evalG1NA <- sum(is.na(x[samplesG1]))/length(x[samplesG1]) 
  evalG2NA <- sum(is.na(x[samplesG2]))/length(x[samplesG2])
  return(!any(evalG1NA >= naLim, evalG2NA >= naLim)) # returning only the good ones
})
dfQuant <- as.data.frame(dfQuant)
df <- dfQuant[which(indOkProts),]
proteinas <- rownames(df) 

# performing comparisons #
df <- apply(df, 2, as.numeric)
rownames(df) <- proteinas
df <- as.data.frame(df)
tabRes <- doTestT(df = df, dfGrupos = dfGrupos, g1 = g1, g2 = g2)
tabRes$ID <- rownames(tabRes) 
tabRes$Description <- data[which(data$Accession_id %in% tabRes$ID), "Description"]

# Adding raw data info to test results #
df$ID <- rownames(df)
dfAll <- merge(tabRes, df, by = "ID", all=T)
dfAll <- dfAll %>% 
  arrange(P.Val) %>%
  select(ID, Description, everything())

# Saving 
xlsx::write.xlsx(dfAll, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = T, 
                 sheetName = "377_100")


# Adding significant proteins based on missing data: only for qualitative analysis #
# Proteins with more than 50% missing data in only one group are considered 
# as significant:
  # Upregulated if missing data > 50% in control group
  # Downregulared if missing data > 50% in T2DM group

# Retrieving proteins with more than 50% missing values in at least one group
df_DEPs <- dfQuant[!(unname(indOkProts)),]

# Estimation of % missign data by group
missingByGroup <- apply(df_DEPs, 1, function(x) {
  evalG1NA <- sum(is.na(x[samplesG1]))/length(x[samplesG1]) 
  evalG2NA <- sum(is.na(x[samplesG2]))/length(x[samplesG2])
  return(c(G1 = evalG1NA, G2 = evalG2NA))
}, simplify = T)
missingByGroup <- as.data.frame(t(missingByGroup))
missingByGroup$ID <- rownames(missingByGroup)
df_DEPs$ID <- rownames(df_DEPs)
dfDEPs <- merge(df_DEPs, missingByGroup, by = "ID")

# Identification of significant proteins and classification based on the direction
# of change using the % missing data
dfDEPs$FCDir <- ifelse(dfDEPs$G1 < 0.5 & dfDEPs$G2 > 0.5, "Down", # más valores perdidos en T2DM que en control
                       ifelse(dfDEPs$G1 > 0.5 & dfDEPs$G2 < 0.5, "Up", "Out")) # más valores perdidos en control que en T2DM

# Out: proteins with more than 50% missing values in both groups
table(dfDEPs$FCDir)
dfDEPs <- dfDEPs %>% filter(FCDir != "Out") %>% dplyr::select(ID, FCDir)

# Adding protein information
dfBasal <- data %>% 
  dplyr::select(Accession_id, Description) %>% 
  rename(ID = Accession_id)
dfDEPs <- merge(dfDEPs, dfBasal, by = "ID", all.x = T)

# Merging cualitative and cuantitative results
tabRes <- tabRes %>% filter(P.Val < 0.05) 
tabRes$FCDir <- ifelse(tabRes$logFC < 0, "Down", 
                       ifelse(tabRes$logFC > 0, "Up", "Ups"))
table(tabRes$FCDir)
tabRes <- tabRes %>% dplyr::select(ID, Description, FCDir)
dfAll <- merge(dfDEPs, tabRes, by = c("ID", "FCDir", "Description"), all = T)
table(duplicated(dfAll$ID)) # checking duplicated proteins





# 432_287 ####
## logFC estimation & p-value adjustment ####
# Loading data #
datos <- xlsx::read.xlsx(file = "rawData/432_287_PubMed.xlsx", 
                         sheetIndex = 1)

# Selecting samples of interest #
tab <- datos %>% select(Protein.accession, Gene.name, Protein.description, 
                        C1:B10, C.B.Ratio, C.B.P.value) # C control / B T2DM

# Estimation of C/B instead of B/C (T2DM/Control)
tab$BCratio <- 1/tab$C.B.Ratio # B/C => 1/(C/B)

# log FC estimation #
tab$logFC <- log(tab$BCratio, base = 2)

# Adjusting p-valkues using the BH correction
tab$pAdjusted <- p.adjust(tab$C.B.P.value, method = "BH")

# Final data
tab <- tab %>% 
  dplyr::select(Protein.accession, Gene.name, Protein.description,
                logFC, BCratio, C.B.P.value, pAdjusted, C1:B10) %>%
  arrange(C.B.P.value)

# Saving #
writexl::write_xlsx(tab, "temp.xlsx")
xlsx::write.xlsx(tab, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = T, 
                 sheetName = "432_287")




# 440_211 ####

# MAIN PROBLEM: IDENTIFYING THE GROUPS IN THE RAW DATA (PROTEINGROUP.TXT)

## First approximation ####
# Initial management of dataframe #
## Loading proteins used in the quantitative analysis (682)
proteinas <- xlsx::read.xlsx(file = "rawData/440_211_Table 2.XLSX", 
                             sheetIndex = 1)
colnames(proteinas) <- proteinas[1,]
proteinas <- proteinas[-1,]
proteinas <- proteinas$`Protein IDs`
## Loading data from proteinGroups.txt
data <- read.table(file = "rawData/440_211_PubMed.txt", 
                   header = T, sep = "\t")
df <- data %>% dplyr::select(Protein.IDs, Protein.names, Gene.names, 
                             starts_with("LFQ.intensity."), Reverse)
# Set 0 to NA 
df[df == 0] <- NA

# Filters #
## Keep only interesting proteins
df <- df %>%
  dplyr::filter(Reverse != "+") %>%
  dplyr::select(-Reverse) 
## Selecting proteins previously reported by the authors in TableS2 #
df <- df %>% filter(Protein.IDs %in% proteinas) # 682 with complete information in 3 samples at least in one group
colnames(df) <- gsub(x = colnames(df), pattern = "LFQ.intensity.", replacement = "")

# New dataframe without protein information columns #
dfQuant <- df %>% dplyr::select(Group1_rep1:Group4_rep8) # 2vs4
rownames(dfQuant) <- df$Protein.IDs

# Design matrix #
dfGrupos <- data.frame(
  Samples = colnames(dfQuant),
  Groups = substr(colnames(dfQuant), start = 1, stop = 6)
)
table(dfGrupos$Groups)

# Log transformation #
dfLog <- as.data.frame(log(dfQuant, base = 2))
proteinas <- rownames(dfLog)

# Imputation from normal distribution #
dfLog  <- as.matrix(dfLog)
ancho <- 0.3
downShift <- 1.8
s = sd(dfLog , na.rm = T)*ancho
m = mean(dfLog , na.rm = T) - downShift*sd(dfLog , na.rm = T)
dfFinal <- as.data.frame(dfLog)
for (i in 1:nrow(dfLog)){
  for (j in 1:ncol(dfLog )){
    dfFinal[i, j] <- ifelse(is.na(dfLog[i,j]), rnorm(n = 1, mean =m, s = s), dfLog[i,j])
  }
}
rownames(dfFinal) <- proteinas

# Test-t with BH correction between all pairs #
tabRes12 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group2", g2 = "Group1")
tabRes23 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group2", g2 = "Group3")
tabRes13 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group1", g2 = "Group3")
tabRes34 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group4", g2 = "Group3")
tabRes14 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group1", g2 = "Group4")
tabRes24 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group2", g2 = "Group4")
## Number of DEPs? 
tabRes12 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
tabRes23 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
tabRes13 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
tabRes34 %>% filter(adj.P.Val<0.05) %>% nrow() # 39
tabRes14 %>% filter(adj.P.Val<0.05) %>% nrow() # 2
tabRes24 %>% filter(adj.P.Val<0.05) %>% nrow() # 1
# DIFFERENT RESULTS FROM THOSE REPORTED IN THE MAIN PAPER

# HBB: only DEP between T2DM vs Control. Exploration:
df %>% filter(Gene.names == "HBB") %>% View
tabRes2[which(rownames(tabRes12) == "P68871"), ] #12
tabRes4[which(rownames(tabRes14) == "P68871"), ] #13
tabRes3[which(rownames(tabRes13) == "P68871"), ] #14
tabRes23[which(rownames(tabRes23) == "P68871"), ] #23
tabRes34[which(rownames(tabRes34) == "P68871"), ] #43
tabRes24[which(rownames(tabRes24) == "P68871"), ] #24



## Second approximation ####
# Keeping only complete data (not affected by imputation) #
dfQuant_Complete <- na.omit(dfQuant)

# log transformation #
dfLog <- as.data.frame(log(dfQuant_Complete, base = 2))
dfFinal <- dfLog

# Test-t with BH correction between all pairs #
tabRes12 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group2", g2 = "Group1")
tabRes23 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group2", g2 = "Group3")
tabRes13 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group1", g2 = "Group3")
tabRes34 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group4", g2 = "Group3")
tabRes14 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group1", g2 = "Group4")
tabRes24 <- doTestT(df = dfFinal, dfGrupos = dfGrupos, g1 = "Group2", g2 = "Group4")
# These results correspond with the ones obtained by Perseus

## Number of DEPs? 
tabRes12 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
tabRes23 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
tabRes13 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
tabRes34 %>% filter(adj.P.Val<0.05) %>% nrow() # 58
tabRes14 %>% filter(adj.P.Val<0.05) %>% nrow() # 1
tabRes24 %>% filter(adj.P.Val<0.05) %>% nrow() # 0
# Any clue to find the identifity of groups.

# HBB: only DEP between T2DM vs Control. What happen with this protein?
df["P68871", ] # Protein without complete data

# Boxplot comparison between Figure S4 and the ones made with the raw data from
# proteinGroup.txt
## Finding the proteins from S4 with complete data
dfFinal$Protein <- rownames(dfFinal)
dfFinal$Genes <- df[df$Protein.IDs %in% rownames(dfFinal), "Gene.names"]
genesS4 <- c("HBB2", "S100P", "S100A11", "S100A13", "CSTB", "NAMPT", "CALML3", 
             "CALM3", "CALML5", "APOBEC3A", "NCCRP1", "CLRX", "CTSL", "GLUL",
             "SERPINF2", "PDIA3", "TXNRD1", "ABHD14B", "MTPN", "PDIA3", 
             "EEF1B2", "PREP", "SH3BGRL", "GSN", "PGD", "NQO2", "CFL1", 
             "GALM", "IMPA1", "EIF5A/EIF5AL1", "PLIN3", "TXNDC17")
proteinas <- dfFinal %>% filter(Genes %in% genesS4) %>% pull(Protein)
geneName <- dfFinal %>% filter(Genes %in% genesS4) %>% pull(Genes)
dfPlot <- as.data.frame(t(dfFinal[proteinas, -c(33:34)]))
dfPlot$Groups <- dfGrupos$Groups
dfPlot$Groups <- factor(dfPlot$Groups, 
                        levels = c("Group4", "Group2", "Group3","Group1"))
dfPlot$Colores <- rep(c("darkgreen", "lightyellow", "orange", "darkred"), each = 8)

# Displaying figures
# S100P - P25815
ggplot(dfPlot, aes(x = Groups, y = P25815, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) +  ylim(c(26, 31)) +
  theme_minimal() + ylab ("S100P - P25815") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

# P04080 - CSTB
ggplot(dfPlot, aes(x = Groups, y = P04080, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() + ylim(c(28, 34)) + ylab ("CSTB - P04080") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

# P06396;CON__Q3SX14;REV__Q6TDU7 - GSN
ggplot(dfPlot, aes(x = Groups, y = `P06396;CON__Q3SX14;REV__Q6TDU7`, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + 
  ylim(c(28, 34)) + ylab ("GSN - P23528;Q9Y281") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

# CFL1 - P23528;Q9Y281
ggplot(dfPlot, aes(x = Groups, y = `P23528;Q9Y281`, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + 
  ylim(c(28, 32)) + ylab ("CFL1 - P23528;Q9Y281") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

# S100A11 - P31949
ggplot(dfPlot, aes(x = Groups, y = P31949, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + 
  ylim(c(28, 34)) + ylab ("S100A11 - P31949") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

# PGD - P52209
ggplot(dfPlot, aes(x = Groups, y = P52209, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + 
  ylim(c(28, 32)) + ylab ("PGD - P52209") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

# MTPN - P58546
ggplot(dfPlot, aes(x = Groups, y = P58546, fill = Groups)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + 
  ylim(c(26, 30)) + ylab ("MTPN - P58546") +
  geom_point(size = 3) + xlab("") +
  scale_fill_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                               "Group3" = "orange", "Group1" = "darkred")) +
  scale_color_manual(values = c("Group4" = "darkgreen", "Group2" = "lightyellow", 
                                "Group3" = "orange", "Group1" = "darkred")) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black", size = 25),
        axis.title = element_text(size = 25), 
        legend.position = "none", 
        text = element_text(family = "Calibri"), 
        axis.line = ggplot2::element_line(linewidth = 0.5, colour = "black"), 
        axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"))

################ Final correspondence of groups:
# cONTROL = GRUPO 4 / T2DM = GRUPO 2 / NPDR = GROUP 3 / PDR = GROUP 1
# Assigment:
dfPlot$Groups <- factor(dfPlot$Groups, 
                        levels = c("Group4", "Group2", "Group3","Group1"),
                        labels = c( "Control", "T2DM", "NPDR","PDR"))

## Data analysis ####
# Loading raw data
data <- read.table(file = "rawData/440_211_PubMed.txt", 
                   header = T, sep = "\t")
df <- data %>% dplyr::select(Protein.IDs, Protein.names, Gene.names, 
                             starts_with("LFQ.intensity."), Reverse)
# Set 0 to NA 
df[df == 0] <- NA

# Filters 
# Keep´only interesting proteins
df <- df %>%
  dplyr::filter(Reverse != "+") %>%
  dplyr::select(-Reverse) 
colnames(df) <- gsub(x = colnames(df), pattern = "LFQ.intensity.", replacement = "")

# selecting samples from group 2 (T2DM) and 4 (control) #
dfQuant <- df %>% dplyr::select(Group4_rep1:Group4_rep8, Group2_rep1:Group2_rep8) # 2vs4
colnames(dfQuant) <- gsub(pattern = "Group2", replacement = "T2DM", x = colnames(dfQuant))
colnames(dfQuant) <- gsub(pattern = "Group4", replacement = "Control", x = colnames(dfQuant))

# Building design matrix #
rownames(dfQuant) <- df$Protein.IDs
dfGrupos <- data.frame(
  Samples = colnames(dfQuant),
  Groups = substr(colnames(dfQuant), start = 1, stop = nchar(colnames(dfQuant)) - 5)
)
table(dfGrupos$Groups)

# Excluding proteins with more than 2 missing values in any group #
g1 <- "Control"
g2 <- "T2DM"
naLim <- 4/8
samplesG1 <- as.character(dfGrupos[dfGrupos$Groups == g1, "Samples"])
samplesG2 <- as.character(dfGrupos[dfGrupos$Groups == g2, "Samples"])
indOkProts <- apply(dfQuant, 1, function(x) {
  evalG1NA <- sum(is.na(x[samplesG1]))/length(x[samplesG1]) 
  evalG2NA <- sum(is.na(x[samplesG2]))/length(x[samplesG2])
  return(!any(evalG1NA >= naLim, evalG2NA >= naLim)) # returning only the good ones
})
dfQuant <- as.data.frame(dfQuant)
dfQuant <- dfQuant[which(indOkProts),]

# Log transformation #
df <- as.data.frame(log(dfQuant, base = 2))

# Test t with BH correction #
tabRes <- doTestT(df = df, dfGrupos = dfGrupos, g1 = g1, g2 = g2)
tabRes$Protein <- rownames(tabRes)
tabRes$Genes <- data[which(data$Protein.IDs %in% tabRes$Protein), "Gene.names"]
tabRes$Description <- data[which(data$Protein.IDs %in% tabRes$Protein), "Protein.names"]

# Adding raw data for metanalisis #
# Add raw data info
df$Protein <- rownames(df)
dfAll <- merge(tabRes, df, by = "Protein", all=T)
dfAll <- dfAll %>% 
  arrange(P.Val) %>%
  select(Protein, Genes, Description, everything())

# Saving #
xlsx::write.xlsx(dfAll, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = T, 
                 sheetName = "440_211")





# 944_401 ####
## Reanalysis ####
# Loading data from proteinGroups.txt #
data0 <- read.table(file = "rawData/944_401_Scopus.txt", 
                    header = T, sep = "\t")
data <- data0 %>% dplyr::select(Protein.IDs, Protein.names, Gene.names, 
                                starts_with("LFQ.intensity."), Reverse, 
                                Potential.contaminant, Only.identified.by.site)
# Set 0 to NA #
data[data == 0] <- NA

# Keep only interesting proteins #
data <- data %>%
  dplyr::filter(Reverse != "+") %>%
  dplyr::filter(Only.identified.by.site!= "+") %>%
  dplyr::filter(Potential.contaminant!= "+") %>%
  dplyr::select(-Reverse, -Potential.contaminant, -Only.identified.by.site) 
colnames(data) <- gsub(x = colnames(data), pattern = "LFQ.intensity.", replacement = "")

# Select only sample information #
dfQuant <- data %>% dplyr::select(`023_01_C6_1_3090`:`197_04_F5_1_3081`) 
table(duplicated(data$Protein.IDs))
rownames(dfQuant) <- data$Protein.IDs

# Desing matrix: replicates and groups #
dfReplicates <- data.frame(
  Samples = colnames(dfQuant), 
  Replicates = sapply(strsplit(colnames(dfQuant), "_"), "[[", 1),
  Groups = rep(c("C", "C", "DM", "DM", "C", "DM", "DM", "C", "C", "DM"), each=4)
)
dfGroups <- data.frame(
  Samples = unique(dfReplicates$Replicates),
  Groups = c("C", "C", "DM", "DM", "C", "DM", "DM", "C", "C", "DM")
)
table(dfGroups$Groups)
table(dfGroups$Groups, dfGroups$Samples)

# Dealing with technical replicates for the sample sample #
# If nº of missing data >=3 => NA
# If nº of missing data <3 => mean
dfRep <- removeDiscordancesReplicates(dfQuant, nReps = 4)
# dfRep <- dfQuant
colnames(dfRep) <- sapply(strsplit(colnames(dfRep), "_"), "[[", 1)

# Log transformation #
dfLog <- as.matrix(log(dfRep, base = 2))

# Dealing with missing data #
# Option 1: keeping only complete data #
dfFinal <- na.omit(as.data.frame(dfLog)) # 995 proteins

# Test t with BH correction #
tabRes <- doTestT(df = dfFinal, dfGrupos = dfGroups, g1 = g1, g2 = g2)
tabRes$Protein <- rownames(tabRes)
tabRes$Genes <- data[which(data$Protein.IDs %in% tabRes$Protein), "Gene.names"]
tabRes$Description <- data[which(data$Protein.IDs %in% tabRes$Protein), "Protein.names"]
# Volcano plot to compare qith Figure 1D
Biostatech::plotVolcano(data = tabRes, xVar = "DiffExpr", 
                        yVar = "P.Val", minFC = 0.75, sizeCirculos = 1.2,  
                        maxAdjP = 0.05, proteins = "Genes")$grafico

# Adding raw data for metanalisis #
colnames(dfFinal) <- paste0(colnames(dfFinal), "_", dfGroups[colnames(dfFinal), "Groups"])
dfFinal$Protein <- rownames(dfFinal)
dfAll <- merge(tabRes, dfFinal, by = "Protein", all=T)
dfAll <- dfAll %>% 
  arrange(P.Val) %>%
  select(Protein, Genes, Description, everything())

# Saving #
xlsx::write.xlsx(dfAll, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = T, 
                 sheetName = "944_401")




# 499_339 ####
## Reanalysis ####
# Loading data from Table S1 #
data0 <- xlsx::read.xlsx(file = "rawData/499_339_PubMed.xls",
                         sheetIndex = 1)
# Selecting interesting information #
data <- data0 %>%
  dplyr::select(Protein, Gene.name, Description, Sample.1:Sample.19)
colnames(data) <- c(colnames(data)[1:3], paste0("ND", 1:4), paste0("DM", 1:15))
data <- data[-1,]
table(duplicated(data$Protein)) # Checking the presence of duplicated proteins #

# Keeping only samples and abundance values #
dfQuant <- data %>% select(ND1:DM15)
rownames(dfQuant) <- data$Protein
dfGrupos <- data.frame(
  Samples = colnames(dfQuant),
  Groups = gsub(colnames(dfQuant), pattern = "[0-9]", replacement = "")
)
g1 <- "ND"
g2 <- "DM"

# Onyl complete data #
df <- na.omit(dfQuant)
proteinas <- rownames(df)
df <- as.data.frame(apply(df, 2, as.numeric))
rownames(df) <- proteinas

# Performing comparisons #
tabRes <- doTestT(df = df, dfGrupos = dfGrupos, g1 = g1, g2 = g2)
tabRes$ID <- rownames(tabRes) 
tabRes$Description <- data[which(data$Protein %in% tabRes$ID), "Description"]
tabRes$Gene.name <- data[which(data$Protein %in% tabRes$ID), "Gene.name"]

# Adding raw data information to results #
df$ID <- rownames(df)
dfAll <- merge(tabRes, df, by = "ID", all=T)
dfAll <- dfAll %>% 
  arrange(P.Val) %>%
  select(ID, Gene.name, Description, everything())

# Saving #
xlsx::write.xlsx(dfAll, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = T, 
                 sheetName = "499_339")






# 1046_653 ####
## Reanalysis ####
# Loading TS1 #
data <- xlsx::read.xlsx(
  file = "rawData/1046_653_Scopus.xlsx", 
  sheetIndex = 1)

# Selecting interesting columns #
df <- data %>% 
  dplyr::select(Fasta.headers, Protein.name, Gene.name, 
                O.P.1_fasting:O.P.9_fasting,
                D.P.1_fasting:D.P.10_fasting)
df$Proteins <- sapply(strsplit(df$Fasta.headers, split = "|", fixed = T), "[[", 2)
table(duplicated(df$Proteins)) # Keeping all proteins
df <- df %>% dplyr::select(Proteins, everything(), -Fasta.headers)

# Design matrix #
dfGrupos <- 
  data.frame(
    Samples = colnames(df)[4:ncol(df)],
    Groups = rep(c("NO", "DO"), 9:10))
g1 <- "NO" # Control
g2 <- "DO" # Diabetes

# Removing duplicated Gene.names
dfQuant <- df[, -c(1:3)]
rownames(dfQuant) <- df$Proteins

# Performing comparisons #
tabRes <- doTestT(df = dfQuant, dfGrupos = dfGrupos, g1 = g1, g2 = g2)
tabRes$Proteins <- rownames(tabRes)
# tabRes$Genes <- df[which(df$Proteins %in% tabRes$Proteins), "Gene.name"]
# tabRes$Description <- df[which(df$Proteins %in% tabRes$Proteins), "Protein.name"]

# Merge raw data with results #
dfAll <- merge(tabRes, df, by = "Proteins", all = T)
dfAll <- dfAll %>% 
  dplyr::select(Proteins, Gene.name, Protein.name, everything())

# Saving #
xlsx::write.xlsx(dfAll, file = "1_reanalysis_RawData.xlsx", 
                 row.names = F, showNA = F, col.names = T, append = T, 
                 sheetName = "1046_653")




#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX IMPORTANT xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
# 1_reanalysis_RawData.xlsx was manually modified to add raw data from studies 
# that did not need a reanalysis (16_128, 231_253, 425_58, 446_99, 511_252) and
# to unify the format across all excel sheets (number of initial columns and names)
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####






























