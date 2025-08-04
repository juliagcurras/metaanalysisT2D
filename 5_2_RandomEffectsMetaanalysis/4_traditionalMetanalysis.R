##############################################################################-

#################       METANALYSIS - MAIN STRUCTURE       ###################-

##############################################################################-


# Julia G Currás - 2025/05/28
rm(list=ls())
graphics.off()
library(dplyr)
library(Biostatech)
library(meta)

# SMD and SE estimation
getSMD <- function(df, proteina){
  # Get data
  df <- df[proteina,]
  dm <- df %>% dplyr::select(starts_with("T2DM_")) %>% as.vector %>% unlist
  nDM <- length(dm)
  control <- df %>% dplyr::select(starts_with("Control_")) %>% as.vector %>% unlist
  nControl <- length(control)
  
  # get SMD and SE pooled
  result <- esc::esc_mean_sd(grp1m = mean(dm), grp1sd = sd(dm), grp1n = nDM,
                             grp2m = mean(control), grp2sd = sd(control), 
                             grp2n = nControl)
  return(list(SMD = result$es, SE = result$se))
}


# Loading curated, raw data ####
lista <- readRDS(file = "2_allRawDataProcessed.rds")

# Loading interesting proteins ####
hojas <- readxl::excel_sheets("3_proteinsForMetanalisis.xlsx")
allDatasets <- list()
for (i in hojas){
  allDatasets[[i]] <- readxl::read_xlsx(path = "3_proteinsForMetanalisis.xlsx", 
                  sheet = i, col_names = T)
  
}

# Loading extra data ####
basal <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", 
                                         sheet = "1. Basal data", col_names = T)
sampleInfo <- basal %>% 
  dplyr::select(ID, Reference, `Type of sample (grouped)`, `Data acquisition mode`) %>%
  rename(Sample = `Type of sample (grouped)`)
sampleInfo$Sample <- as.character(sampleInfo$Sample)
sampleInfo$Sample <- ifelse(sampleInfo$Sample == "Pancreas biopsy (pancreatic island cells)", 
                                    "Pancreas", 
                                    ifelse(sampleInfo$Sample == "Extracellular vesicles (blood)", "Blood (EV)", 
                                    ifelse(sampleInfo$Sample == "Serum", "Blood (Serum)", 
                                    ifelse(sampleInfo$Sample == "Plasma", "Blood (Plasma)", 
                                           sampleInfo$Sample))))


# Meta-analysis ####
allMeta <- list()

for(j in hojas){
  # j <- hojas[2]
  data <- as.data.frame(allDatasets[[j]])
  proteinas <- unique(data$ProteinID)
  nombre <- gsub(j, pattern = "[^A-Za-z]", replacement = "")
  # proteina <- proteinas[1]
  for (proteina in proteinas){
    # Select info
    ids <- data[which(data$ProteinID == proteina), "ID"]
    selectedList <- lista[ids]
    
    # Get SDM and SE
    result <- lapply(selectedList, getSMD, proteina)
    df <- data.frame(
      ID = names(result),
      SMD = sapply(result, "[[", 1),
      SE = sapply(result, "[[", 2),
      row.names = NULL
    )
    
    df <- merge(df, sampleInfo, by = "ID", all.x = T)
    
    # Metanalysis 
    m.gen <- metagen(TE = df$SMD,
                     seTE = df$SE,
                     studlab = df$Reference,
                     data = df,
                     sm = "SMD",
                     common = FALSE,
                     random = TRUE,
                     method.tau = "REML",
                     method.random.ci = "HK",
                     title = proteina)
    
    # Save
    finalName <- paste0(nombre, "_", proteina)
    allMeta[[finalName]] <- m.gen
  }
}


## Saving ####
newAppend <- F
for (protein in names(allMeta)){
  res <- allMeta[[protein]]
  # Extract info studies
  df_estudios <- data.frame(
    Study = res$studlab,
    SMD = res$TE,
    SD_pooled = res$seTE,
    IC95_inf = res$lower,
    IC95_sup = res$upper,
    TypeOfSample = res$data$Sample,
    DataAcquisition = res$data$`Data acquisition mode`
  )
  # Extract aggregated results (SMD)
  df_global <- data.frame(
    Study = "Global (random effects)",
    SMD = res$TE.random,
    SD_pooled = res$seTE.random,
    IC95_inf = res$lower.random,
    IC95_sup = res$upper.random
  )
  # Extract heterogenity measures
  df_hetero <- data.frame(
    Study = "Heterogeneity",
    SMD = NA,
    SD_pooled = NA,
    IC95_inf = NA,
    IC95_sup = NA,
    tau2 = res$tau^2,
    se_tau2 = res$se.tau,
    I2 = res$I2,
    lower_I2 = res$lower.I2,
    upper_I2 = res$upper.I2,
    pval_Q = res$pval.Q
  )
  
  # Bind everything
  df_final <- dplyr::bind_rows(
    df_estudios,
    df_global,
    df_hetero
  )
  
  # Saving
  xlsx::write.xlsx(x = df_final, file = "4_resultMetanalysis.xlsx", 
                   sheetName = protein, col.names = T, row.names = F, 
                   append = newAppend, showNA = F)
  newAppend <- T
}

## Code for display results ####
p1 <- meta::forest(allMeta$All_P62937, 
             sortvar = Sample,
             leftcols = c("studlab", "Sample", "effect", "ci", "w.random"),
             rightcols = FALSE,
             leftlabs = c("Author", "Sample", "SMD", "95%-CI", "Weigth"),
             prediction = TRUE, 
             print.I2 = TRUE,
             print.I2.ci = TRUE,
             print.tau2 = TRUE,
             print.tau2.ci = TRUE)
p1




















