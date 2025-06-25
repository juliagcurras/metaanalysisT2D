##############################################################################-

#################       METANALYSIS - DATASETS       ########################-

##############################################################################-


# Julia G Currás - 2025/05/28
rm(list=ls())
graphics.off()
library(dplyr)

rawData <- c("16_128", "284_268", "377_100", "425_58", "432_287", "446_99", 
             "440_211", "511_252", "944_401", "231_253", "499_339", "1046_653")



# Loading data #### 
directorio <- "G:/Mi unidad/_General/DoctoradoIndustrial2022/Publicaciones/METAANALISIS/Data/3_Extraction/"
archivos <- c("proteinData_1_31.xlsx", "proteinData_29_56.xlsx", "proteinData_57_84.xlsx")

dfProts <- list()
for (j in archivos){
  cat(j, ": ")
  ruta <- paste0(directorio, j)
  hojas <- readxl::excel_sheets(ruta)[-1]
  hojas <- hojas[hojas %in% rawData]
  for (i in hojas){
    cat(i, "| ")
    df <- readxl::read_xlsx(path = ruta, sheet = i, col_names = T)
    df <- as.data.frame(df)
    if (i == "511_252"){
      df <- df[, c(2, 9:ncol(df))]
    } else {
      rownames(df) <- df$ProteinID
      df <- df[, c(9:ncol(df))]
    }
    dfProts[[i]] <- df
  }
  cat("\n")
}
dfProtsOriginal <- dfProts
dfProts <- dfProtsOriginal

# Processing by study ####
## 16_128 ####
data <- dfProts$`16_128`
colnames(data)
colnames(data) <- c(paste0("Control_", colnames(data)[1:3]), 
                    paste0("T2DM_", colnames(data)[4:6]))
colnames(data)
dfProts[["16_128"]] <- data


## 284_268 ####
data <- dfProts$`284_268`
colnames(data)
colnames(data) <- paste0(colnames(data), "__", 1:ncol(data))
colnames(data) <- gsub(colnames(data), pattern = "\\.\\.\\.[0-9]+_", replacement = "")
colnames(data)
data <- dfProts$`284_268`


## 377_100 ####
data <- dfProts$`377_100`
colnames(data)
colnames(data) <- c(paste0("Control_", colnames(data)[22:ncol(data)]), 
                    paste0("T2DM_", colnames(data)[1:21]))
colnames(data)
dfProts[["377_100"]] <- data


## 425_58 ####
data <- dfProts$`425_58`
colnames(data)
colnames(data) <- gsub(colnames(data), replacement = "T2DM", pattern = "Diabetes")
colnames(data)
dfProts[["425_58"]] <- data


## 432_287 ####
data <- dfProts$`432_287`
colnames(data)
colnames(data) <- c(paste0("Control_", colnames(data)[1:10]), 
                    paste0("T2DM_", colnames(data)[11:ncol(data)]))
colnames(data)
dfProts[["432_287"]] <- data


## 440_211 ####
data <- dfProts$`440_211`
View(data)
colnames(data) # ok


## 511_252 ####
data <- dfProts$`511_252`
View(data)
colnames(data) # ok

# Protein IDs
dfList <- readRDS(file = "../../Data/4_Metanalisis/dfs_quantitative.rds")
df <- dfList$complete %>% 
  filter(ID == "511_252") %>% 
  dplyr::select(ProteinID, ProteinName)
data <- merge(data, df, by = "ProteinName", all = T)

data <- data %>% 
  filter(!is.na(ProteinID)) %>%
  dplyr::select(ProteinID, everything()) %>%
  dplyr::select(-ProteinName)
rownames(data) <- data$ProteinID
data$ProteinID <- NULL

# Group info
colnames(data)
colnames(data) <- c(paste0("Control_", colnames(data)[1:7]), 
                    paste0("T2DM_", colnames(data)[8:ncol(data)]))
colnames(data) <- gsub(colnames(data), replacement = "", pattern = "_DM")
colnames(data)
dfProts[["511_252"]] <- data



## 231_253 ####
data <- dfProts$`231_253`
# Group info
colnames(data)
colnames(data) <- c(paste0("Control_", colnames(data)[1:12]), 
                    paste0("T2DM_", colnames(data)[13:ncol(data)]))
colnames(data)
dfProts[["231_253"]] <- data


## 446_99 ####
data <- dfProts$`446_99`
colnames(data) # ok

data <- data %>% dplyr::select(Ac24:Dx94)
colnames(data) <- c(paste0("Control_", colnames(data)[1:9]),
                    paste0("T2DM_", colnames(data)[10:ncol(data)]))
colnames(data) # ok
dfProts[["446_99"]] <- data



## 499_339 ####
data <- dfProts$`499_339`
colnames(data) # ok
colnames(data) <- c(paste0("Control_", colnames(data)[1:4]),
                    paste0("T2DM_", colnames(data)[5:ncol(data)]))
colnames(data) # ok
dfProts[["499_339"]] <- data


## 944_401 ####
data <- dfProts$`944_401`
colnames(data) 
originalColnames <- colnames(data) 
colnames(data) <- gsub(x = colnames(data), pattern = "_C", replacement = "")
colnames(data) <- gsub(x = colnames(data), pattern = "_DM", replacement = "")
colnames(data) 
colnames(data) <- c("Control_023", "Control_024", "T2DM_034", "T2DM_057",
                    "Control_086", "T2DM_099", "T2DM_118", "Control_125",
                    "Control_146", "T2DM_197")
originalColnames
colnames(data)
dfProts[["944_401"]] <- data



## 1046_653 ####
data <- dfProts$`1046_653`
colnames(data) # ok
colnames(data) <- c(paste0("Control_", colnames(data)[1:9]),
                    paste0("T2DM_", colnames(data)[10:ncol(data)]))
colnames(data) # ok
dfProts[["1046_653"]] <- data


# SAVING ####
saveRDS(dfProts, file = "../../Data/3_Extraction/rawData/curatedDataframes.rds")




































