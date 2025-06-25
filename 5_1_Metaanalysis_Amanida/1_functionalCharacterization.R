##############################################################################-

#############   Funcitonal characterization of biomarkers    #################-

##############################################################################-


# Julia G Currás - 2025/05/21

# Setup ####
rm(list=ls())
graphics.off()
pathTofigures <- ""
pathToData <- ""
library(dplyr)
library(rbioapi)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GOSemSim)
library(stringr)
library(factoextra)



# Qualitative data ####
## Loading data ####
df <- readxl::read_xlsx(path = "../4_Processing/allProcessedProteinData.xlsx", 
                        sheet = "Qualitative (>2 studies)", col_names = T)
allFrec <- table(df$ProteinID)
sum(table(allFrec))

## Saving interesting proteins (identified > 7 studies) ####
prots <- names(allFrec[allFrec > 7])

## Interaction analysis - STRING ####
# Map #
targetID <- rbioapi::rba_string_map_ids(prots, species = 9606)
# Net # 
rbioapi::rba_string_network_image(ids = targetID$stringId, 
                                  required_score = 700,
                                  image_format = "image",
                                  node_labels_font_size = 20, 
                                  save_image = file.path("STRING_qualitative.png"),
                                  species = 9606) # only query proteins

int_net <- rba_string_interactions_network(ids = targetID$stringId, # todas
                                           species = 9606,
                                           required_score = 700) # medium confident
int_test <- rba_string_enrichment_ppi(ids = targetID$stringId, species = 9606, 
                                      required_score = 700)

output <- list(tableNet = int_net, tableTest = int_test)


## Functional analysis ###
### GO ORA analysis ####
#### With universe ####
resultGO <- enrichGO(
  gene = prots,
  OrgDb = "org.Hs.eg.db",
  keyType = "UNIPROT",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = unique(df$ProteinID),
  qvalueCutoff = 0.2,
  minGSSize = 3,
  maxGSSize = 800,
  readable = FALSE,
  pool = FALSE
)

resultGO@result %>% View
resultGO@result <- resultGO@result %>% arrange(FoldEnrichment)
# saveRDS(resultGO, file = paste0(pathToData, "functionalAnalysis/GO_ORA_qualitative.rds"))
xlsx::write.xlsx(x = as.data.frame(resultGO@result), 
                 col.names = T, row.names = F, append = F,
                 file = "results_Qualitative.xlsx", 
                 sheetName = "ORA GO (proteins in >7 studies)")


### KEGG ORA analysis ####
#### With universe ####
resultKEGG <- enrichKEGG(
  gene = prots,
  organism = "hsa",
  keyType = "uniprot",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe =  unique(df$ProteinID),
  minGSSize = 3,
  maxGSSize = 800,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
resultKEGG@result %>% View
# saveRDS(resultKEGG, file = paste0(pathToData, "functionalAnalysis/KEGG_ORA_qualitative.rds"))
xlsx::write.xlsx(x = as.data.frame(resultKEGG@result), 
                 col.names = T, row.names = F, append = T,
                 file = "results_Qualitative.xlsx", 
                 sheetName = "ORA KEGG (proteins in >7 studies)")








## Semantic similarity analysis ###
resultGO_original <- resultGO
### CC ####
colors <- c("#AB83CF", "#73C4C8", "#00437B", "#342C5F", "#117BA0")
# Automatically with enrichplor
resultGO@result <- resultGO@result %>% 
  filter(ONTOLOGY == "CC") %>% 
  filter(p.adjust < 0.05)
semsim <- enrichplot::pairwise_termsim(resultGO)
enrichplot::treeplot(semsim, nCluster = 5, nWords =3, 
                     showCategory = nrow(resultGO@result), 
                     group_color = colors, hexpand = 0.5)


### MF ####
resultGO <- resultGO_original
colors <- c("#675BA5", "#C992DC", "#312959", "#096C9D", "#73C4C8")
# Automatically with enrichplor
resultGO@result <- resultGO@result %>% 
  filter(ONTOLOGY == "MF") %>% 
  filter(p.adjust < 0.05)
semsim <- enrichplot::pairwise_termsim(resultGO)
enrichplot::treeplot(semsim, nCluster = 5, nWords =3, 
                     showCategory = nrow(resultGO@result), 
                     group_color = colors, hexpand = 0.5)


### BP ####
resultGO <- resultGO_original 
colors <- c("#004B85", "#096C9D", "#4B3F87", "#004780", "#C7C8E5", "#249FA8",
                "#C3CEE5", "#06669C", "#342C5F", "#8CCBD2")
# Automatically with enrichplor
resultGO@result <- resultGO@result %>% 
  filter(ONTOLOGY == "BP") %>% 
  filter(p.adjust < 0.05)
semsim <- enrichplot::pairwise_termsim(resultGO)
enrichplot::treeplot(semsim, nCluster = 10, nWords =3, 
                     # showCategory = nrow(resultGO@result), 
                     group_color = colors, hexpand = 0.5)






