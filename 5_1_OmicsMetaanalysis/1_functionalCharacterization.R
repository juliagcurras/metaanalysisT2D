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
                        sheet = "Qualitative (protein>1 study)", col_names = T)
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
xlsx::write.xlsx(x = as.data.frame(int_test), 
                 col.names = T, row.names = F, append = F,
                 file = "results_Qualitative.xlsx", 
                 sheetName = "Interaction pairs (STRING)")



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
xlsx::write.xlsx(x = as.data.frame(resultGO@result), 
                 col.names = T, row.names = F, append = T,
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
xlsx::write.xlsx(x = as.data.frame(resultKEGG@result), 
                 col.names = T, row.names = F, append = T,
                 file = "results_Qualitative.xlsx", 
                 sheetName = "ORA KEGG (proteins in >7 studies)")








## Semantic similarity analysis ###
resultGO_original <- resultGO

replicateClustering <- function(resultGO, k = 5, nwords = 2,
                                metodo = "ward.D"){
  #### Replicate clustering from treeplot using the similarity matrix computed in pairwise()
  
  # .- clustering -.#
  # Retrieve similarity matrix
  sim <- resultGO@termsim
  # Convert to dissimilarity/distance matrix
  d <- as.dist(1 - sim)
  # Hierarchical clustering (same as treeplot())
  hc <- hclust(d, method = metodo)
  # Cortar el dendrograma
  clusters <- cutree(hc, k = k)
  # Crear el dataframe
  df_clusters <- data.frame(
    Description = names(clusters),
    Cluster = clusters
  )
  # Adding GO ID
  df_clusters <- merge(
    df_clusters,
    resultGO@result[, c("ID", "Description")],
    by = "Description"
  )
  
  # .- Naming clusters -.#
  # Step 1: remove empty words and tokenin descriptions
  cluster_terms <- df_clusters %>%
    unnest_tokens(word, Description) %>%
    filter(!word %in% stop_words$word,  # remove common words
           !str_detect(word, "^\\d+$")) # remove numbers
  # Step 2: count frequency of words by cluster
  word_freq <- cluster_terms %>%
    count(Cluster, word, sort = TRUE)
  # Step 3: most frequent word by cluster
  top_words <- word_freq %>%
    group_by(Cluster) %>%
    top_n(nwords, n) %>% 
    arrange(Cluster, desc(n)) %>%
    summarise(ClusterLabel = paste(word[1:nwords], collapse = " "), .groups = "drop")
  
  # Step 4: Joint words for a final label
  df_clusters <- df_clusters %>%
    left_join(top_words, by = "Cluster") %>%
    dplyr::select(ID, Description, Cluster, ClusterLabel) %>% 
    arrange(Cluster)
  word_freq <- word_freq %>% arrange(Cluster)
  
  # Return
  return(list(
    clusters = df_clusters,
    wordCount = word_freq))
}



### CC ####
colors <- c("#AB83CF", "#73C4C8", "#00437B", "#342C5F", "#117BA0")
# Automatically with enrichplor
resultGO@result <- resultGO@result %>% 
  filter(ONTOLOGY == "CC") %>%
  filter(qvalue < 0.05) %>%
  filter(p.adjust < 0.05)
resultGO <- enrichplot::pairwise_termsim(resultGO)
p1 <- enrichplot::treeplot(resultGO, nCluster = 5, nWords =3, 
                     showCategory = nrow(resultGO@result), 
                     group_color = colors, hexpand = 0.5)
ggsave(filename = "CC_treeplot.svg", 
       plot = p1, width = 15, height = 12)
# Retrieving table with clustering with enrichplot
resClus <- replicateClustering(resultGO, k = 7)
# Saving
xlsx::write.xlsx(x = as.data.frame(resClus$clusters), 
                 file = "results_Qualitative.xlsx",
                 sheetName = "CC clustering", col.names = T, row.names = F, 
                 append = T, showNA = F
)
xlsx::write.xlsx(x = resClus$wordCount, 
                 file = "results_Qualitative.xlsx",
                 sheetName = "CC word count", col.names = T, row.names = F, 
                 append = T, showNA = F
)



### MF ####
resultGO <- resultGO_original
colors <- c("#675BA5", "#C992DC", "#312959", "#096C9D", "#73C4C8")
# Automatically with enrichplor
resultGO@result <- resultGO@result %>% 
  filter(ONTOLOGY == "MF") %>% 
  filter(qvalue < 0.05) %>%
  filter(p.adjust < 0.05)
resultGO <- enrichplot::pairwise_termsim(resultGO)
p2 <- enrichplot::treeplot(resultGO, nCluster = 5, nWords =3, 
                           showCategory = nrow(resultGO@result), 
                           group_color = colors, hexpand = 0.5)
ggsave(filename = "MF_treeplot.svg", 
       plot = p2, width = 15, height = 12)
# Retrieving table with clustering with enrichplot
resClus <- replicateClustering(resultGO, k = 5)
# Saving
xlsx::write.xlsx(x = as.data.frame(resClus$clusters), 
                 file = "results_Qualitative.xlsx",
                 sheetName = "MF clustering", col.names = T, row.names = F, 
                 append = T, showNA = F
)
xlsx::write.xlsx(x = resClus$wordCount, 
                 file = "results_Qualitative.xlsx",
                 sheetName = "MF word count", col.names = T, row.names = F, 
                 append = T, showNA = F
)


### BP ####
resultGO <- resultGO_original 
colors <- c("#004B85", "#8CCBD2", "#4B3F87", "#004780", "#C7C8E5")
# Automatically with enrichplor
resultGO@result <- resultGO@result %>% 
  filter(ONTOLOGY == "BP") %>% 
  filter(p.adjust < 0.05)%>%
  filter(qvalue < 0.05)
resultGO <- clusterProfiler::simplify(resultGO)
resultGO <- enrichplot::pairwise_termsim(resultGO)
p3 <- enrichplot::treeplot(resultGO, nCluster = 5, nWords =3, 
                           showCategory = nrow(resultGO@result),
                           group_color = colors, hexpand = 0.5)
ggsave(filename = "BP_treeplot.svg", 
       plot = p3, width = 20, height = 25)
# Retrieving table with clustering with enrichplot
resClus <- replicateClustering(resultGO, k = 5, nwords = 4)
xlsx::write.xlsx(x = as.data.frame(resClus$clusters), 
                 file = "results_Qualitative.xlsx",
                 sheetName = "BP clustering", col.names = T, row.names = F, 
                 append = T, showNA = F
)
xlsx::write.xlsx(x = resClus$wordCount, 
                 file = "results_Qualitative.xlsx",
                 sheetName = "BP word count", col.names = T, row.names = F, 
                 append = T, showNA = F
)






