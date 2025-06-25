########################################################-

############   POTENTIAL RELEVANT STUDIES    ###########

########################################################-

# Julia G Currás - 14/02/2025

rm(list=ls())
graphics.off()
search_directory <- ""
setwd(search_directory)
library(dplyr)



# Manually ####
## Open data ####
scopus <-  read.csv(file = paste0(search_directory, "/1_scopus.csv"), header = T)
pubmed <- read.csv(file = paste0(search_directory, "/2_Pubmed.csv"), header = T)
pubmed_2 <- read.csv(file = paste0(search_directory, "/2_PubMed_FromPubmed_format.csv"), 
                   header = T, check.names = F) # From Pubmedtocsv shiny app https://jgcurras.shinyapps.io/pubmedtocsv/
wos <- read.csv(file = paste0(search_directory, "/3_WOS.csv"), header = T, sep = ";")

### Extract abstract from pubmed ####
# pubmed_abstract <- read.table(file = paste0(search_directory, "2_PubMed_abstract.txt"))


## Select common data from all dataframes ####
sort(colnames(scopus))
sort(colnames(pubmed))
sort(colnames(wos))

commonCols <- c(
  "ID", # nº_database
  "Authors", # Authors
  "Title", # Title, Title, Article.Title
  "Abstract", # extract abstract from 2_PubMed_Pubmed.txt or 2_Pubmed_abstract.txt
  "Keywords", # Author.Keywords, - , Author.Keywords
  "Affiliations", # Affiliations, - , Affiliations # extract AD tag from 2_PubMed_Pubmed.txt for all articles
  "Year", # Year, Publication.Year, Publication.Year 
  "Journal", # Source, Journal.Book, Source.Title 
  "DOI"
  )




## Extract the columns ####
### SCOUPS ####
dfScopus <- scopus %>% 
  dplyr::select(Authors, Title, Abstract, Author.Keywords, 
                Affiliations, Year, Source, DOI) %>%
  dplyr::rename(
    Keywords = Author.Keywords,
    Journal = Source) %>%
  dplyr::mutate(
    ID = row_number(), 
    Database = "Scopus"
    ) %>%
  dplyr::select(ID, everything())


### PUBMED ####
dfPubmed <- pubmed %>%
  dplyr::select(Authors, Title, Publication.Year, Journal.Book, 
                DOI, PMID) %>%
  dplyr::rename(
    Year = Publication.Year,
    Journal =  Journal.Book
  ) %>%
  dplyr::mutate(
    ID = row_number(),
    Database = "PubMed") %>%
  dplyr::select(ID, everything())

dfAbstract <- pubmed_2 %>% 
  select(`PubMed Identifier  (PMID)`, 
         `Abstract  (AB)`, `MeSH Terms  (MH)`, 
         `Affiliation  (AD)`) %>% 
  rename(
    PMID = `PubMed Identifier  (PMID)`, 
    Abstract = `Abstract  (AB)`, 
    Keywords = `MeSH Terms  (MH)`,
    Affiliations = `Affiliation  (AD)`)
dfPubmed <- merge(dfAbstract, dfPubmed, by = "PMID", all = T)
dfPubmed <- dfPubmed %>% 
  select(-PMID, ID, Authors, Title, Abstract, )


### WOS ####
dfWOS <- wos %>%
  dplyr::select(Authors, Article.Title, Abstract, Author.Keywords, Affiliations,
                Publication.Year, Source.Title, DOI) %>%
  dplyr::rename(
    Title = Article.Title, 
    Keywords = Author.Keywords,
    Year = Publication.Year,
    Journal = Source.Title
  ) %>%
  dplyr::mutate(
    ID = row_number(),
    Database = "WOS"
    ) %>%
  dplyr::select(ID, everything())




## Merge everything ####
data <- merge(dfScopus, dfWOS, by = colnames(dfScopus), all =T)
data <- merge(data, dfPubmed, by = colnames(dfPubmed), all = T)
data <- data %>%
  # dplyr::arrange(Database) %>%
  mutate(ID = paste0(ID, "_", Database)) %>%
  dplyr::select(ID, everything()) %>%
  dplyr::arrange(Database)
data[data == ""] <- NA # Just in case




## Find duplicated references ####
### With DOI ####
table(table(data$DOI, useNA = 'always') > 1, useNA = 'always') # 657 únicos (28 son NA en DOI), 611 repetidos 
sum(is.na(data$DOI)) # 28
# Venn diagram
library(VennDiagram)
vectorVenn <- list(
  Scopus = data[which(data$Database == "Scopus"), "DOI"],
  PubMed = data[which(data$Database == "PubMed"), "DOI"],
  WOS = data[which(data$Database == "WOS"), "DOI"]
)

vectorVenn <- lapply(vectorVenn, na.omit) # se eliminan 56
venn.diagram(vectorVenn, filename = "doi.png", # Hay valores perdidos en DOI!
             fill = Biostatech::colorPalette(n=3), # 505 en las 3 bases de datos
             cat.cex = 1.5)


### With Title ####
dfAux <- data %>%
  filter(is.na(DOI))
table(table(dfAux$Title) > 1)
table(table(data$Title) > 1)

table(table(data$Title) > 1, useNA = 'always') # 960 unique + 554 duplicates/triplicates

vectorVenn <- list(
  Scopus = data[which(data$Database == "Scopus"), "Title"],
  PubMed = data[which(data$Database == "PubMed"), "Title"],
  WOS = data[which(data$Database == "WOS"), "Title"]
)

vectorVenn <- lapply(vectorVenn, na.omit)
venn.diagram(vectorVenn, filename = "title.png", 
             fill = Biostatech::colorPalette(n=9)[c(1, 5,9)],
             cat.cex = 1.5) # 353 en las 3 bases de datos
sum(is.na(data$Title))


### Removing... ####
dataFinal <- data %>%
  dplyr::distinct(DOI, .keep_all = TRUE) # 1273 potentially relevant studies
# 28 entries did not have DOI information: most of them were reviews, an the others,
# book chapters; so these entries were excluded.

dataFinal <- dataFinal %>% filter(!is.na(DOI))




## SAVING ####
xlsx::write.xlsx(dataFinal, file = "finalSearchDocuments.xlsx", 
                 col.names = T, row.names = T, showNA = F)

































