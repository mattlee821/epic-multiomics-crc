rm(list = ls())

# enviroment ====
#install.packages("readstata13")
library(readstata13)
library(dplyr)
library(data.table)

# epic main data ====
epic <- read.table("data/epic_phenofile.txt", header = T, sep = "\t")
epic <- epic[,c("Idepic", "Match_Caseset")]

# proteins ====
data <- fread("data/proteomics/metaboprep_release_2023_03_30/filtered_data/epic_proteomics_2023_03_30_Filtered_metabolite_data.txt")
colnames(data)[1] <- "Idepic"
data <- left_join(data, epic, by = "Idepic")
## exclude all individuals without a matched caseset
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts < 2], ] # 4
## exclude all matched casesets with more than two people
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts > 2], ] # 0
## ID
proteins_ID <- data$Idepic

# metabolites ====
data <- fread("data/metabolomics/positive/metaboprep_release_2023_03_30/filtered_data/metabolites-positive_2023_03_30_Filtered_metabolite_data.txt")
colnames(data)[1] <- "Idepic"
data <- left_join(data, epic, by = "Idepic")
## exclude all individuals without a matched caseset
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts < 2], ] # 53
## exclude all matched casesets with more than two people
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts > 2], ] # 0
## ID
metabolites_pos_ID <- data$Idepic

# metabolites ====
data <- fread("data/metabolomics/negative/metaboprep_release_2023_03_30/filtered_data/metabolites-negative_2023_03_30_Filtered_metabolite_data.txt")
colnames(data)[1] <- "Idepic"
data <- left_join(data, epic, by = "Idepic")
## exclude all individuals without a matched caseset
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts < 2], ] # 8
## exclude all matched casesets with more than two people
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts > 2], ] # 0
## ID
metabolites_neg_ID <- data$Idepic

# shared ID ====
shared_ID <- intersect(proteins_ID, metabolites_pos_ID)
shared_ID <- intersect(shared_ID, metabolites_neg_ID)

# make QC'd files ====
epic <- read.table("data/epic_phenofile.txt", header = T, sep = "\t")
epic <- epic[epic$Idepic %in% shared_ID, ]
write.table(epic, "data/epic_phenofile_processed.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## proteins 
data <- fread("data/proteomics/metaboprep_release_2023_03_30/filtered_data/epic_proteomics_2023_03_30_Filtered_metabolite_data.txt")
colnames(data)[1] <- "Idepic"
data <- data[data$Idepic %in% shared_ID, ]
data <- left_join(data, epic, by = "Idepic")
write.table(data, "data/proteomics/proteins_phenofile_processed.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## metabolites pos 
data <- fread("data/metabolomics/positive/metaboprep_release_2023_03_30/filtered_data/metabolites-positive_2023_03_30_Filtered_metabolite_data.txt")
colnames(data)[1] <- "Idepic"
data <- data[data$Idepic %in% shared_ID, ]
data <- left_join(data, epic, by = "Idepic")
write.table(data, "data/metabolomics/positive/metabolites-positive_phenofile_processed.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## metabolites neg 
data <- fread("data/metabolomics/negative/metaboprep_release_2023_03_30/filtered_data/metabolites-negative_2023_03_30_Filtered_metabolite_data.txt")
colnames(data)[1] <- "Idepic"
data <- data[data$Idepic %in% shared_ID, ]
data <- left_join(data, epic, by = "Idepic")
write.table(data, "data/metabolomics/negative/metabolites-negative_phenofile_processed.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

