rm(list=ls())

# enviroment ====
library(haven)
library(readstata13)
library(dplyr)

# data ====
## epic main data
epic <- read.dta13("../000_datasets/EPIC/clrt_caco.dta")
epic_ID <- epic$Idepic
  
## proteomics
proteins <- read.dta13("../000_datasets/EPIC/proteomics/proteomics_olink.dta")
proteins <- select(proteins, -Country, -Center, -Idepic_Bio, -Plate_Id, -Sample_Typ)
proteins_ID <- proteins$Idepic

## metabolomics
metabolites_pos <- read_sas("../000_datasets/EPIC/metabolomics/untg_rp_pos_metabo_clrt_2023_feb.sas7bdat")
metabolites_pos <- select(metabolites_pos, -Batch, -Idepic_Bio, -Idepic_Samp, -Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt)
metabolites_pos_ID <- metabolites_pos$Idepic

metabolites_neg <- read_sas("../000_datasets/EPIC/metabolomics/untg_rp_neg_metabo_clrt_2023_feb.sas7bdat")
metabolites_neg <- select(metabolites_neg, -Batch, -Idepic_Bio, -Idepic_Samp, -Match_Caseset, -Country, -Center, -Cncr_Caco_Clrt)
metabolites_neg_ID <- metabolites_neg$Idepic

## genetic data


## shared ID
shared_ID <- intersect(epic_ID, proteins_ID)
shared_ID <- intersect(shared_ID, metabolites_pos_ID)
shared_ID <- intersect(shared_ID, metabolites_neg_ID)

# make combined dataframe ====
data <- left_join(epic, proteins, by = "Idepic")
data <- left_join(data, metabolites_pos, by = "Idepic")
data <- left_join(data, metabolites_neg, by = "Idepic")
## format 
data <- data[data$Idepic %in% shared_ID, ]
## exclude inelligible tumours
data <- subset(data, Typ_Tumo != "Exc-Morphology excluded") # 11
data <- subset(data, Typ_Tumo != "Exc-Morphology ineligible") # 2
## exclude all individuals without a matched caseset
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts < 2], ] # 30
## exclude all matched casesets with more than two people
level_counts <- table(data$Match_Caseset)
data <- data[!data$Match_Caseset %in% names(level_counts)[level_counts > 2], ] # 14
## ID
shared_ID <- data$Idepic

# make pre QC epic file ====
epic <- epic[epic$Idepic %in% shared_ID, ]
## exclude inelligible tumours
epic <- subset(epic, Typ_Tumo != "Exc-Morphology excluded") # 11
epic <- subset(epic, Typ_Tumo != "Exc-Morphology ineligible") # 2
## exclude all individuals without a matched caseset
level_counts <- table(epic$Match_Caseset)
epic <- epic[!epic$Match_Caseset %in% names(level_counts)[level_counts < 2], ] # 30
## exclude all matched casesets with more than two people
level_counts <- table(epic$Match_Caseset)
epic <- epic[!epic$Match_Caseset %in% names(level_counts)[level_counts > 2], ] # 14
## save
write.table(epic, "data/epic_phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make pre QC protein file ====
proteins <- proteins[proteins$Idepic %in% shared_ID, ]
proteins <- select(proteins, contains(c("idepic","olk")))
proteins <- select(proteins, -contains(c("outdq", "bio")))
proteins <- proteins[!duplicated(proteins),]
write.table(proteins, "data/proteomics/proteins_raw.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make pre QC metab-pos file ====
metabolites_pos <- metabolites_pos[metabolites_pos$Idepic %in% shared_ID, ]
metabolites_pos <- select(metabolites_pos, contains(c("idepic","Untg_Rp")), -contains("bio"))
metabolites_pos <- metabolites_pos[!duplicated(metabolites_pos),]
write.table(metabolites_pos, "data/metabolomics/metabolites-pos_raw.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make pre QC metab-pos file ====
metabolites_neg <- metabolites_neg[metabolites_neg$Idepic %in% shared_ID, ]
metabolites_neg <- select(metabolites_neg, contains(c("idepic","Untg_Rp")), -contains("bio"))
metabolites_neg <- metabolites_neg[!duplicated(metabolites_neg),] 
write.table(metabolites_neg, "data/metabolomics/metabolites-neg_raw.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
