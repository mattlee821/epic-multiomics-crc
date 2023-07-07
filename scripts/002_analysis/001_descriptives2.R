rm(list=ls())
set.seed(821)

# environment ====
library(functions)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)

# proteomics ====
data <- fread("data/proteomics/proteins_phenofile_processed.txt")
columns <- colnames(select(data, contains("olk")))
proteomics <- data %>% 
  select(Idepic, Cncr_Caco_Clrt, Sex,
         any_of(columns))

# metabolomics positive ====
data <- fread("data/metabolomics/positive/metabolites-positive_phenofile_processed.txt")
columns <- colnames(select(data, contains("untg")))
metabolomics_pos <- data %>% 
  select(Idepic, 
         any_of(columns))

# metabolomics negative ====
data <- fread("data/metabolomics/negative/metabolites-negative_phenofile_processed.txt")
columns <- colnames(select(data, contains("untg")))
metabolomics_neg <- data %>% 
  select(Idepic, 
         any_of(columns))

# combine ====
data <- left_join(proteomics,metabolomics_pos, by = "Idepic")
data <- left_join(data,metabolomics_neg, by = "Idepic")
data <- data[,-1]
rm(list = setdiff(ls(), "data"))

# table ====
## combined
summary_combined <- data %>%
  summarize(across(where(is.numeric), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))))
summary_combined <- summary_combined[,-c(1:4)]
a <- data.frame(Cncr_Caco_Clrt = "combined",
                Sex = "combined")
summary_combined <- cbind(a,summary_combined)

## cancer
summary_cancer <- data %>%
  group_by(Cncr_Caco_Clrt) %>%
  summarize(across(where(is.numeric), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))), .groups = "drop")
summary_cancer <- summary_cancer[,-3]
colnames(summary_cancer)[2] <- "Sex"
summary_cancer$Sex <- "combined"
summary_cancer <- summary_cancer %>%
  mutate(Cncr_Caco_Clrt = ifelse(Cncr_Caco_Clrt == 0, "control", "case"))

## sex
summary_sex <- data %>%
  group_by(Sex) %>%
  summarize(across(where(is.numeric), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))), .groups = "drop")
summary_sex <- summary_sex[,-3]
colnames(summary_sex)[2] <- "Cncr_Caco_Clrt"
summary_sex$Cncr_Caco_Clrt <- "combined"
summary_sex <- summary_sex[,c(2,1,3:ncol(summary_sex))]
summary_sex <- summary_sex %>%
  mutate(Sex = ifelse(Sex == 1, "male", "female"))

## cancer and sex 
summary_cancer_sex <- data %>%
  group_by(Cncr_Caco_Clrt, Sex) %>%
  summarize(across(where(is.numeric), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))), .groups = "drop")

summary_cancer_sex <- summary_cancer_sex %>%
  mutate(Cncr_Caco_Clrt = ifelse(Cncr_Caco_Clrt == 0, "control", "case")) %>%
  mutate(Sex = ifelse(Sex == 1, "male", "female"))

## combine 
summary <- rbind(summary_combined, summary_cancer, summary_sex, summary_cancer_sex)
summary <- as.data.frame(t(summary))
colnames(summary) <- c("combined", "control", "case", "male", "female", "male_control", "female_control", "male_case", "female_case")
summary <- select(summary, combined, control, case, female, female_control, female_case, male, male_control, male_case)
summary <- summary[-c(1,2),]
summary <- rownames_to_column(summary, var = "trait")

write.table(summary, "analysis/tables/descriptives_metabolomics_proteomics.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
