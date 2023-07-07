# script to identify the 'best' model
rm(list=ls())
set.seed(821)

# environment ====
# install_github("mattlee821/functions", force = T)
library(functions)
library(data.table)
library(dplyr)
library(AICcmodavg)
library(survival)
library(purrr)

# proteins ====
data <- fread("data/proteomics/proteins_phenofile_processed.txt")
# columns of interest ====
columns_main <- c("Idepic", "Match_Caseset", "Cncr_Caco_Clrt")
proteins <- colnames(select(data, contains("olk")))

colnames(select(data, contains("olk_il"), -contains("bio")))
bmi <- c("Idepic", "Bmi_C") # BMI adjusted for clothing, BMI computed - (kg/m2, continuous) no info on what computed means
alcohol <- c("Idepic", "Alc_Re") # alcohol use at recruitment (g/d, continuous)
smoking <- c("Idepic", "Smoke_Stat") # smoking status (never, former, current, unknown) - info from Rhea 
education <- c("Idepic", "L_School") # highest school level (1-5) - no info on what 1-5 is
physicsal_activity <- c("Idepic", "Pa_Index") # Cambridge physical activity index (inactive, moderately inactive, moderately active, active, unknown) - info from Rhea
diet_energy <- c("Idepic", "QE_ENERGY") # energy intake (kcal/d, continuous)
fibre <- c("Idepic", "QE_FIBT") # DQ ENDB Total dietary fibre (g) - from Rhea
meat_red <- c("Idepic", "QgE0701") # red meat consumption (g/d, continuous) - from Rhea
meat_processed <- c("Idepic", "QgE0704") # processed meat consumption (g/d, continuous) - from Rhea
covariates <- c(bmi, alcohol, smoking, education, physicsal_activity, diet_energy, fibre, meat_red, meat_processed)
covariates <- covariates[!covariates %in% "Idepic"]

rm(list = setdiff(ls(), c("data", "columns_main", "proteins", "covariates")))

# data for analysis ====
columns <- c(covariates, proteins)
data_analysis <- data %>% 
  select(Idepic, Match_Caseset, Cncr_Caco_Clrt, 
         any_of(columns))

# model test ====
## randomly select 5 proteins and test for the best model based on provided covariates
proteins_random <- sample(proteins, size = 5)

model_test_results <- list()

system.time(
  for (i in 1:length(proteins_random)){
    
    model_test_results[[i]] <- model_test(data = data_analysis,
                                          exposure = proteins_random[[i]],
                                          outcome = "Cncr_Caco_Clrt",
                                          covariates = covariates,
                                          match_ID = "Match_Caseset",
                                          model = "clogit")
  }
)

## pull out the 5 rows with the lowest values in each dataframe
lowest_rows <- map(model_test_results, function(df) {
  slice(arrange(df, AICc), 1:5)
})
lowest_rows <- bind_rows(lowest_rows)

## remove `exposure ~ outcome`, `+`, and `strata(Match_Caseset)` from column
lowest_rows$covariates <- sub("^.*?~\\s", "", lowest_rows$Modnames)
lowest_rows$covariates <- sub("^.*?\\s", "", lowest_rows$covariates)
lowest_rows$covariates <- sub("^.*?\\s", "", lowest_rows$covariates)
lowest_rows$covariates <- as.factor(lowest_rows$covariates)

## identify the combination of covariates which appear most often across the X proteins tested
freq_table <- table(lowest_rows$covariates)
sorted_levels <- sort(freq_table, decreasing = TRUE)
proteins_largest_factor <- names(sorted_levels)[1] # "Bmi_C + Alc_Re + Smoke_Stat + L_School + Pa_Index + QE_ENERGY + QgE0701 + strata(Match_Caseset)"

# metabolites positive ====
data <- fread("data/metabolomics/positive/metabolites-positive_phenofile_processed.txt")
# columns of interest ====
columns_main <- c("Idepic", "Match_Caseset", "Cncr_Caco_Clrt")
metabolites <- colnames(select(data, contains("untg")))

colnames(select(data, contains("untg")))
bmi <- c("Idepic", "Bmi_C") # BMI adjusted for clothing, BMI computed - (kg/m2, continuous) no info on what computed means
alcohol <- c("Idepic", "Alc_Re") # alcohol use at recruitment (g/d, continuous)
smoking <- c("Idepic", "Smoke_Stat") # smoking status (never, former, current, unknown) - info from Rhea 
education <- c("Idepic", "L_School") # highest school level (1-5) - no info on what 1-5 is
physicsal_activity <- c("Idepic", "Pa_Index") # Cambridge physical activity index (inactive, moderately inactive, moderately active, active, unknown) - info from Rhea
diet_energy <- c("Idepic", "QE_ENERGY") # energy intake (kcal/d, continuous)
fibre <- c("Idepic", "QE_FIBT") # DQ ENDB Total dietary fibre (g) - from Rhea
meat_red <- c("Idepic", "QgE0701") # red meat consumption (g/d, continuous) - from Rhea
meat_processed <- c("Idepic", "QgE0704") # processed meat consumption (g/d, continuous) - from Rhea
covariates <- c(bmi, alcohol, smoking, education, physicsal_activity, diet_energy, fibre, meat_red, meat_processed)
covariates <- covariates[!covariates %in% "Idepic"]

# data for analysis ====
columns <- c(covariates, metabolites)
data_analysis <- data %>% 
  select(Idepic, Match_Caseset, Cncr_Caco_Clrt, 
         any_of(columns))

# model test ====
## randomly select 5 metabolites and test for the best model based on provided covariates
metabolites_random <- sample(metabolites, size = 5)

model_test_results <- list()

system.time(
  for (i in 1:length(metabolites_random)){
    
    model_test_results[[i]] <- model_test(data = data_analysis,
                                          exposure = metabolites_random[[i]],
                                          outcome = "Cncr_Caco_Clrt",
                                          covariates = covariates,
                                          match_ID = "Match_Caseset",
                                          model = "clogit")
  }
)

## pull out the 5 rows with the lowest values in each dataframe
lowest_rows <- map(model_test_results, function(df) {
  slice(arrange(df, AICc), 1:5)
})
lowest_rows <- bind_rows(lowest_rows)

# remove `exposure ~ outcome`, `+`, and `strata(Match_Caseset)` from column
lowest_rows$covariates <- sub("^.*?~\\s", "", lowest_rows$Modnames)
lowest_rows$covariates <- sub("^.*?\\s", "", lowest_rows$covariates)
lowest_rows$covariates <- sub("^.*?\\s", "", lowest_rows$covariates)
lowest_rows$covariates <- as.factor(lowest_rows$covariates)

# identify the combination of covariates which appear most often across the X metabolites tested
freq_table <- table(lowest_rows$covariates)
sorted_levels <- sort(freq_table, decreasing = TRUE)
metabolites_positive_largest_factor <- names(sorted_levels)[1] # "Bmi_C + Alc_Re + L_School + Pa_Index + QgE0701 + strata(Match_Caseset)"


# metabolites negative ====
data <- fread("data/metabolomics/negative/metabolites-negative_phenofile_processed.txt")
# columns of interest ====
columns_main <- c("Idepic", "Match_Caseset", "Cncr_Caco_Clrt")
metabolites <- colnames(select(data, contains("untg")))

colnames(select(data, contains("untg")))
bmi <- c("Idepic", "Bmi_C") # BMI adjusted for clothing, BMI computed - (kg/m2, continuous) no info on what computed means
alcohol <- c("Idepic", "Alc_Re") # alcohol use at recruitment (g/d, continuous)
smoking <- c("Idepic", "Smoke_Stat") # smoking status (never, former, current, unknown) - info from Rhea 
education <- c("Idepic", "L_School") # highest school level (1-5) - no info on what 1-5 is
physicsal_activity <- c("Idepic", "Pa_Index") # Cambridge physical activity index (inactive, moderately inactive, moderately active, active, unknown) - info from Rhea
diet_energy <- c("Idepic", "QE_ENERGY") # energy intake (kcal/d, continuous)
fibre <- c("Idepic", "QE_FIBT") # DQ ENDB Total dietary fibre (g) - from Rhea
meat_red <- c("Idepic", "QgE0701") # red meat consumption (g/d, continuous) - from Rhea
meat_processed <- c("Idepic", "QgE0704") # processed meat consumption (g/d, continuous) - from Rhea
covariates <- c(bmi, alcohol, smoking, education, physicsal_activity, diet_energy, fibre, meat_red, meat_processed)
covariates <- covariates[!covariates %in% "Idepic"]

# data for analysis ====
columns <- c(covariates, metabolites)
data_analysis <- data %>% 
  select(Idepic, Match_Caseset, Cncr_Caco_Clrt, 
         any_of(columns))

# model test ====
## randomly select 5 metabolites and test for the best model based on provided covariates
metabolites_random <- sample(metabolites, size = 5)

model_test_results <- list()

system.time(
  for (i in 1:length(metabolites_random)){
    
    model_test_results[[i]] <- model_test(data = data_analysis,
                                          exposure = metabolites_random[[i]],
                                          outcome = "Cncr_Caco_Clrt",
                                          covariates = covariates,
                                          match_ID = "Match_Caseset",
                                          model = "clogit")
  }
)

## pull out the 5 rows with the lowest values in each dataframe
lowest_rows <- map(model_test_results, function(df) {
  slice(arrange(df, AICc), 1:5)
})
lowest_rows <- bind_rows(lowest_rows)

# remove `exposure ~ outcome`, `+`, and `strata(Match_Caseset)` from column
lowest_rows$covariates <- sub("^.*?~\\s", "", lowest_rows$Modnames)
lowest_rows$covariates <- sub("^.*?\\s", "", lowest_rows$covariates)
lowest_rows$covariates <- sub("^.*?\\s", "", lowest_rows$covariates)
lowest_rows$covariates <- as.factor(lowest_rows$covariates)

# identify the combination of covariates which appear most often across the X metabolites tested
freq_table <- table(lowest_rows$covariates)
sorted_levels <- sort(freq_table, decreasing = TRUE)
metabolites_negative_largest_factor <- names(sorted_levels)[1] # "Bmi_C + Alc_Re + L_School + Pa_Index + QgE0701 + strata(Match_Caseset)"


# save ====
table <- data.frame(proteins = proteins_largest_factor, 
                    metabolites_positive = metabolites_positive_largest_factor, 
                    metabolites_negative = metabolites_negative_largest_factor)
write.table(table, "data/aicc-model.txt")
