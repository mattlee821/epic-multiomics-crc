rm(list=ls())
set.seed(821)

# environment ====
library(survival)
library(gtsummary)
library(ggplot2)
library(dplyr)
library(knitr)
library(patchwork)
library(tidyr)
library(purrr)
library(rlang)
library(data.table)

# data ====
data <- read.table("data/metabolomics/negative/metabolites-negative_phenofile_processed.txt", header = T, sep = "\t")
aicc <- read.table("data/aicc-model.txt")
aicc <- aicc[1,2]

# columns of interest ====
columns_main <- c("Idepic", "Match_Caseset", "Cncr_Caco_Clrt")
trait <- colnames(select(data, contains("untg")))

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
columns <- c(covariates, trait)
data_analysis <- data %>% 
  select(Idepic, Match_Caseset, Cncr_Caco_Clrt, Sex, Siteclrt, Length_Bld, 
         any_of(columns))
coln <- ncol(data_analysis) - length(trait) + 1
data_analysis$Length_Bld <- data_analysis$Length_Bld/365
data_analysis <- subset(data_analysis, Length_Bld > 2) # removes 42
level_counts <- table(data_analysis$Match_Caseset)
data_analysis <- data_analysis[!data_analysis$Match_Caseset %in% names(level_counts)[level_counts < 2], ]
data_analysis_male <- subset(data_analysis, Sex == 1) # 16
data_analysis_female <- subset(data_analysis, Sex == 2) # 26
data_analysis_combined <- data_analysis

# site specific ====
proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
distal <- c("C186", "C187")
overlapping <- "C188"
unspecified <- "C189"

data_analysis_proximal <- data_analysis[data_analysis$Siteclrt %in% proximal, ]
data_analysis_proximal <- data_analysis_proximal$Match_Caseset
data_analysis_proximal_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_proximal, ]
data_analysis_proximal_male <- subset(data_analysis_proximal_combined, Sex == 1)
data_analysis_proximal_female <- subset(data_analysis_proximal_combined, Sex == 2)

data_analysis_distal <- data_analysis[data_analysis$Siteclrt %in% distal, ]
data_analysis_distal <- data_analysis_distal$Match_Caseset
data_analysis_distal_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_distal, ]
data_analysis_distal_male <- subset(data_analysis_distal_combined, Sex == 1)
data_analysis_distal_female <- subset(data_analysis_distal_combined, Sex == 2)

data_analysis_overlapping <- data_analysis[data_analysis$Siteclrt %in% overlapping, ]
data_analysis_overlapping <- data_analysis_overlapping$Match_Caseset
data_analysis_overlapping_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_overlapping, ]
data_analysis_overlapping_male <- subset(data_analysis_overlapping_combined, Sex == 1)
data_analysis_overlapping_female <- subset(data_analysis_overlapping_combined, Sex == 2)

data_analysis_unspecified <- data_analysis[data_analysis$Siteclrt %in% unspecified, ]
data_analysis_unspecified <- data_analysis_unspecified$Match_Caseset
data_analysis_unspecified_combined <- data_analysis[data_analysis$Match_Caseset %in% data_analysis_unspecified, ]
data_analysis_unspecified_male <- subset(data_analysis_unspecified_combined, Sex == 1)
data_analysis_unspecified_female <- subset(data_analysis_unspecified_combined, Sex == 2)

# median followup time ====
data_analysis_combined_below <- subset(data_analysis_combined, Length_Bld < median(data_analysis_combined$Length_Bld))
data_analysis_combined_above <- subset(data_analysis_combined, Length_Bld > median(data_analysis_combined$Length_Bld))
data_analysis_proximal_combined_below <- subset(data_analysis_proximal_combined, Length_Bld < median(data_analysis_proximal_combined$Length_Bld))
data_analysis_proximal_combined_above <- subset(data_analysis_proximal_combined, Length_Bld > median(data_analysis_proximal_combined$Length_Bld))
data_analysis_distal_combined_below <- subset(data_analysis_distal_combined, Length_Bld < median(data_analysis_distal_combined$Length_Bld))
data_analysis_distal_combined_above <- subset(data_analysis_distal_combined, Length_Bld > median(data_analysis_distal_combined$Length_Bld))

data_analysis_male_below <- subset(data_analysis_male, Length_Bld < median(data_analysis_male$Length_Bld))
data_analysis_male_above <- subset(data_analysis_male, Length_Bld > median(data_analysis_male$Length_Bld))
data_analysis_proximal_male_below <- subset(data_analysis_proximal_male, Length_Bld < median(data_analysis_proximal_male$Length_Bld))
data_analysis_proximal_male_above <- subset(data_analysis_proximal_male, Length_Bld > median(data_analysis_proximal_male$Length_Bld))
data_analysis_distal_male_below <- subset(data_analysis_distal_male, Length_Bld < median(data_analysis_distal_male$Length_Bld))
data_analysis_distal_male_above <- subset(data_analysis_distal_male, Length_Bld > median(data_analysis_distal_male$Length_Bld))

data_analysis_female_below <- subset(data_analysis_female, Length_Bld < median(data_analysis_female$Length_Bld))
data_analysis_female_above <- subset(data_analysis_female, Length_Bld > median(data_analysis_female$Length_Bld))
data_analysis_proximal_female_below <- subset(data_analysis_proximal_female, Length_Bld < median(data_analysis_proximal_female$Length_Bld))
data_analysis_proximal_female_above <- subset(data_analysis_proximal_female, Length_Bld > median(data_analysis_proximal_female$Length_Bld))
data_analysis_distal_female_below <- subset(data_analysis_distal_female, Length_Bld < median(data_analysis_distal_female$Length_Bld))
data_analysis_distal_female_above <- subset(data_analysis_distal_female, Length_Bld > median(data_analysis_distal_female$Length_Bld))

# sex analysis ====
# sex-combined below ====
results_list <- list()
data_analysis <- data_analysis_combined_below
## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
sex_combined <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
sex_combined$sex <- "sex-combined"
sex_combined$cancer <- "overall"

# male below ====
results_list <- list()
data_analysis <- data_analysis_male_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
male <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
male$sex <- "male"
male$cancer <- "overall"

# female below ====
results_list <- list()
data_analysis <- data_analysis_female_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
female <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
female$sex <- "female"
female$cancer <- "overall"

# combine analyses ====
results_sex_below <- rbind(sex_combined,male,female)
results_sex_below <- results_sex_below[order(results_sex_below$exposure),]
results_sex_below$followup <- "below"

# sex-combined above ====
results_list <- list()
data_analysis <- data_analysis_combined_above
## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
sex_combined <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
sex_combined$sex <- "sex-combined"
sex_combined$cancer <- "overall"

# male above ====
results_list <- list()
data_analysis <- data_analysis_male_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
male <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
male$sex <- "male"
male$cancer <- "overall"

# female above ====
results_list <- list()
data_analysis <- data_analysis_female_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
female <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
female$sex <- "female"
female$cancer <- "overall"

# combine analyses ====
results_sex_above <- rbind(sex_combined,male,female)
results_sex_above <- results_sex_above[order(results_sex_above$exposure),]
results_sex_above$followup <- "above"

# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_combined_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
proximal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
proximal$sex <- "sex-combined"
proximal$cancer <- "proximal"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_combined_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
distal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
distal$sex <- "sex-combined"
distal$cancer <- "distal"

## combine analyses ====
results_combined_below <- rbind(proximal, distal)
results_combined_below <- results_combined_below[order(results_combined_below$exposure),]

# male ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_male_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
proximal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
proximal$sex <- "male"
proximal$cancer <- "proximal"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_male_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
distal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
distal$sex <- "male"
distal$cancer <- "distal"

## combine analyses ====
results_male_below <- rbind(proximal, distal)
results_male_below <- results_male_below[order(results_male_below$exposure),]

# female ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_female_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
proximal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
proximal$sex <- "female"
proximal$cancer <- "proximal"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_female_below

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
distal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
distal$sex <- "female"
distal$cancer <- "distal"

## combine analyses ====
results_female_below <- rbind(proximal, distal)
results_female_below <- results_female_below[order(results_female_below$exposure),]

# combine site below ====
results_site_below <- rbind(results_combined_below, results_male_below, results_female_below)
results_site_below <- results_site_below[order(results_site_below$exposure),]
results_site_below$followup <- "below"

# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_combined_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
proximal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
proximal$sex <- "sex-combined"
proximal$cancer <- "proximal"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_combined_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
distal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
distal$sex <- "sex-combined"
distal$cancer <- "distal"

## combine analyses ====
results_combined_above <- rbind(proximal, distal)
results_combined_above <- results_combined_above[order(results_combined_above$exposure),]

# male ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_male_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
proximal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
proximal$sex <- "male"
proximal$cancer <- "proximal"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_male_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
distal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
distal$sex <- "male"
distal$cancer <- "distal"

## combine analyses ====
results_male_above <- rbind(proximal, distal)
results_male_above <- results_male_above[order(results_male_above$exposure),]

# female ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_female_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
proximal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
proximal$sex <- "female"
proximal$cancer <- "proximal"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_female_above

## model 1 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model1 <- do.call(rbind, results_list)
results_model1$model <-  "model 1"

## model 2 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model2 <- do.call(rbind, results_list)
results_model2$model <-  "model 2"

## model 3 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model3 <- do.call(rbind, results_list)
results_model3$model <-  "model 3"

## model 4 ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + 
                                              Bmi_C + Alc_Re + Smoke_Stat +
                                              L_School + Pa_Index + QE_ENERGY +
                                              QE_FIBT + QgE0701 + QgE0704 +
                                              strata(Match_Caseset), 
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model4 <- do.call(rbind, results_list)
results_model4$model <- "model 4"

## model 5 - ACIc identified covariates ====
for(i in coln:ncol(data_analysis)) { 
  
  exposure <- colnames(data_analysis)[i]
  #print(exposure)
  
  eval(parse(text = paste("temporary <- clogit(Cncr_Caco_Clrt ~ ",exposure," + ", aicc, ",
                                              data = data_analysis, method = 'exact', na.action = 'na.exclude')", sep =""))) 
  #print(temporary)
  
  temp_df <- data.frame(
    exposure = names(temporary[["coefficients"]][1]),
    n = summary(temporary)[["n"]],
    nevent = summary(temporary)[["nevent"]],
    b = summary(temporary)[["coefficients"]][1,1],
    se = summary(temporary)[["coefficients"]][1,3],
    p = summary(temporary)[["coefficients"]][1,5]
  )
  
  j=i-(1-coln)
  results_list[[j]] <- temp_df
}  
results_model5 <- do.call(rbind, results_list)
results_model5$model <- "model 5"

### join models ====
distal <- rbind(results_model1, results_model2, results_model3, results_model4, results_model5)
distal$sex <- "female"
distal$cancer <- "distal"

## combine analyses ====
results_female_above <- rbind(proximal, distal)
results_female_above <- results_female_above[order(results_female_above$exposure),]

# combine site above ====
results_site_above <- rbind(results_combined_above, results_male_above, results_female_above)
results_site_above <- results_site_above[order(results_site_above$exposure),]
results_site_above$followup <- "above"

# combine and save all ====
results <- rbind(results_sex_below, results_sex_above, results_site_below, results_site_above)
results <- results[order(results$exposure),]
results$analysis <- "2year follow-up"
write.table(results, "analysis/002_analysis/metabolomics/negative/002_analysis_sex_site-specific_median-followup_removed-2y-followup.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
