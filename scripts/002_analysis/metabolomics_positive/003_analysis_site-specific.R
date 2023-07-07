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
data <- read.table("data/metabolomics/positive/metabolites-positive_phenofile_processed.txt", header = T, sep = "\t")
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
  select(Idepic, Match_Caseset, Cncr_Caco_Clrt, Sex, Siteclrt,
         any_of(columns))
coln <- ncol(data_analysis) - length(trait) + 1
data_analysis_male <- subset(data_analysis, Sex == 1)
data_analysis_female <- subset(data_analysis, Sex == 2)
data_analysis_combined <- data_analysis

# site-specific data ====
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

# sex-combined ====
results_list <- list()
data_analysis <- data_analysis_combined
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
sex_combined$followup <- "NA"

# sex-specific analysis ====

# male ====
results_list <- list()
data_analysis <- data_analysis_male

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
male$followup <- "NA"

# female ====
results_list <- list()
data_analysis <- data_analysis_female

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
female$followup <- "NA"

# combine analyses ====
overall <- rbind(sex_combined,male,female)

# site specific ====
# sex combined  ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_combined

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
proximal$followup <- "NA"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_combined

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
distal$followup <- "NA"

## combine analyses ====
results_combined <- rbind(proximal, distal)
results_combined <- results_combined[order(results_combined$exposure),]

# male ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_male

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
proximal$followup <- "NA"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_male

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
distal$followup <- "NA"

## combine analyses ====
results_male <- rbind(proximal, distal)
results_male <- results_male[order(results_male$exposure),]

# female ====
# proximal ====
results_list <- list()
data_analysis <- data_analysis_proximal_female

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
proximal$followup <- "NA"

# distal ====
results_list <- list()
data_analysis <- data_analysis_distal_female

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
distal$followup <- "NA"

## combine analyses ====
results_female <- rbind(proximal, distal)
results_female <- results_female[order(results_female$exposure),]

# combine and save ====
results <- rbind(overall, results_combined, results_male, results_female)
results <- results[order(results$exposure),]
write.table(results, "analysis/002_analysis/metabolomics/positive/001_analysis_sex_site-specific.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
