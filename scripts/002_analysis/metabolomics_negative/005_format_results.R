rm(list=ls())
set.seed(821)

# data ====
a <- read.table("analysis/002_analysis/metabolomics/negative/001_analysis_sex_site-specific.txt", header = T, sep = "\t")
b <- read.table("analysis/002_analysis/metabolomics/negative/002_analysis_sex_site-specific_median-followup.txt", header = T, sep = "\t")
data <- bind_rows(a,b)

# format ====
## OR and CI
data$OR <- exp(data$b)
data$lower_ci <- exp(data$b - (1.96 * data$se))
data$upper_ci <- exp(data$b + (1.96 * data$se))
data$analysis <- "main"

# save
write.table(data, "analysis/002_analysis/metabolomics/negative/results_formatted.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# sensitivity ====
a <- read.table("analysis/002_analysis/metabolomics/negative/001_analysis_sex_site-specific_removed-2y-followup.txt", header = T, sep = "\t")
b <- read.table("analysis/002_analysis/metabolomics/negative/002_analysis_sex_site-specific_median-followup_removed-2y-followup.txt", header = T, sep = "\t")
data <- bind_rows(a,b)

# format ====
## OR and CI
data$OR <- exp(data$b)
data$lower_ci <- exp(data$b - (1.96 * data$se))
data$upper_ci <- exp(data$b + (1.96 * data$se))

# save
write.table(data, "analysis/002_analysis/metabolomics/negative/results_2y-followup_formatted.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
