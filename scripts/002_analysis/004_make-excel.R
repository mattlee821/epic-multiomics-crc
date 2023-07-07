rm(list=ls())

files <- list.files(path = "analysis/002_analysis/", pattern = "results", recursive = T, full.names = T, include.dirs = T)
data <- lapply(files, read.table, header = T, sep = "\t")


# Install and load the 'openxlsx' package
install.packages("openxlsx")
library(openxlsx)

sheet_names <- c("met_neg_2yfollowup", "met_neg",
                 "met_pos_2yfollowup", "met_pos",
                 "prot_2yfollowup", "prot")

# Specify the file path and name for the Excel file
filepath <- "analysis/002_analysis/epic_crc_proteomics_metabolomics_results.xlsx"

# Create a new workbook
wb <- createWorkbook()

for (i in seq_along(data)) {
  addWorksheet(wb, sheetName = sheet_names[i])
  writeData(wb, sheet = i, x = data[[i]])
}

# Save the workbook as an Excel file
saveWorkbook(wb, filepath)

