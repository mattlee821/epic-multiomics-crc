library(metaboprep)

# proteomics ====
output_dir_path = paste0("/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/data/proteomics/metaboprep_release_2023_03_30/")
rdfile = paste0(output_dir_path, "ReportData.Rdata")
setwd("/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/scripts/001_data_formatting/002_metaboprep/proteomics/")
generate_report(full_path_2_Rdatafile = rdfile, dir_4_report = output_dir_path,
                path_2_Rmd_template = "metaboprep_Report_v0.Rmd")

# metabolites - positive ====
output_dir_path = paste0("/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/data/metabolomics/positive/metaboprep_release_2023_03_30/")
rdfile = paste0(output_dir_path, "ReportData.Rdata")
setwd("/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/scripts/001_data_formatting/002_metaboprep/metabolomics/positive/")
generate_report(full_path_2_Rdatafile = rdfile, dir_4_report = output_dir_path,
                path_2_Rmd_template = "metaboprep_Report_v0.Rmd")

# metabolites - negative ====
output_dir_path = paste0("/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/data/metabolomics/negative/metaboprep_release_2023_03_30/")
rdfile = paste0(output_dir_path, "ReportData.Rdata")
setwd("/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/scripts/001_data_formatting/002_metaboprep/metabolomics/negative/")
generate_report(full_path_2_Rdatafile = rdfile, dir_4_report = output_dir_path,
                path_2_Rmd_template = "metaboprep_Report_v0.Rmd")
