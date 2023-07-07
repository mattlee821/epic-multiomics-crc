PROTEINS=/home/leem//001_projects/epic_metabolomics_proteomics/scripts/001_data_formatting/002_metaboprep/proteomics/
METABOLITES_POS=/home/leem//001_projects/epic_metabolomics_proteomics/scripts/001_data_formatting/002_metaboprep/metabolomics/positive/
METABOLITES_NEG=/home/leem//001_projects/epic_metabolomics_proteomics/scripts/001_data_formatting/002_metaboprep/metabolomics/negative/

# if running on HPC likely to get an error at the end when generating a report. to make report, copy newly made folder locally and run:
# output_dir_path = paste0("FULL/PATH/TO/DIR/OF/CHOICE/")
# rdfile = paste0(output_dir_path, "ReportData.Rdata")
# generate_report( full_path_2_Rdatafile = rdfile, dir_4_report = output_dir_path )
# when finished, copy over the new html report from metaboprep directory AND in the script location the new .md file and the figures/

cd ${PROTEINS}
Rscript ${PROTEINS}run_metaboprep_pipeline.R ${PROTEINS}parameter_file.txt

cd ${METABOLITES_POS}
Rscript ${METABOLITES_POS}run_metaboprep_pipeline.R ${METABOLITES_POS}parameter_file.txt

cd ${METABOLITES_NEG}
Rscript ${METABOLITES_NEG}run_metaboprep_pipeline.R ${METABOLITES_NEG}parameter_file.txt
