library(dplyr)
library(tidyproteomics)

load_all()
str_break <- "\n------\n"

test_path <- "~/Local/gitdata/tidyproteomics_data/"
test_files <- c(
  # "FragPipe_19.1/combined_peptide.tsv",
  "FragPipe_19.1/combined_protein.tsv",
  # "ProteomeDiscoverer_2.5/p97KD_HCT116_peptides.xlsx",
  # "ProteomeDiscoverer_2.5/p97KD_HCT116_proteins.xlsx",
  # "MaxQuant_1.6.10.43/txt/evidence.txt",
  # "MaxQuant_1.6.10.43/txt/proteinGroups.txt",
  # "Skyline/dda_tryp_msstats_peptides.csv",
  # "DIA_NN_1.8.1/tryp-report.tsv",
  "PXD004163/miR_Proteintable.tsv"
)

platforms <- c("FragPipe", "ProteomeDiscoverer", "MaxQuant", "Skyline", "DIA-NN", 'mzTab')

data_list <- list()
for( i in 1:length(test_files) ){

  file_name <- test_files[i]
  platform <- NULL
  path <- NULL

  for(platform in platforms){ if(grepl(platform, sub("_", "-", file_name))){break} }
  if(grepl('PXD004163', file_name)){
    platform <- 'UserDef_PXD004163'
    path <- sub("miR_Proteintable", "TMTOpenMS_proteins", paste0(test_path, file_name))
  }

  analyte = 'peptides'
  if(grepl('protein', file_name)) {analyte = 'proteins'}

  cat(str_break)
  file_names <- paste0(test_path, file_name)
  data_list[[i]] <- import(file_names, platform, analyte, path)

  if(data_list[[i]]$quantitative %>% filter(!is.na(abundance_raw)) %>% nrow() <= 1){
    cli::cli_abort("Issue with {platform}:{analyte} > {file_name}")
  }
}
