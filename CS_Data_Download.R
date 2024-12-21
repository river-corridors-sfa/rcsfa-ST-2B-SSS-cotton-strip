# Load devtools package
library(devtools)

# Specify the GitHub URL of the script
github_url <- "https://raw.githubusercontent.com/river-corridors-sfa/rcsfa-data_processing_for_publication/refs/heads/main/Data_Package_ESS-DIVE/download_from_ESS-DIVE_landing_page/ESS-DIVE_Download_R/script_ess_dive_file_download_function.R"

# Source the script from GitHub
source_url(github_url)

csv_files_from_data_package <- download_and_read_data(
  target_url = "https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-f8b4bf15fc21108-20241107T161454014", # the url from ess-dive
  filename = "v3_SSS_Data_Package.zip", # your choice for how the downloaded zip will be named
  downloads_folder = getwd())
