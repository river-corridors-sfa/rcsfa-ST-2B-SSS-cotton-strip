# ==============================================================================
#
# Calculated weighted average for each site
# 
# Need to run "SSS_Sonar_Depth_Processing.Rmd" for each site before running this 
# code because the outputs from that script are not included in this data package
#
# ==============================================================================
#
# Author: Brieanne Forbes (brieanne.forbes@pnnl.gov)
# 20 September 2024

# remove all files
rm(list=ls(all=TRUE))

# ==============================================================================

library(tidyverse)

#set working directory to this folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

# ============================ find and read files =============================

files <- list.files('.', '_mean', full.names = T)

IDs <- c('SSS014', 'SSS025', 'SSS026', 'SSS047')

combine <- tibble(Parent_ID = as.character(),
                  weighted_mean_cm = as.numeric(),
                  total_bins = as.numeric())

for (ID in IDs) {
  
  ID_files <- files[grepl(ID, files)] # select files for one site
  
  means <- ID_files %>%
    map(~ read_csv(.x) %>%
          rename(Mean = 1,
                 Bins_count = 2)) %>%
    bind_rows()
  
  #calculate a weighted mean. Weighted by the amount of spatial (1m x 1m) bins 
  weighted_mean_cm <- weighted.mean(means$Mean, means$Bins_count) * 100 # converting from m to cm
  
  combine <- combine %>%
    add_row(Parent_ID = ID, 
            weighted_mean_cm = weighted_mean_cm,
            total_bins = sum(means$Bins_count))
  
}

# the values in the "combine" dataframe are reported in "v*_SSS_Water_Depth_Summary.csv"
