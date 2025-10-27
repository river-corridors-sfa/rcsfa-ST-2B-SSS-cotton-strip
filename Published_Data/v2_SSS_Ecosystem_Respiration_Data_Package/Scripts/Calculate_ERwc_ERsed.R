# ==============================================================================
#
# Calculate SSS ERwc and ERsed
#
# ==============================================================================
#
# Author: Brieanne Forbes 
# 23 January 2025
#
# ==============================================================================

library(tidyverse)

rm(list=ls(all=T))

# ================================= Set WD and read files ================================

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./../")

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
ERwc_slope <- read_csv('./Published_Data/v3_SSS_Data_Package/Sensor_Data/miniDOTManualChamber/Plots_and_Summary_Statistics/v3_SSS_MC_Water_DO_Temp_Summary.csv', comment='#', na = c('N/A', -9999)) %>%
  select(Parent_ID,Site_ID,Dissolved_Oxygen_1_Slope,
         Dissolved_Oxygen_2_Slope,Dissolved_Oxygen_3_Slope,
         Dissolved_Oxygen_1_NRMSE,Dissolved_Oxygen_2_NRMSE,Dissolved_Oxygen_3_NRMSE)

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
depth <- read_csv('./Published_Data/v3_SSS_Data_Package/v3_SSS_Water_Depth_Summary.csv', comment='#', na = c('N/A', -9999)) %>%
  select(Parent_ID, Site_ID, Average_Depth)


ERtotal <- read_csv('./Stream_Metabolizer/Outputs/v2_SSS_combined_SM_results.csv', na = '-9999') %>%
  select(Parent_ID, Site_ID, `ERdailymeanmean_gO2/m2day`, `GPPdailymeanmean_gO2/m2day`, daysofdata)

# =============================== calculate ERwc and ERsed ===================================

# ERwc conversion
# mg/L = g/m3; multiply by m (depth) to get g/m2
# min -> day = multiply by (60min * 24hr)

all_data <- ERwc_slope %>%
  full_join(depth) %>%
  full_join(ERtotal) %>%
  mutate(Average_Depth_m = Average_Depth/100, # convert depth from cm to m
         Dissolved_Oxygen_Slope_Mean_mg_per_L_per_min = rowMeans(select(., contains("Slope")), na.rm = TRUE), # calculate mean slope
         ERwc_g_per_m2_per_day = Dissolved_Oxygen_Slope_Mean_mg_per_L_per_min*(60*24)*Average_Depth_m) %>% # convert ERwc
  rename(Total_Ecosystem_Respiration = "ERdailymeanmean_gO2/m2day",
         Water_Column_Respiration = ERwc_g_per_m2_per_day,
         Gross_Primary_Production = "GPPdailymeanmean_gO2/m2day",
         ERtot_Days_of_Data = daysofdata) %>%
  mutate(Sediment_Respiration = case_when(Water_Column_Respiration > 0 ~ round(Total_Ecosystem_Respiration - 0, 2),
                                          TRUE ~ round(Total_Ecosystem_Respiration - Water_Column_Respiration, 2)), # calculate ERsed; If ERwc is positive, set ERwc to 0 before calculating
         Sediment_Respiration = case_when(Sediment_Respiration > 0 ~ 0,
                                          TRUE ~ Sediment_Respiration),
         Total_Ecosystem_Respiration = round(Total_Ecosystem_Respiration, 2),
         Water_Column_Respiration = signif(Water_Column_Respiration, 3),
         Gross_Primary_Production = round(Gross_Primary_Production, 2),
         ERtot_Days_of_Data = case_when(is.na(Total_Ecosystem_Respiration) ~ NA,
                                                  TRUE ~ ERtot_Days_of_Data), # remove days of data where we don't report values because model failed
         Sediment_Respiration = case_when(Total_Ecosystem_Respiration > 0 ~ NA, # remove ERsed when ERtot is positive 
                                                  TRUE ~ Sediment_Respiration),
         Sediment_Respiration_Contribution = round((Sediment_Respiration/Total_Ecosystem_Respiration), 2),
         Water_Column_Respiration_Contribution = case_when(is.na(Sediment_Respiration_Contribution) ~ NA, # when ERsed contribution in NA, make ERwc contribution NA
                                                           Sediment_Respiration_Contribution == 0 ~ 1, # when ERsed contribution is 0, make ERwc contribution 1. This happens when ERwc > ERtot
                                                           Total_Ecosystem_Respiration < 0 & Water_Column_Respiration < 0 ~ (Water_Column_Respiration/Total_Ecosystem_Respiration), # when ERtot and ERwc are negative, calculate as normal
                                                           Water_Column_Respiration > 0 ~ 0, # when ERwc is positive, make ERwc contribution 0
                                                           TRUE ~ NA),
         Water_Column_Respiration_Contribution = round(Water_Column_Respiration_Contribution, 2))

final_data <- all_data %>%
  select(Parent_ID, Site_ID, Total_Ecosystem_Respiration, Water_Column_Respiration, Sediment_Respiration, Gross_Primary_Production, ERtot_Days_of_Data, Sediment_Respiration_Contribution, Water_Column_Respiration_Contribution)

header <- tibble(x1 = c('# HeaderRows_9',
                        '# HeaderRows_Format: Column_Header; Unit; InstallationMethod_ID; Instrument_Summary',
                        '# Total_Ecosystem_Respiration; grams_O2_per_square_meter_per_day; N/A; Estimated via streamMetabolizer model.',
                        '# Water_Column_Respiration; grams_O2_per_square_meter_per_day; ManualChamberDark_02; PME MiniDOT Logger with optical dissolved oxygen sensor (fluorescence quenching).',
                        '# Sediment_Respiration; grams_O2_per_square_meter_per_day; N/A; Calculated using Ecosystem Respiration from the streamMetabolizer model and Water Column Respiration from PME MiniDOT logger manual chamber.',
                        '# Gross_Primary_Production; grams_O2_per_square_meter_per_day; N/A; Estimated via streamMetabolizer model.',
                        '# Sediment_Respiration_Contribution; N/A; N/A; Calculated with Total_Ecosystem_Respiration and Sediment_Respiration',
                        '# Water_Column_Respiration_Contribution; N/A; N/A; Calculated with Total_Ecosystem_Respiration and Water_Column_Respiration.')
)

write_csv(header, './v2_SSS_Water_Sediment_Total_Respiration_GPP.csv', append = F, col_names = F)

write_csv(final_data,'./v2_SSS_Water_Sediment_Total_Respiration_GPP.csv', na = '-9999', append = T, col_names = T)
