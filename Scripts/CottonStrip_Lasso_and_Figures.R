# ==============================================================================
#
# LASSO analysis for SSS Cotton Strips
#
# Status: Complete
#
# ==============================================================================
#
# Author: Brieanne Forbes 
# 27 Oct 2025
#
# ==============================================================================
library(tidyverse) 
library(corrplot)
library(glmnet)
library(ggpmisc)

rm(list=ls(all=T))

# Setting wd to parent folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd('./../')

# =================================== user input ===============================

response_variable <- c('scale_cube_Mean_Decay_Rate_per_day', 'scale_cube_Mean_degree_decay_rate')

# =================================== find files ===============================
 # degree_decay_rate = Kdd; Decay_Rate_per_day = Kcc
decay_temp <- read_csv('./Outputs/Decay_Data.csv') %>%
  select(Parent_ID, sum_mean_daily_temp, degree_decay_rate, Decay_Rate_per_day )

decay_temp_means <- decay_temp%>%
  group_by(Parent_ID) %>%
  summarise(Mean_sum_mean_daily_temp = mean(sum_mean_daily_temp),
            Mean_degree_decay_rate = mean(degree_decay_rate),
            Mean_Decay_Rate_per_day = mean(Decay_Rate_per_day),
  ) %>%
  ungroup()



er_gpp <- read_csv('./Published_Data/v2_SSS_Ecosystem_Respiration_Data_Package/v2_SSS_Water_Sediment_Total_Respiration_GPP.csv',
                   comment = '#', na = '-9999') %>%
  select(Parent_ID, Site_ID, Sediment_Respiration, Total_Ecosystem_Respiration, Water_Column_Respiration, Gross_Primary_Production)

d50 <- read_csv('./Published_Data/v2_SSS_Ecosystem_Respiration_Data_Package/v2_SSS_ER_d50_TotalOxygenConsumed.csv', 
                comment = '#', na = '-9999') %>%
  select(Parent_ID, D50_m)

slope_vel_dis <- read_csv('./Published_Data/v2_SSS_Ecosystem_Respiration_Data_Package/Stream_Metabolizer/Inputs/v2_SSS_Slope_Discharge_Velocity.csv',
                          comment = '#', na = '-9999') %>%
  select(Site_ID, Slope, Discharge, Velocity)

# NHD+, streamcat, NLCD, ET0 extracted geospatial variables https://github.com/river-corridors-sfa/Geospatial_variables
geospatial <- read_csv('https://github.com/river-corridors-sfa/Geospatial_variables/raw/refs/heads/main/v4_RCSFA_Extracted_Geospatial_Data_2025-01-31.csv') %>%
  select(site, totdasqkm, pctmxfst2019ws,pctconif2019ws,pctdecid2019ws, AridityWs, pctcrop2019ws, pcthay2019ws, pctshrb2019ws, minelevsmo) %>%
  filter(site %in% er_gpp$Site_ID)%>%
  mutate(PctFst = pctmxfst2019ws + pctdecid2019ws + pctconif2019ws,
         PctAg = pctcrop2019ws + pcthay2019ws) %>%
  rename(Site_ID = site)%>%
  select(Site_ID, totdasqkm, PctFst, AridityWs, PctAg, pctshrb2019ws, minelevsmo)

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
tss <- read_csv('./Published_Data/v3_SSS_Data_Package/Sample_Data/SSS_Water_TSS.csv',
                skip = 2, na = c('', 'N/A', '-9999')) %>%
  filter(!is.na(Sample_Name)) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         '00530_TSS_mg_per_L' = case_when(str_detect(`00530_TSS_mg_per_L`, 'LOD') ~ 0.12, # replace below LOD values with half LOD (LOD = 0.24)
                                          TRUE ~ as.numeric(`00530_TSS_mg_per_L`)))%>%
  select(Parent_ID, contains('TSS')) 

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
depth <- read_csv('./Published_Data/v3_SSS_Data_Package/v3_SSS_Water_Depth_Summary.csv',
                  comment = '#', na = c('', 'N/A', '-9999')) %>%
  select(Parent_ID, Average_Depth)

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
water_npoc_tn <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Water_NPOC_TN.csv',
                    skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         '00602_TN_mg_per_L_as_N' = case_when(str_detect(`00602_TN_mg_per_L_as_N`, 'LOD') ~ 0.013, # replace below LOD values with half LOD (LOD = 0.026)
                                              str_detect(`00602_TN_mg_per_L_as_N`, 'Standard') ~ 0.05, # replace below standard values with half standard (standard = 0.1)
                                              TRUE ~ as.numeric(`00602_TN_mg_per_L_as_N`)),
         `00681_NPOC_mg_per_L_as_C` = as.numeric(`00681_NPOC_mg_per_L_as_C`)) %>%
  select(Parent_ID, contains('NPOC'), contains('TN')) %>%
  group_by(Parent_ID) %>%
  summarise(Mean_00602_TN_mg_per_L_as_N = round(mean(`00602_TN_mg_per_L_as_N`), 2),
            Mean_00681_NPOC_mg_per_L_as_C = round(mean(`00681_NPOC_mg_per_L_as_C`), 2),
  ) %>%
  ungroup()

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
sed_npoc_tn <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/CM_SSS_Sediment_NPOC_TN.csv',
                          skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         'Extractable_TN_mg_per_kg' = case_when(str_detect(`Extractable_TN_mg_per_kg`, 'Standard') ~ 0.05, # replace below standard values with half standard (standard = 0.1)
                                              TRUE ~ as.numeric(`Extractable_TN_mg_per_kg`)),
         Extractable_NPOC_mg_per_kg = as.numeric(Extractable_NPOC_mg_per_kg)) %>%
  select(Parent_ID, Extractable_NPOC_mg_per_kg, Extractable_TN_mg_per_kg) %>%
  group_by(Parent_ID) %>%
  summarise(Mean_Extractable_TN_mg_per_kg = round(mean(`Extractable_TN_mg_per_kg`), 2),
            Mean_Extractable_NPOC_mg_per_kg = round(mean(`Extractable_NPOC_mg_per_kg`), 2)) %>%
  ungroup()


#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
cn <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/CM_SSS_Sediment_CN.csv',
                          skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         `01395_C_percent_per_mg` = as.numeric(`01395_C_percent_per_mg`),
         `01397_N_percent_per_mg` = as.numeric(`01397_N_percent_per_mg`)) %>%
  select(Parent_ID, `01395_C_percent_per_mg`, `01397_N_percent_per_mg`)

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
ions <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/CM_SSS_Water_Ions.csv',
               skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}")) %>%
  select(Parent_ID, `00915_Ca_mg_per_L`, `00940_Cl_mg_per_L`, `00945_SO4_mg_per_L_as_SO4`)%>%
  group_by(Parent_ID) %>%
  summarise(Mean_00915_Ca_mg_per_L  = round(mean(as.numeric(`00915_Ca_mg_per_L` ), na.rm = T), 2),
            Mean_00940_Cl_mg_per_L  = round(mean(as.numeric(`00940_Cl_mg_per_L` ), na.rm = T), 2),
            Mean_00945_SO4_mg_per_L_as_SO4  = round(mean(as.numeric(`00945_SO4_mg_per_L_as_SO4` ), na.rm = T), 2)) %>%
  ungroup()


#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
resp <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/v2_CM_SSS_Sediment_Normalized_Respiration_Rates.csv',
                 skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}")) %>%
  group_by(Parent_ID) %>%
  summarise(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment  = round(mean(as.numeric(Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment ), na.rm = T), 2)) %>%
  ungroup()


#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
iron <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Fe.csv',
                 skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}")) %>%
  group_by(Parent_ID) %>%
  summarise(Mean_Fe_mg_per_kg = round(mean(as.numeric(Fe_mg_per_kg), na.rm = T), 2)) %>%
  ungroup()

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
sand <- read_csv('./Published_Data/v5_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Grain_Size.csv',
               skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         Percent_Tot_Sand = as.numeric(Percent_Tot_Sand)) %>%
  select(Parent_ID, Percent_Tot_Sand)

# =============================== combine data ===============================

all_data <- decay_temp_means %>%
  full_join(er_gpp, by = 'Parent_ID')%>%
  full_join(d50, by = 'Parent_ID')%>%
  full_join(slope_vel_dis, by = 'Site_ID')%>%
  full_join(geospatial, by = 'Site_ID')%>%
  full_join(tss, by = 'Parent_ID')%>%
  full_join(depth, by = 'Parent_ID')%>%
  full_join(water_npoc_tn, by = 'Parent_ID')%>%
  full_join(sed_npoc_tn, by = 'Parent_ID')%>%
  full_join(cn, by = 'Parent_ID')%>%
  full_join(ions, by = 'Parent_ID')%>%
  full_join(resp, by = 'Parent_ID')%>%
  full_join(iron, by = 'Parent_ID')%>%
  full_join(sand, by = 'Parent_ID') %>%
  filter(!is.na(Sediment_Respiration)) %>%
  filter(!is.na(Mean_degree_decay_rate)) # dropping sites that are missing ERsed or decay data

#check that there are no NA values
# if returns any obs, there are NA values 
# if returns 0 obs, no NA values 
# returning 0 x 33 so no NA values 
all_data %>%
  filter(if_any(everything(), is.na))

# ======================= assess co-correlation ===============================

long_data <-  all_data %>% 
  pivot_longer(cols = -c(Site_ID, Parent_ID), names_to = "variable", values_to = "value")

# ggplot() + 
#   geom_histogram(long_data, mapping = aes(x = value)) + 
#   facet_wrap(~ variable, scales = "free") +

## ======== Variable Name Mapping ===========
# Define a tibble mapping original, scaled, and plot labels
variable_names <- tibble(
  original = c("Mean_sum_mean_daily_temp", "Mean_degree_decay_rate", "Mean_Decay_Rate_per_day", "Sediment_Respiration",
               "Total_Ecosystem_Respiration", "Water_Column_Respiration", "Gross_Primary_Production", "D50_m", "Slope",
               "Discharge", "Velocity", "totdasqkm", "PctFst", "AridityWs", "PctAg", "pctshrb2019ws", "minelevsmo",
               "00530_TSS_mg_per_L", "Average_Depth", "Mean_00602_TN_mg_per_L_as_N", "Mean_00681_NPOC_mg_per_L_as_C",
               "Mean_Extractable_NPOC_mg_per_kg", "Mean_Extractable_TN_mg_per_kg", "01395_C_percent_per_mg", "01397_N_percent_per_mg",
               "Mean_00915_Ca_mg_per_L", "Mean_00940_Cl_mg_per_L", "Mean_00945_SO4_mg_per_L_as_SO4", "Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment",
               "Mean_Fe_mg_per_kg", "Percent_Tot_Sand"),
  labels = c("Summed Temperature", "Decay Rate per degree day (Kdd)", "Decay Rate per chronological day (Kcd)", "Sediment Respiration",
             "Total Ecosystem Respiration", "Water Column Respiration", "Gross Primary Production", "D50", "Slope", "Discharge",
             "Velocity", "Drainage Area", "Percent Forest Cover", "Aridity", "Percent Agricultural Cover", "Percent Shrub Cover",
             "Minimum Elevation", "TSS", "Depth", "Water TDN", "Water NPOC", "Sediment NPOC", "Sediment TN", "Percent C", "Percent N",
             "Calcium", "Chloride", "Sulfate", "Respiration (Lab)", "Iron", "Percent Sand") 
) %>%
  mutate(cubed = paste0("cube_",original))

## ======== Cube Root Transformation ===========

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

cube_data <- all_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(.fn = ~ paste0("cube_", .x), .cols = where(is.numeric))

## ======== Histogram of Cubed Data ===========
long_cube_data <- cube_data %>%
  pivot_longer(cols = -c(Site_ID, Parent_ID), names_to = "variable", values_to = "value")

# ggplot() +
#   geom_histogram(long_cube_data, mapping = aes(x = value)) +
#   facet_wrap(~ variable, scales = "free") +
#   theme_minimal()


## ======== Pearson Correlation Matrix with Cube Transformation ===========
renamed_cube_data <- cube_data %>%
  select(-Site_ID, -Parent_ID) %>%
  rename_with(~ ifelse(!is.na(match(., variable_names$cubed)),
                       variable_names$labels[match(., variable_names$cubed)],
                       .), .cols = names(cube_data %>% select(-Site_ID, -Parent_ID)))

pearson_cubed <- cor(renamed_cube_data, method = "pearson", use = "complete.obs")

rdylbu_colors <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090",
                   "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")

png(file = "./Figures/FigureS1_Pearson_Correlation_Matrix_Cubed.png", width = 12, height = 12, units = "in", res = 300)
corrplot(pearson_cubed, type = "upper", method = "number", 
         col = colorRampPalette(rdylbu_colors)(200), tl.col = "black", tl.cex = 0.5, 
         number.cex = 0.5, cl.cex = 1.25, mar = c(0, 0, 2, 0), bg = "black")
dev.off()

# ======== LASSO  ============

## scale data
scale_cube_variables = as.data.frame(scale(cube_data %>% select(-Parent_ID, -Site_ID)))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x))

## Loop through LASSO to get average over a lot of seeds ####
for (variable in response_variable) {
  

num_seeds = 100
seeds = sample(1:500, num_seeds)

## Set response variable (scale_cube_Mean_Decay_Rate_per_day/scale_cube_Mean_degree_decay_rate) and scale
yvar <- data.matrix(scale_cube_variables %>% pull(variable))
round(mean(yvar), 4)
sd(yvar)

# list for storing LASSO iterations
norm_coeffs = list()
lasso_coefs_pull = list()
r2_scores = numeric(num_seeds)

## Set predictor variables; exclude response variable(s)

x_cube_variables = scale_cube_variables %>%
  select(-response_variable)

if(variable == 'scale_cube_Mean_degree_decay_rate'){
  
  x_cube_variables = x_cube_variables %>%
    select( -scale_cube_Mean_sum_mean_daily_temp)
}

xvars <- data.matrix(x_cube_variables)


for (i in 1:num_seeds) {
  
  seed = seeds[i]
  set.seed(seed)
  
  lasso = cv.glmnet(xvars, yvar, alpha = 1, nfolds = 5,
                    standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                    #,standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                    # , standardize = TRUE, standardize.response = FALSE, intercept = FALSE
  )
  
  best_lambda <- lasso$lambda.min
  #best_lambda
  #plot(lasso)
  
  best_lasso_model <- glmnet(xvars, yvar, alpha = 1, lambda = best_lambda, family = "gaussian",
                             standardize = FALSE, standardize.response = FALSE, intercept = FALSE
                             #  , standardize = TRUE, standardize.response = TRUE, intercept = FALSE
                             #, standardize = TRUE, standardize.response = FALSE, intercept = FALSE
  )
  
  
  lasso_coefs = as.matrix(coef(best_lasso_model, s = best_lambda))
  
  lasso_coefs_pull[[as.character(seed)]] = lasso_coefs[-1, , drop = FALSE]
  
  norm_coeffs_scale = lasso_coefs/max(abs(lasso_coefs[-1]))
  
  norm_coeffs[[as.character(seed)]] = norm_coeffs_scale[-1, , drop = FALSE]
  
  y_pred = predict(best_lasso_model, newx = xvars, s = best_lambda)
  
  sst = sum((yvar - mean(yvar))^2)
  sse = sum((y_pred - yvar)^2)
  r2_scores[i] = 1 - (sse / sst)
  
}

lasso_coef_mat = as.data.frame(do.call(cbind, lasso_coefs_pull)) 
colnames(lasso_coef_mat) <- paste0("s", seq_len(ncol(lasso_coef_mat)))
# Make DF of all LASSO results with mean and std. dev  
lasso_coef_means = lasso_coef_mat %>% 
  mutate(RowNames = rownames(lasso_coef_mat)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1"))),
         cv = sd/mean) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)%>% 
  relocate(cv, .after = sd) %>%
  add_column(response_variable = variable)

norm_coeffs_matrix = do.call(cbind, norm_coeffs)

mean_coeffs = as.data.frame(norm_coeffs_matrix, row.names = rownames(norm_coeffs_matrix))
colnames(mean_coeffs) <- paste0("s", seq_len(ncol(mean_coeffs)))

norm_lasso_coef_means = mean_coeffs %>% 
  mutate(RowNames = rownames(mean_coeffs)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1"))),
         cv = sd/mean) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)%>% 
  relocate(cv, .after = sd) %>% 
  add_column(response_variable = variable)

results_r2 = as.data.frame(r2_scores) 
mean(results_r2$r2_scores)
sd(results_r2$r2_scores)

if(match(variable, response_variable) == 1){
  
  lasso_coef_means_all <- lasso_coef_means
  norm_lasso_coef_means_all <- norm_lasso_coef_means
  mean_r2_all <- tibble(mean_r2 = mean(results_r2$r2_scores),
                        sd = sd(results_r2$r2_scores),
                        response_variable = variable)
  
} else{
  
  lasso_coef_means_all <- lasso_coef_means_all %>%
    add_row(lasso_coef_means)
  
  norm_lasso_coef_means_all <- norm_lasso_coef_means_all %>%
    add_row(norm_lasso_coef_means)
  
  mean_r2_all <- mean_r2_all %>%
    add_row(mean_r2 = mean(results_r2$r2_scores),
           sd = sd(results_r2$r2_scores),
           response_variable = variable)
}



}

# ================================ investigate cv ==============================

all_results_long <- bind_rows(
  lasso_coef_means_all %>%
    select(RowNames, mean, sd, cv, response_variable) %>%
    add_column(type = 'Not_Normalized'),
  
  norm_lasso_coef_means_all %>%
    select(RowNames, mean, sd, cv, response_variable) %>%
    add_column(type = 'Normalized')
) %>%
  mutate(cv = round(cv, 3))

# absolute cv vs absolute mean; all norm/not norm + response variables
ggplot(data = all_results_long, aes(x = abs(cv), y = abs(mean))) + 
  geom_point()+
  theme_bw()

# absolute cv vs absolute mean; all norm/not norm + response variables; filtered cv <= 1
ggplot(data = all_results_long %>% filter(abs(cv) <= 1), aes(x = abs(cv), y = abs(mean))) + 
  geom_point()+
  theme_bw()


# absolute cv vs absolute mean; pivoted by norm/not norm and response variable 
cv_plot <- ggplot(data = all_results_long %>% mutate(response_variable = case_when(response_variable == 'scale_cube_Mean_degree_decay_rate' ~ 'Kdd',
                                                                        response_variable == 'scale_cube_Mean_Decay_Rate_per_day' ~ 'Kcd'),
                                          type = str_replace(type, 'Not_Normalized', 'Not Normalized' )), aes(x = abs(cv), y = abs(mean))) + 
  geom_point()+
  facet_grid(response_variable ~ type)+
  theme_bw()


# ggsave(
#   paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()), "_Mean_vs_CV.png"),
#   cv_plot,
#   width = 8,
#   height = 8,
#   units = 'in',
#   dpi = 300
# )

# absolute cv vs absolute mean; pivoted by norm/not norm and response variable; filtered cv <= 1
ggplot(data = all_results_long %>% filter(abs(cv) <= 1), aes(x = abs(cv), y = abs(mean))) + 
  geom_point()+
  facet_grid(response_variable ~ type)+
  theme_bw()

#cv histo; all norm/not norm + response variables
ggplot(data = all_results_long, aes(x = abs(cv))) + 
  geom_histogram()+
  theme_bw()

#cv histo; all norm/not norm + response variables; filtered cv <= 1
ggplot(data = all_results_long %>% filter(abs(cv) <= 1), aes(x = abs(cv))) + 
  geom_histogram()+
  theme_bw()

#cv histo; pivoted by norm/not norm and response variable
ggplot(data = all_results_long, aes(x = abs(cv))) + 
  geom_histogram()+
  facet_grid(response_variable ~ type)+
  theme_bw()

#cv histo; pivoted by norm/not norm and response variable; filtered cv <= 1
ggplot(data = all_results_long %>% filter(abs(cv) <= 1), aes(x = abs(cv))) + 
  geom_histogram()+
  facet_grid(response_variable ~ type)+
  theme_bw()

#cv rank; norm decay rate per day
filtered_data1 <- all_results_long %>%
  filter(response_variable == 'scale_cube_Mean_Decay_Rate_per_day',
         type == 'Normalized') %>%
  arrange(abs(cv)) %>%  # Ensure order before setting factor
  mutate(RowNames = factor(RowNames, levels = unique(RowNames)))  # Explicit ordering
ggplot(filtered_data1, aes(x = abs(cv), y = RowNames)) + 
  geom_point() +
  theme_bw()+
  ggtitle("norm decay rate per day")

#cv rank; not norm decay rate per day
filtered_data2 <- all_results_long %>%
  filter(response_variable == 'scale_cube_Mean_Decay_Rate_per_day',
         type == 'Not_Normalized') %>%
  arrange(abs(cv)) %>%  # Ensure order before setting factor
  mutate(RowNames = factor(RowNames, levels = unique(RowNames)))  # Explicit ordering

ggplot(filtered_data2, aes(x = abs(cv), y = RowNames)) + 
  geom_point() +
  theme_bw()+
  ggtitle("not norm decay rate per day")

#cv rank; norm degree decay rate
filtered_data3 <- all_results_long %>%
  filter(response_variable == 'scale_cube_Mean_degree_decay_rate',
         type == 'Normalized') %>%
  arrange(abs(cv)) %>%  # Ensure order before setting factor
  mutate(RowNames = factor(RowNames, levels = unique(RowNames)))  # Explicit ordering

ggplot(filtered_data3, aes(x = abs(cv), y = RowNames)) + 
  geom_point() +
  theme_bw()+
  ggtitle("norm degree decay rate")

#cv rank; not norm degree decay rate
filtered_data4 <- all_results_long %>%
  filter(response_variable == 'scale_cube_Mean_degree_decay_rate',
         type == 'Not_Normalized') %>%
  arrange(abs(cv)) %>%  # Ensure order before setting factor
  mutate(RowNames = factor(RowNames, levels = unique(RowNames)))  # Explicit ordering

ggplot(filtered_data4, aes(x = abs(cv), y = RowNames)) + 
  geom_point() +
  theme_bw()+
  ggtitle("not norm degree decay rate")

# ================================ create out table ============================

output <- all_results_long %>%
  mutate(mean = signif(mean, 3),
         abs_mean = abs(mean),
         sd = signif(sd, 3),
         cv = signif(cv, 3),
         abs_cv = abs(cv),
         cv = case_when(is.na(cv) ~ '',
                        TRUE ~ as.character(cv)),
         RowNames = str_remove(RowNames, 'scale_cube_'),
         response_variable = case_when(response_variable == 'scale_cube_Mean_degree_decay_rate' ~ 'Kdd',
                                       response_variable == 'scale_cube_Mean_Decay_Rate_per_day' ~ 'Kcd'),
         type = str_replace(type, 'Not_Normalized', 'Not Normalized' )) %>%
  left_join(variable_names %>% select(original, labels), by = c('RowNames' = 'original')) %>%
  rename(Predictor = labels)%>% 
  # group_by(response_variable, type) %>%
  # arrange(desc(abs_mean)) %>%
  # ungroup() %>%
  # select(response_variable, type, Predictor, mean, sd, cv)
  clipr::write_clip()
  
out_r2 <- mean_r2_all %>%
  mutate(mean_r2 = signif(mean_r2, 3),
         sd = signif(sd, 3))%>%
  rename('Mean R2' = mean_r2) %>%
  mutate('Response Variable' = case_when(response_variable == 'scale_cube_Mean_degree_decay_rate' ~ 'Kdd',
                                       response_variable == 'scale_cube_Mean_Decay_Rate_per_day' ~ 'Kcd')) %>%
  select('Response Variable', 'Mean R2', sd) %>%
  clipr::write_clip()

# ================================ create plots ===============================

## ==========  Decay vs ERs ==========

# Generating the 6-panel plot
pdf("Outputs/Decay_vs_ERs.pdf", width = 14, height = 9)

# Setting plotting parameters
par(pty = "s", mfrow = c(2, 3), oma = c(4, 4, 4, 4), mar = c(6, 6, 2, 1), mgp = c(4.5, 1, 0))

# Defining a function to create each plot
create_plot <- function(x, y, xlab, ylab, label) {
  mod_to_plot <- lm(as.formula(paste(y, "~", x)), data = cube_data)
  mod_sum <- summary(mod_to_plot)
  plot(cube_data[[x]], cube_data[[y]], xlab = xlab, ylab = ylab, cex.lab = 2.8, cex.axis = 2.0)
  
  # Add the regression line
  abline(mod_to_plot, lwd = 3)
  
  # Add R-squared and p-value
  mtext(bquote(R^2 == .(round(mod_sum$r.squared, digits = 2))), side = 3, adj = 0.85, line = -2, cex = 1.2)
  mtext(paste0("p = ", round(mod_sum$coefficients[2, 4], digits = 3)), side = 3, adj = 0.85, line = -3.5, cex = 1.2)
  
  # Add label
  mtext(label, cex = 2, side = 1, adj = 0.05, line = -1.5)
}

# Creating plots for per day decay rates
create_plot("cube_Total_Ecosystem_Respiration", "cube_Mean_Decay_Rate_per_day", "", expression((K[cd])^{1/3}), "A")
create_plot("cube_Sediment_Respiration", "cube_Mean_Decay_Rate_per_day", "", "", "B")
create_plot("cube_Water_Column_Respiration", "cube_Mean_Decay_Rate_per_day", "", "", "C")

# Creating plots for per degree day decay rates
create_plot("cube_Total_Ecosystem_Respiration", "cube_Mean_degree_decay_rate", expression(ER[tot]~(g~O[2]~m^-2~day^-1)^{1/3}), expression((K[dd])^{1/3}), "D")
create_plot("cube_Sediment_Respiration", "cube_Mean_degree_decay_rate", expression(ER[sed]~(g~O[2]~m^-2~day^-1)^{1/3}), "", "E")
create_plot("cube_Water_Column_Respiration", "cube_Mean_degree_decay_rate", expression(ER[wc]~(g~O[2]~m^-2~day^-1)^{1/3}), "", "F")

dev.off()



# Performing multiple regressions
day_mult_mod <- lm(cube_Mean_Decay_Rate_per_day ~ cube_Water_Column_Respiration + cube_Sediment_Respiration, data = cube_data)
day_mult_inter_mod <- lm(cube_Mean_Decay_Rate_per_day ~ cube_Water_Column_Respiration * cube_Sediment_Respiration, data = cube_data)
degree_mult_mod <- lm(cube_Mean_degree_decay_rate ~ cube_Water_Column_Respiration + cube_Sediment_Respiration, data = cube_data)
degree_mult_inter_mod <- lm(cube_Mean_degree_decay_rate ~ cube_Water_Column_Respiration * cube_Sediment_Respiration, data = cube_data)

day_ERtot_mod <- lm(cube_Mean_Decay_Rate_per_day ~ cube_Total_Ecosystem_Respiration, data = cube_data)
degree_ERtot_mod <- lm(cube_Mean_degree_decay_rate ~ cube_Total_Ecosystem_Respiration, data = cube_data)

# AIC values for model comparison
AIC(day_ERtot_mod)
AIC(day_mult_mod)
AIC(day_mult_inter_mod)

AIC(degree_ERtot_mod)
AIC(degree_mult_mod)
AIC(degree_mult_inter_mod)

## Plotting regression of rates vs. drainage area
pdf(file = "Outputs/Decay_vs_Drainage.pdf", height = 14, width = 10)

par(pty = "s", mfrow = c(2, 1), oma = c(4, 4, 4, 4), mar = c(5, 5, 2, 1), mgp = c(4, 0.75, 0))

# First plot with label "A"
mod_to_plot1 <- lm(cube_Mean_Decay_Rate_per_day ~ cube_totdasqkm, data = cube_data)
mod_sum1 <- summary(mod_to_plot1)
plot(cube_data$cube_totdasqkm, cube_data$cube_Mean_Decay_Rate_per_day, xlab = bquote(Drainage ~ Area ~ (km^-2)^{1/3}), ylab = expression((K[cd])^{1/3}), cex.lab = 3, cex.axis = 1.5)
abline(mod_to_plot1, lwd = 3)
mtext(bquote(R^2 == .(round(mod_sum1$r.squared, digits = 2))), side = 1, adj = 0.95, line = -3.5, cex = 2)
mtext(paste0("p = ", round(mod_sum1$coefficients[2, 4], digits = 3)), side = 1, adj = 0.95, line = -2, cex = 2)
mtext("A", cex = 3, side = 1, adj = 0.05, line = -1.5)

# Second plot with label "B"
mod_to_plot2 <- lm(cube_Mean_degree_decay_rate ~ cube_totdasqkm, data = cube_data)
mod_sum2 <- summary(mod_to_plot2)
plot(cube_data$cube_totdasqkm, cube_data$cube_Mean_degree_decay_rate, xlab = bquote(Drainage ~ Area ~ (km^-2)^{1/3}), ylab = expression((K[dd])^{1/3}), cex.lab = 3, cex.axis = 1.5)
abline(mod_to_plot2, lwd = 3)
mtext(bquote(R^2 == .(round(mod_sum2$r.squared, digits = 2))), side = 1, adj = 0.95, line = -3.5, cex = 2)
mtext(paste0("p = ", round(mod_sum2$coefficients[2, 4], digits = 3)), side = 1, adj = 0.95, line = -2, cex = 2)
mtext("B", cex = 3, side = 1, adj = 0.05, line = -1.5)

dev.off()


## ==========  Decay Rate Hists Scatters ==========

# function to generate scatter plot with adjacent histograms
# original version modified and pulled from https://www.r-bloggers.com/2011/06/example-8-41-scatterplot-with-marginal-histograms/
# some hard coding in this function to tweak the plot aesthetics

scatterhist = function(x, y, xlab.use="", ylab.use=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE,breaks=20)
  yhist = hist(y, plot=FALSE,breaks=20)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(6,6,1,1))
  plot(x,y,cex.axis=1.5,cex=1.5,xlab="",ylab="",cex.lab=2.5)
  points(lowess(y~x),typ="l",col=4,lwd=2)
  par(mar=c(0,5,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(5,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab.use, side=1, line=1, outer=TRUE, adj=0, cex = 2,
        at=1.75 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab.use, side=2, line=0.5, outer=TRUE, adj=0,cex = 2, 
        at=(1.5 * (mean(y) - min(y))/(max(y) - min(y))))
}

# generate the plot

pdf("Outputs/Decay_Scatter_Hists.pdf")
scatterhist(x = decay_temp$degree_decay_rate,y = decay_temp$Decay_Rate_per_day, xlab.use = expression(K[dd]), ylab.use = expression(K[cd]))
dev.off()

# regress decay per day against summed temperature

scatter.no.hist = function(x, y, xlab.use="", ylab.use=""){ # this keeps everything the same for sizing, but doesn't plot the histograms
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE,breaks=20)
  yhist = hist(y, plot=FALSE,breaks=20)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(6,6,1,1))
  plot(x,y,cex.axis=1.5,cex=1.5,xlab="",ylab="",cex.lab=2.5)
  points(lowess(y~x),typ="l",col=4,lwd=2)
  par(mar=c(0,5,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,plot=FALSE)
  par(mar=c(5,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE,plot=FALSE)
  par(oma=c(3,3,0,0))
  mtext(xlab.use, side=1, line=1, outer=TRUE, adj=0, cex = 2,
        at=.15 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab.use, side=2, line=0.5, outer=TRUE, adj=0,cex = 2, 
        at=(1.5 * (mean(y) - min(y))/(max(y) - min(y))))
}

pdf("Outputs/Decay_v_Temp.pdf")
scatter.no.hist(x = decay_temp$sum_mean_daily_temp,y = decay_temp$Decay_Rate_per_day,xlab.use = expression(Summed~Temperature~(degree*C)),ylab.use = expression(K[cd]))
dev.off()

## ==========  Hists: YRB vs global ==========

# compare statistical distributions of cotton strip decay rates between those observed in the YRB and those from a global dataset

global.cotton.Kcd = read_csv('https://github.com/dmcostello/CELLDEX2018/raw/refs/heads/master/str_k.csv')
range(global.cotton.Kcd$k)

global.cotton.Kdd = read_csv('https://github.com/dmcostello/CELLDEX2018/raw/refs/heads/master/str_k_dd.csv')
range(global.cotton.Kdd$k)


yrb.trim = decay_temp[,c("degree_decay_rate","Decay_Rate_per_day")]
range(yrb.trim$degree_decay_rate)
range(yrb.trim$Decay_Rate_per_day)

# generate kernels
global.Kcd.kern = density(global.cotton.Kcd$k,from = 0)
global.Kdd.kern = density(global.cotton.Kdd$k,from = 0)
yrb.Kdd.kern = density(yrb.trim$degree_decay_rate,from = 0)
yrb.Kcd.kern = density(yrb.trim$Decay_Rate_per_day,from = 0)

# make plots

pdf("Outputs/Global_YRB_Kern.pdf",width = 10,height = 5)

par(pty="s",mfrow=c(1,2))


plot(global.Kdd.kern,typ="l",main="",xlab = "Kdd",lwd=2,xlim=c(0,0.025),ylim=c(0,250),cex.lab=2,cex.axis=1.5)
points(yrb.Kdd.kern,typ="l",col=2,lwd=2)
legend(x = 0.012,y=250,col = c(1,2),lwd=2,legend = c('Global','YRB'),cex = 1.25)

plot(global.Kcd.kern,typ="l",main="",xlab = "Kcd",lwd=2,xlim=c(0,0.25),cex.lab=2,cex.axis=1.5)
points(yrb.Kcd.kern,typ="l",col=2,lwd=2)


dev.off()

