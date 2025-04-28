# ==============================================================================
#
# LASSO analysis for SSS Cotton Strips
#
# Status: In progress
#
# Note: Will need to download ERsed dp and change paths once published
#
# ==============================================================================
#
# Author: Brieanne Forbes 
# 22 April 2025
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

# =================================== user input ===============================
# toggle between response variables to run LASSO on each 

response_variable <- 'scale_cube_Mean_Decay_Rate_per_day'
# response_variable <- 'scale_cube_Mean_degree_decay_rate'

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



er_gpp <- read_csv('../SSS_metabolism/v2_SSS_Water_Sediment_Total_Respiration_GPP.csv',
                   comment = '#', na = '-9999') %>%
  select(Parent_ID, Site_ID, Sediment_Respiration, Total_Ecosystem_Respiration, Water_Column_Respiration, Gross_Primary_Production)

d50 <- read_csv('../SSS_metabolism/v2_SSS_ER_d50_TotalOxygenConsumed.csv', 
                comment = '#', na = '-9999') %>%
  select(Parent_ID, D50_m)

slope_vel_dis <- read_csv('../SSS_metabolism/Stream_Metabolizer/Inputs/v2_SSS_Slope_Discharge_Velocity.csv',
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
tss <- read_csv('./data/v3_SSS_Data_Package/Sample_Data/SSS_Water_TSS.csv',
                skip = 2, na = c('', 'N/A', '-9999')) %>%
  filter(!is.na(Sample_Name)) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         '00530_TSS_mg_per_L' = case_when(str_detect(`00530_TSS_mg_per_L`, 'LOD') ~ 0.12, # replace below LOD values with half LOD (LOD = 0.24)
                                          TRUE ~ as.numeric(`00530_TSS_mg_per_L`)))%>%
  select(Parent_ID, contains('TSS')) 

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
depth <- read_csv('./data/v3_SSS_Data_Package/v3_SSS_Water_Depth_Summary.csv',
                  comment = '#', na = c('', 'N/A', '-9999')) %>%
  select(Parent_ID, Average_Depth)

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
water_npoc_tn <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Water_NPOC_TN.csv',
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
sed_npoc_tn <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/CM_SSS_Sediment_NPOC_TN.csv',
                          skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         'Extractable_TN_mg_per_kg' = case_when(str_detect(`Extractable_TN_mg_per_kg`, 'Standard') ~ 0.05, # replace below standard values with half standard (standard = 0.1)
                                              TRUE ~ as.numeric(`Extractable_TN_mg_per_kg`))) %>%
  select(Parent_ID, contains('NPOC'), Extractable_TN_mg_per_kg) %>%
  group_by(Parent_ID) %>%
  summarise(Extractable_NPOC_mg_per_kg = round(mean(`Extractable_TN_mg_per_kg`), 2),
            Mean_Extractable_NPOC_mg_per_kg = round(mean(`Extractable_NPOC_mg_per_kg`), 2),
  ) %>%
  ungroup()


#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
cn <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/CM_SSS_Sediment_CN.csv',
                          skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         `01395_C_percent_per_mg` = as.numeric(`01395_C_percent_per_mg`),
         `01397_N_percent_per_mg` = as.numeric(`01397_N_percent_per_mg`)) %>%
  select(Parent_ID, `01395_C_percent_per_mg`, `01397_N_percent_per_mg`)

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
ions <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/CM_SSS_Water_Ions.csv',
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
resp <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/v2_CM_SSS_Sediment_Normalized_Respiration_Rates.csv',
                 skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}")) %>%
  group_by(Parent_ID) %>%
  summarise(Mean_Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment  = round(mean(as.numeric(Normalized_Respiration_Rate_mg_DO_per_H_per_L_sediment ), na.rm = T), 2)) %>%
  ungroup()


#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
iron <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Fe.csv',
                 skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}")) %>%
  group_by(Parent_ID) %>%
  summarise(Mean_Fe_mg_per_kg = round(mean(as.numeric(Fe_mg_per_kg), na.rm = T), 2)) %>%
  ungroup()

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
sand <- read_csv('./data/v5_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Grain_Size.csv',
               skip = 2, na = c('', 'N/A', '-9999'))%>%
  filter(!is.na(Sample_Name),
         str_detect(Sample_Name, 'SSS')) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         Percent_Tot_Sand = as.numeric(Percent_Tot_Sand)) %>%
  select(Parent_ID, Percent_Tot_Sand)

# =============================== cube function ===============================

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

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
# retruning 0 x 33 so no NA values 
all_data %>%
  filter(if_any(everything(), is.na))

# ======================= assess co-correlation ===============================

long_data <-  all_data %>% 
  pivot_longer(cols = -c(Site_ID, Parent_ID), names_to = "variable", values_to = "value")

ggplot() + 
  geom_histogram(long_data, mapping = aes(x = value)) + 
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()

## ======== Spearman correlation before transformations ============

spearman <- cor(all_data %>% select(-Site_ID, -Parent_ID), method = "spearman", use = "complete.obs")

# png(file = paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()),"_Scale_Spearman_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

# corrplot(spearman,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Spearman Correlation")

# dev.off()

spear.panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  r = (cor(x, y, method = c("spearman")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)}
  text(0.5, 0.5, txt, cex = cex.cor * (1 + abs(r))/2)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y, method = "spearman")
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex.cor, col=2)
  
}

panel.smooth <- function(x, y) {
  points(x, y, pch = 19, col = rgb(0.1, 0.2, 0.5, alpha = 0.3))
  abline(lm(y ~ x), col = 'blue', lty = 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1))
  
  h <- hist(x, plot = FALSE, breaks = "FD")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  
  rect(breaks[-nB], 0, breaks[-1], y, col="grey", border="white", ...)
}

# png(file = paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()),"_Pairs_Spearman_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

# pairs(all_data %>% select(-Site_ID, -Parent_ID),
#       lower.panel = panel.smooth, 
#       upper.panel = spear.panel.cor, 
#       diag.panel = panel.hist,
#       labels = colnames(all_data %>% select(-Site_ID, -Parent_ID)),
#       cex.labels = 0.8) 

# dev.off()

## ======== Pearson correlation before transformations ============
# function for pearson corr matrix

pear.panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y, method = c("pearson")))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) {cex.cor <- 0.8/strwidth(txt)} 
  text(0.5, 0.5, txt, cex = cex.cor * (1 + abs(r))/2)
  
  # if(missing(cex.cor)) {cex <- 1.2/strwidth(txt)} else {cex = cex.cor}
  # text(0.5, 0.5, txt, cex = cex * sin(sqrt(abs(r))))
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
  #text(0.5, 0.5, txt, cex = cex * r)
  text(.5, .8, Signif, cex=cex.cor, col=2)
  
}

# png(file = paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()),"_Pairs_Pearson_Correlation_Matrix.png"), width = 12, height = 12, units = "in", res = 300)

# pairs(all_data %>% select(-Site_ID, -Parent_ID),
#       lower.panel = panel.smooth, 
#       upper.panel = pear.panel.cor, 
#       diag.panel = panel.hist,
#       labels = colnames(all_data %>% select(-Site_ID, -Parent_ID)),
#       cex.labels = 0.8) 

# dev.off()

## ======== Cube root ======

cube_data <-  all_data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(where(is.numeric), .fn = ~ paste0("cube_", .x)) 

long_cube_data <-  cube_data %>%
  pivot_longer(cols = -c(Site_ID, Parent_ID), names_to = "variable", values_to = "value")

ggplot() +
  geom_histogram(long_cube_data, mapping = aes(x = value)) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()

### ======== Spearman correlation with cube transformation ============

spearman <- cor(cube_data %>% select(-Site_ID, -Parent_ID), method = "spearman", use = "complete.obs")

# png(file = paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()),"_Scale_Spearman_Correlation_Matrix_Cubed.png"), width = 12, height = 12, units = "in", res = 300)

# corrplot(spearman,type = "upper", method = "number", tl.col = "black", tl.cex = 1.6, cl.cex = 1.25,  title = "Spearman Correlation")

# dev.off()

# png(file = paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()),"_Pairs_Spearman_Correlation_Matrix_Cubed.png"), width = 12, height = 12, units = "in", res = 300)

# pairs(cube_data %>% select(-Site_ID, -Parent_ID),
#       lower.panel = panel.smooth, 
#       upper.panel = spear.panel.cor, 
#       diag.panel = panel.hist,
#       labels = colnames(cube_data %>% select(-Site_ID, -Parent_ID)),
#       cex.labels = 0.8) 

# dev.off()

### ======== Pearson correlation cube transformation ============
# function for pearson corr matrix

# png(file = paste0("./Figures/LASSO_Analysis/", as.character(Sys.Date()),"_Pairs_Pearson_Correlation_Matrix_Cubed.png"), width = 12, height = 12, units = "in", res = 300)
# 
# pairs(cube_data %>% select(-Site_ID, -Parent_ID),
#       lower.panel = panel.smooth, 
#       upper.panel = pear.panel.cor, 
#       diag.panel = panel.hist,
#       labels = colnames(cube_data %>% select(-Site_ID, -Parent_ID)),
#       cex.labels = 0.8) 
# 
# dev.off()

# ======== LASSO  ============

## scale data
scale_cube_variables = as.data.frame(scale(cube_data %>% select(-Parent_ID, -Site_ID)))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x))

## Loop through LASSO to get average over a lot of seeds ####

num_seeds = 100
seeds = sample(1:500, num_seeds)


## Set response variable (scale_cube_Mean_Decay_Rate_per_day/scale_cube_Mean_degree_decay_rate) and scale
yvar <- data.matrix(scale_cube_variables %>% pull(response_variable))
round(mean(yvar), 4)
sd(yvar)

# list for storing LASSO iterations
norm_coeffs = list()
lasso_coefs_pull = list()
r2_scores = numeric(num_seeds)

## Set predictor variables and scale
exclude_col = c("scale_cube_Mean_Decay_Rate_per_day", 'scale_cube_Mean_degree_decay_rate')

x_cube_variables = scale_cube_variables %>%
  select(-exclude_col)

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
colnames(lasso_coef_mat) = make.names(colnames(lasso_coef_mat), unique = T)
# Make DF of all LASSO results with mean and std. dev  
lasso_coef_means = lasso_coef_mat %>% 
  mutate(RowNames = rownames(lasso_coef_mat)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)

norm_coeffs_matrix = do.call(cbind, norm_coeffs)

mean_coeffs = as.data.frame(norm_coeffs_matrix, row.names = rownames(norm_coeffs_matrix))

colnames(mean_coeffs) = make.names(colnames(mean_coeffs), unique = T)

mean_coeffs_df = mean_coeffs %>% 
  mutate(RowNames = rownames(mean_coeffs)) %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(contains("s1"))), 
         sd = sd(c_across(contains("s1")))) %>% 
  relocate(mean, .before = s1) %>% 
  relocate(sd, .before = s1) %>% 
  relocate(RowNames, .before = mean)

results_r2 = as.data.frame(r2_scores) 
mean(results_r2$r2_scores)
sd(results_r2$r2_scores)

clipr::write_clip(mean_coeffs_df%>% 
  mutate(abs_mean = abs(mean),
         mean = signif(mean, 5) ,
         sd = signif(sd, 5)) %>% 
  arrange(desc(abs_mean))%>% 
    select(RowNames, mean, sd) ) 

clipr::write_clip(lasso_coef_means%>% 
                    mutate(abs_mean = abs(mean),
                           mean = signif(mean, 5) ,
                           sd = signif(sd, 5)) %>% 
                    arrange(desc(abs_mean))%>% 
                    select(RowNames, mean, sd) ) 

# ================================ create plots ===============================

## ==========  Decay vs ERs ==========

# Generating the 6-panel plot
pdf("Outputs/Decay_vs_ERs.pdf", width = 14, height = 9)

# Setting plotting parameters
par(pty = "s", mfrow = c(2, 3), oma = c(4, 4, 4, 4), mar = c(5, 5, 2, 1), mgp = c(3.5, 1, 0))

# Defining a function to create each plot
create_plot <- function(x, y, xlab, ylab, label) {
  mod_to_plot <- lm(as.formula(paste(y, "~", x)), data = cube_data)
  mod_sum <- summary(mod_to_plot)
  plot(cube_data[[x]], cube_data[[y]], xlab = xlab, ylab = ylab, cex.lab = 2, cex.axis = 1.5)
  
  # Add the regression line
  abline(mod_to_plot, lwd = 3)
  
  # Add R-squared and p-value
  mtext(bquote(R^2 == .(round(mod_sum$r.squared, digits = 2))), side = 3, adj = 0.85, line = -2)
  mtext(paste0("p = ", round(mod_sum$coefficients[2, 4], digits = 3)), side = 3, adj = 0.85, line = -3.5)
  
  
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

# Plotting regression of rates vs. drainage area
pdf(file = "Outputs/Decay_vs_Drainage.pdf", height = 14, width = 10)

par(pty = "s", mfrow = c(2, 1), oma = c(4, 4, 4, 4), mar = c(5, 5, 2, 1), mgp = c(3.5, 1, 0))

mod_to_plot1 <- lm(cube_Mean_Decay_Rate_per_day ~ cube_totdasqkm, data = cube_data)
mod_sum1 <- summary(mod_to_plot1)
plot(cube_data$cube_totdasqkm, cube_data$cube_Mean_Decay_Rate_per_day, xlab = bquote(Drainage ~ Area ~ (km^-2)^{1/3}), ylab = expression((K[cd])^{1/3}), cex.lab = 2, cex.axis = 1.5)
abline(mod_to_plot1, lwd = 3)
mtext(bquote(R^2 == .(round(mod_sum1$r.squared, digits = 2))), side = 1, adj = 0.95, line = -3.5)
mtext(paste0("p = ", round(mod_sum1$coefficients[2, 4], digits = 3)), side = 1, adj = 0.95, line = -2)
mtext("A", cex = 2, side = 1, adj = 0.05, line = -1.5)

mod_to_plot2 <- lm(cube_Mean_degree_decay_rate ~ cube_totdasqkm, data = cube_data)
mod_sum2 <- summary(mod_to_plot2)
plot(cube_data$cube_totdasqkm, cube_data$cube_Mean_degree_decay_rate, xlab = bquote(Drainage ~ Area ~ (km^-2)^{1/3}), ylab = expression((K[dd])^{1/3}), cex.lab = 2, cex.axis = 1.5)
abline(mod_to_plot2, lwd = 3)
mtext(bquote(R^2 == .(round(mod_sum2$r.squared, digits = 2))), side = 1, adj = 0.95, line = -3.5)
mtext(paste0("p = ", round(mod_sum2$coefficients[2, 4], digits = 3)), side = 1, adj = 0.95, line = -2)
mtext("B", cex = 2, side = 1, adj = 0.05, line = -1.5)

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
