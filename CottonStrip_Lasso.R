# ==============================================================================
#
# LASSO analysis for SSS Cotton Strips
#
# Status: In progress
#
# Note: Lines saving investigation plots are commented out 
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
setwd("./..")

# =================================== find files ===============================

er_gpp <- read_csv('./v2_SSS_Water_Sediment_Total_Respiration_GPP.csv',
                   comment = '#', na = '-9999') %>%
  select(Parent_ID, Site_ID, Sediment_Respiration, Gross_Primary_Production)

# NHD+, streamcat, NLCD, ET0 extracted geospatial variables https://github.com/river-corridors-sfa/Geospatial_variables
geospatial <- read_csv('https://github.com/river-corridors-sfa/Geospatial_variables/raw/refs/heads/main/v4_RCSFA_Extracted_Geospatial_Data_2025-01-31.csv') %>%
  select(site, totdasqkm, pctmxfst2019ws,pctconif2019ws,pctdecid2019ws, AridityWs, pctcrop2019ws, pcthay2019ws, pctshrb2019ws) %>%
  filter(site %in% er_gpp$Site_ID)%>%
  mutate(PctFst = pctmxfst2019ws + pctdecid2019ws + pctconif2019ws,
         PctAg = pctcrop2019ws + pcthay2019ws) %>%
  rename(Site_ID = site)%>%
  select(Site_ID, totdasqkm, PctFst, AridityWs, PctAg, pctshrb2019ws)

d50 <- read_csv('./v2_SSS_ER_d50_TotalOxygenConsumed.csv', 
                comment = '#', na = '-9999') %>%
  select(Parent_ID, D50_m)

slope_vel_dis <- read_csv('./Stream_Metabolizer/Inputs/v2_SSS_Slope_Discharge_Velocity.csv',
                          comment = '#', na = '-9999') %>%
  select(Site_ID, Slope, Discharge, Velocity)

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
tss <- read_csv('./Published_Data/v3_SSS_Data_Package/Sample_Data/SSS_Water_TSS.csv',
                skip = 2, na = c('', 'N/A', '-9999')) %>%
  filter(!is.na(Sample_Name)) %>%
  mutate(Parent_ID = str_extract(Sample_Name, "^.{1,6}"),
         '00530_TSS_mg_per_L' = case_when(str_detect(`00530_TSS_mg_per_L`, 'LOD') ~ 0.12, # replace below LOD values with half LOD (LOD = 0.24)
                                          TRUE ~ as.numeric(`00530_TSS_mg_per_L`)))%>%
  select(Parent_ID, contains('TSS')) 

#downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1923689
npoc_tn <- read_csv('./Published_Data/v4_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Water_NPOC_TN.csv',
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

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
depth <- read_csv('./Published_Data/v3_SSS_Data_Package/v3_SSS_Water_Depth_Summary.csv',
                  comment = '#', na = c('', 'N/A', '-9999')) %>%
  select(Parent_ID, Average_Depth)

# downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566
hobo_temp <- read_csv('./Published_Data/v3_SSS_Data_Package/Sensor_Data/DepthHOBO/Plots_and_Summary_Statistics/v3_SSS_Water_Press_Temp_Summary.csv',
                      comment = '#', na = c('', 'N/A', '-9999')) %>%
  select(Parent_ID, Temperature_Mean)

# =============================== cube function ===============================

cube_root <- function(x) sign(x) * (abs(x))^(1/3)

# =============================== combine data ===============================

all_data <- er_gpp %>%
  full_join(geospatial, by = 'Site_ID')%>%
  full_join(slope_vel_dis, by = 'Site_ID')%>%
  full_join(d50, by = 'Parent_ID')%>%
  full_join(tss, by = 'Parent_ID')%>%
  full_join(npoc_tn, by = 'Parent_ID')%>%
  full_join(depth, by = 'Parent_ID')%>%
  full_join(hobo_temp, by = 'Parent_ID') %>%
  filter(!is.na(Sediment_Respiration)) 

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

## ===== run function to automatically determine best variable ====
cube_pearson <- cor(cube_data %>% select(-Site_ID, -Parent_ID), method = "pearson")

pearson_df <- as.data.frame(cube_pearson)

row_names_pearson <- rownames(pearson_df)

pearson_df$Variable <- row_names_pearson

pearson_melted <- reshape2::melt(pearson_df, id.vars = "Variable") %>% 
  filter(value != 1) %>% # remove everything correlated with self
  mutate(value = abs(value)) %>% # do this so it removes in order, and doesn't leave out high negative correlations
  filter(!grepl("cube_Sediment_Respiration", Variable)) # %>% # remove ERsed, don't want it to be removed 

# pull out ersed correlations only
ersed_melted <- pearson_melted %>% 
  filter(grepl("cube_Sediment_Respiration", variable)) 

choose_melted <- pearson_melted %>% 
  filter(!grepl("cube_Sediment_Respiration", variable)) %>%
  left_join(ersed_melted, by = "Variable") %>% 
  rename(Variable_1 = Variable) %>% 
  rename(Variable_2 = variable.x) %>% 
  rename(Correlation = value.x) %>% 
  rename(Variable_1_ERsed_Correlation = value.y) %>% 
  select(-c(variable.y)) %>% 
  left_join(ersed_melted, by = c("Variable_2" = "Variable")) %>% 
  rename(Variable_2_ERsed_Correlation = value) %>% 
  select(-c(variable))

loop_melt = choose_melted %>% 
  arrange(desc(Correlation))

#Pearson correlation coefficient to remove above
correlation = 0.7

## Start loop to remove highly correlated (> 0.7)
ersed_filter = function(loop_melt) {
  
  rows_to_keep = rep(TRUE, nrow(loop_melt))
  
  for (i in seq_len(nrow(loop_melt))) {
    
    if (!rows_to_keep[i]) next
    
    row = loop_melt[i, ]
    
    if (row$Correlation < correlation) next
    
    if(row$Variable_1_ERsed_Correlation >= row$Variable_2_ERsed_Correlation) {
      
      var_to_keep = row$Variable_1
      var_to_remove = row$Variable_2
      
    } else {
      
      var_to_keep = row$Variable_2
      var_to_remove = row$Variable_1
      
    }
    
    loop_melt$Variable_to_Keep[i] = var_to_keep
    loop_melt$Variable_to_Remove[i] = var_to_remove
    
    for (j in seq(i + 1, nrow(loop_melt))) {
      
      if(loop_melt$Variable_1[j] == var_to_remove || loop_melt$Variable_2[j] == var_to_remove) {
        
        rows_to_keep[j] = FALSE
        
      }
      
    }
    
    
  }
  
  return(loop_melt[rows_to_keep, ])
  
}

filtered_data <-  ersed_filter(loop_melt) 

# pull out variables to remove
removed_variables <-  filtered_data %>% 
  distinct(Variable_to_Remove)

# pull out all variables 
all_variables <-  ersed_melted %>% 
  select(c(Variable))

# remove variables from all variables to get variables to keep for LASSO 
kept_variables <- ersed_melted %>%
  filter(!Variable %in% removed_variables$Variable) %>%
  pull(Variable) %>%
  unique() 

col_to_keep = c(kept_variables, "cube_Sediment_Respiration")

cube_variables = cube_data %>%
  # select(all_of(col_to_keep)) # decided to not remove any variables 
  select(-Parent_ID, -Site_ID) # remove character columns before lasso

# ======== LASSO  ============

## scale data
scale_cube_variables = as.data.frame(scale(cube_variables))%>% 
  rename_with(where(is.numeric), .fn = ~ paste0("scale_", .x))

## Loop through LASSO to get average over a lot of seeds ####

num_seeds = 100
seeds = sample(1:500, num_seeds)

## Set response variable (cube_Sediment_Respiration) and scale
yvar <- data.matrix(scale_cube_variables$scale_cube_Sediment_Respiration)
round(mean(yvar), 4)
sd(yvar)

# list for storing LASSO iterations
norm_coeffs = list()
lasso_coefs_pull = list()
r2_scores = numeric(num_seeds)

## Set predictor variables and scale
exclude_col = "scale_cube_Sediment_Respiration"

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




