# ==============================================================================
#
# Maps for ERsed manuscript 
#
# Status: Complete 
# 
# ==============================================================================
#
# Author: Brieanne Forbes 
# 23 January 2025
#
# ==============================================================================
library(tidyverse) #keep it tidy
library(raster) # work with rasters, NOTE: masks dplyr::select
library(janitor) # clean_names()
library(ggthemes) # theme_map()
library(ggnewscale) # set multiple color scales
library(ggspatial) # add north arrow and scale bar
library(nhdplusTools) # get watershed boundary/flowlines
library(elevatr) # pull elevation maps
library(sf) # tidy spatial
library(spData)
library(cowplot)
library(viridis)

rm(list=ls(all=T))

# Setting wd to parent folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./..")

# ================================= User inputs ================================

# downloaded from https://data.ess-dive.lbl.gov/view/doi%3A10.15485%2F1971251
metadata_file <- './Published_Data/v3_RCSFA_Geospatial_Data_Package/v3_RCSFA_Geospatial_Site_Information.csv'

data_file <- './v2_SSS_Water_Sediment_Total_Respiration_GPP.csv'

shp_dir <- './Figures/Map_Layers/YakimaRiverBasin_Boundary'

modelled_ER <- './v2_SSS_ER_d50_TotalOxygenConsumed.csv'

common_crs = 4326

# ============================== read in and merge =============================

metadata <- read_csv(metadata_file) %>%
  dplyr::select(Site_ID, Latitude, Longitude)

data <- read_csv(data_file, comment = '#', na = '-9999') %>%
  mutate(Total_Ecosystem_Respiration = case_when(Total_Ecosystem_Respiration > 0 ~ NA,
                                                        TRUE ~ Total_Ecosystem_Respiration))

merge <- data %>%
  left_join(metadata)
  
# ============================ read in YRB shp file ============================

YRB_shp <- list.files(shp_dir, 'shp', full.names = T)

YRB_boundary <- read_sf(YRB_shp) %>%
  st_transform(common_crs)

# ============================ convert to sf object ============================

sites <- st_as_sf(merge, coords = c('Longitude','Latitude'), crs = common_crs)

# ======================== pull NHD data and elevation =========================

YRB_flowlines <- get_nhdplus(AOI = YRB_boundary$geometry, streamorder = 3)

elevation_raw <- get_elev_raster(YRB_boundary$geometry, z = 10)

elevation_crop <- mask(elevation_raw, YRB_boundary)

elevation <- as.data.frame(elevation_crop, xy = T) %>% 
  as_tibble() %>% 
  rename("long" = x, 
         "lat" = y, 
         "elevation" = 3) %>% #column index -> name (changing resolution changes colname)
  filter(!is.na(elevation))


# ======================== create map of observed ER Tot =======================

merge_ER <- data %>% 
  full_join(read_csv(modelled_ER, comment = '#', na = '-9999')) %>%
  left_join(metadata) %>%
  arrange(Total_Ecosystem_Respiration)

ER_sf <- merge_ER %>% 
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_tot_obs_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf, aes(color = Total_Ecosystem_Respiration, size = Total_Ecosystem_Respiration), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3)+
  scale_color_viridis(option = 'B', begin = 0.3)+ 
  scale_size(range = c(3, 8), trans = 'reverse')+
  new_scale_color()+
  geom_sf(data = ER_sf %>% filter(is.na(Total_Ecosystem_Respiration)), aes(color = 'grey'), size = 2.5, show.legend = T) +
  scale_color_manual(values = c("grey" = "grey60")) +
  theme_map() + 
  labs(x = "", y = "", color = "Total Ecosystem\nRespiration\n(g O2 m2 day-1)") + 
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"), 
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))


ggsave('./Figures/Intermediate_Files/ERtot_Map.pdf',
       ER_tot_obs_map,
       width = 10,
       height = 5
)

# ======================= create map of predicted ER Hz =======================

merge_ER <- merge_ER %>%
  arrange(Total_Oxygen_Consumed_g_per_m2_per_day)

ER_sf <- merge_ER %>% 
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_tot_pred_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf, aes(color = Total_Oxygen_Consumed_g_per_m2_per_day, size = Total_Oxygen_Consumed_g_per_m2_per_day), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3)+
  scale_color_viridis(option = 'B', begin = 0.3)+
  scale_size(range = c(2, 8), trans = 'reverse')+
  theme_map() + 
  labs(x = "", y = "", color = "Total Oxygen\nConsumed\n(g O2 m2 day-1)") + 
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"), 
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))


ggsave('./Figures/Intermediate_Files/ERhz_Map.pdf',
       ER_tot_pred_map,
       width = 10,
       height = 5
)

# =============== create map of predicted ER Hz (z scores) =====================

merge_ER <- merge_ER %>%
  mutate(Total_Oxygen_Consumed_g_per_m2_per_day_Z = c(scale(Total_Oxygen_Consumed_g_per_m2_per_day, center = TRUE, scale = TRUE))) %>%
  arrange(Total_Oxygen_Consumed_g_per_m2_per_day_Z) %>%
  filter(!is.na(Total_Ecosystem_Respiration))

ER_sf <- merge_ER %>%
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_tot_pred_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf, aes(color = Total_Oxygen_Consumed_g_per_m2_per_day_Z, size = Total_Oxygen_Consumed_g_per_m2_per_day_Z), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3, limits = c(-4.5, 1))+
  scale_color_viridis(option = 'B', begin = 0.3, limits = c(-4.5, 1))+
  scale_radius(trans = 'reverse', range = c(2, 6))+
  theme_map() +
  labs(x = "", y = "", color = "Normalized Total\nOxygen Consumed") +
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"),
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))


ggsave('./Figures/Intermediate_Files/ERhz_Map_ZScores.pdf',
       ER_tot_pred_map,
       width = 10,
       height = 5
)

# ======================= create map of observed ER sed =======================

merge_ER <- merge_ER %>% 
  arrange(Sediment_Respiration)

ER_sf <- merge_ER %>% 
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_sed_obs_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf, aes(color = Sediment_Respiration, size = Sediment_Respiration), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3)+
  scale_color_viridis(option = 'B', begin = 0.3)+ 
  scale_size(range = c(3, 8), trans = 'reverse')+
  new_scale_color()+
  geom_sf(data = ER_sf %>% filter(is.na(Sediment_Respiration)), aes(color = 'grey'), size = 2.5, show.legend = T) +
  scale_color_manual(values = c("grey" = "grey60")) +
  theme_map() + 
  labs(x = "", y = "", color = "Sediment Respiration\n(g O2 m2 day-1)") + 
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"), 
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))


ggsave('./Figures/Intermediate_Files/ERsed_Map.pdf',
       ER_sed_obs_map,
       width = 10,
       height = 5
)

# =================== create map of observed ER Sed (z scores) =================

merge_ER <- merge_ER %>%
  mutate(Sediment_Respiration_Z = c(scale(Sediment_Respiration, center = TRUE, scale = TRUE))) %>%
  arrange(Sediment_Respiration_Z) %>%
  filter(!is.na(Total_Ecosystem_Respiration))

ER_sf <- merge_ER %>%
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_sed_obs_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf, aes(color = Sediment_Respiration_Z, size = Sediment_Respiration_Z), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3, limits = c(-4.5, 1))+
  scale_color_viridis(option = 'B', begin = 0.3, limits = c(-4.5, 1))+
  scale_radius(range = c(2, 8), trans = 'reverse')+
  theme_map() +
  labs(x = "", y = "", color = "Normalized Sediment\nEcosystem Respiration") +
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"),
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))


ggsave('./Figures/Intermediate_Files/ERsed_Map_ZScores.pdf',
       ER_sed_obs_map,
       width = 10,
       height = 5
)

# ======================= create map of observed ER wc =======================
 
merge_ER <- merge_ER %>%
  arrange(Water_Column_Respiration)

ER_sf <- merge_ER %>% 
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_wc_obs_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf, aes(color = Water_Column_Respiration, size = Water_Column_Respiration), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3)+
  scale_color_viridis(option = 'B', begin = 0.3)+ 
  scale_size(range = c(3, 8), trans = 'reverse')+
  theme_map() + 
  labs(x = "", y = "", color = "Water Column\nRespiration\n(g O2 m2 day-1)") + 
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"), 
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))

ggsave('./Figures/Intermediate_Files/ERwc_Map.pdf',
       ER_wc_obs_map,
       width = 10,
       height = 5
)

# ======================= create map of ERsed contribution =====================

merge_ER <-  merge_ER %>%
  rename(ERsed_Contribution = Sediment_Respiration_Contribution) %>%
  arrange(ERsed_Contribution)

ER_sf <- merge_ER %>% 
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_sed_contrib_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf %>% filter(!is.na(Total_Ecosystem_Respiration)), aes(color = ERsed_Contribution, size = Sediment_Respiration), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3, direction = -1)+
  scale_color_viridis(option = 'B', begin = 0.3, direction = -1)+ 
  scale_size(range = c(3, 8), trans = 'reverse')+
  theme_map() + 
  labs(x = "", y = "", color = "ERsed Fractional\nContribution to ERtot") + 
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"), 
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))

ggsave('./Figures/Intermediate_Files/ERsed_Contribution.pdf',
       ER_sed_contrib_map,
       width = 10,
       height = 5
)

# ======================= create map of ERwc contribution =====================

merge_ER <-  merge_ER %>%
  rename(ERwc_Contribution = Water_Column_Respiration_Contribution) %>%
  arrange(ERwc_Contribution)

ER_sf <- merge_ER %>% 
  st_as_sf(coords = c('Longitude','Latitude'), crs = common_crs)

ER_wc_contrib_map <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.8)+
  new_scale_fill()+
  geom_sf(data = ER_sf %>% filter(!is.na(Total_Ecosystem_Respiration)), aes(color = ERwc_Contribution, size = Water_Column_Respiration), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.3, direction = -1)+
  scale_color_viridis(option = 'B', begin = 0.3, direction = -1)+ 
  scale_size(range = c(3, 8), trans = 'reverse')+
  theme_map() + 
  labs(x = "", y = "", color = "ERwc Fractional\nContribution to ERtot") + 
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"), 
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(2, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))

ggsave('./Figures/Intermediate_Files/ERwc_Contribution.pdf',
       ER_wc_contrib_map,
       width = 10,
       height = 5
)

