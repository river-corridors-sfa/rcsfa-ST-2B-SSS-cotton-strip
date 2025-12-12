# ==============================================================================
#
# Make map for Cotton Strips
#
# Status: complete
#
# ==============================================================================
#
# Author: Brieanne Forbes, James Stegen
# 20 Dec 2024
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
library(rstudioapi)
library(viridis)


rm(list=ls(all=T))

# Setting wd 
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('./../')

# ================================= User inputs ================================

metadata_file <- './Published_Data/v3_SSS_Data_Package/v2_SSS_Field_Metadata.csv'

data_file <- './Outputs/Decay_Data.csv'

summary <-  read_csv("./Published_Data/v3_SSS_Data_Package/Sample_Data/SSS_CottonStrip_TensileStrength_DecayRate_Summary.csv", skip = 2)%>% 
  mutate(Sample_Name = str_extract(Sample_Name,"SSS[0-9]{3}")) %>% 
  filter(!is.na(Sample_Name)) %>% 
  dplyr::select(Sample_Name, Mean_Decay_Rate_per_day)


yrb_shp_dir <- './Published_Data/v2_SSS_Ecosystem_Respiration_Data_Package/Figures/Map_Layers/YakimaRiverBasin_Boundary/'

cluster_shp_dir <- './Published_Data/v2_SSS_Ecosystem_Respiration_Data_Package/Figures/Map_Layers/YRB_Cluster/'

common_crs = 4326

# ============================== read in and merge =============================

metadata <- read_csv(metadata_file) %>%
  dplyr::select(Parent_ID, Site_ID, Latitude, Longitude)

data <- read_csv(data_file) %>% 
  group_by(Parent_ID) %>% 
  mutate(mean_Decay_Rate_per_day = mean(Decay_Rate_per_day))

data <- data %>%
  dplyr::select(Parent_ID, mean_degree_decay_rate,mean_Decay_Rate_per_day) %>% 
  distinct(Parent_ID, .keep_all = TRUE)

merge <- data %>%
  left_join(metadata, by = c('Parent_ID' = 'Parent_ID')) 
  
# ============================ read in YRB shp file ============================

YRB_shp <- list.files(yrb_shp_dir, 'Basin.shp', full.names = T)

YRB_boundary <- read_sf(YRB_shp) %>%
  st_transform(common_crs)

# ============================ convert to sf object ============================

sites <- st_as_sf(merge, coords = c('Longitude','Latitude'), crs = common_crs)

# ======================== pull NHD data and elevation =========================

YRB_flowlines <- get_nhdplus(AOI = YRB_boundary$geometry, streamorder = 3)

elevation_raw <- get_elev_raster(YRB_boundary, z = 10)

elevation_crop <- mask(elevation_raw, YRB_boundary)

elevation <- as.data.frame(elevation_crop, xy = T) %>%
  as_tibble() %>%
  rename("long" = x,
         "lat" = y,
         "elevation" = 3) %>% #column index > name (changing resolution changes colname)
  filter(!is.na(elevation))

# ============================ read in cluster shp file ========================

cluster_shp <- list.files(cluster_shp_dir, 'shp', full.names = T)

cluster <- read_sf(cluster_shp) %>%
  st_transform(common_crs)


# ========================= create insert map ======================

data("us_states", package = "spData")
us_states_4326 = st_transform(us_states, crs = 4326)

wa <- us_states_4326 %>% filter(NAME == "Washington")

insert <- ggplot() +
  geom_sf(data = us_states_4326, fill = "white") + 
  geom_sf(data = wa, fill = "black",colour = "black")+
  geom_sf(data = YRB_boundary, colour = "red", fill = 'red') +
  labs(x = "", y = "")+
  theme_map()

# ========================= create map of Decay Rate per Degree Day (elevation) ======================

kdd_map <- ggplot()+
  geom_sf(data = YRB_boundary, fill = "white")+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.6)+
  new_scale_fill()+
  geom_sf(data = sites, aes(color = mean_degree_decay_rate, size = mean_degree_decay_rate), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.2, direction = -1)+
  scale_color_viridis(option = 'B', begin = 0.2, direction = -1)+
  scale_size(range = c(1.5, 6))+
  theme_map() +
  labs(x = "", y = "", color = bquote(K[dd])) +
  guides(size = "none")+
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"),
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(1.5, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20")) +
  theme(legend.position = c(0.75,0.7))

full <- ggdraw() +
  draw_plot(kdd_map) 

ggsave('./Outputs/SPS_kdd_Map.pdf',
       full,
       width = 8,
       height = 5
)

rm('kdd_map')

# ========================= create map of Decay Rate per Day (elevation) ======================

kcd_map <- ggplot()+
  geom_sf(data = YRB_boundary, fill = "white")+
  geom_raster(data = elevation, aes(long, lat, fill = elevation), show.legend = F, alpha = 0.4)+
  scale_fill_gradient(low = 'white', high = 'black')+
  geom_sf(data = YRB_flowlines, color = "royalblue", alpha = 0.6)+
  new_scale_fill()+
  geom_sf(data = sites, aes(color = mean_Decay_Rate_per_day, size = mean_Decay_Rate_per_day), show.legend = T) +
  scale_fill_viridis(option = 'B', begin = 0.2, direction = -1)+
  scale_color_viridis(option = 'B', begin = 0.2, direction = -1)+
  scale_size(range = c(1.5, 6))+
  theme_map() +
  labs(x = "", y = "", color = bquote(K[cd])) +
  guides(size = "none")+
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"),
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(1.5, "in"),
    # pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20")) +
  theme(legend.position = c(0.75,0.7))

full <- ggdraw() +
  draw_plot(kcd_map) 

ggsave('./Outputs/SPS_kcd_Map.pdf',
       full,
       width = 8,
       height = 5
)


# ========================= create map using CRB cluster and show field sites ======================

Cotton_map_cluster <- ggplot()+
  geom_sf(data = YRB_boundary)+
  geom_sf(data = cluster, aes(fill = as.factor(ClusterNum), color = as.factor(ClusterNum)), show.legend = T)+
  scale_fill_manual(values = c('#1a9850', 'steelblue1', '#91cf60', '#8c510a', '#d9ef8b', '#f6e8c3'))+
  scale_color_manual(values = c('#1a9850', 'steelblue1', '#91cf60', '#8c510a', '#d9ef8b', '#f6e8c3'))+
  geom_sf(data = YRB_flowlines, color = "royalblue")+
  new_scale_fill()+
  new_scale_color()+
  geom_sf(data = sites, show.legend = F, size = 3) +
  theme_map() +
  ggspatial::annotation_scale(
    location = "br",
    pad_x = unit(0.5, "in"),
    bar_cols = c("black", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(1.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20"))

full_cluster <- ggdraw() +
  draw_plot(Cotton_map_cluster) +
  draw_plot(insert, x = 0.4, y = 0.4, width = 0.3, height = 0.3)

ggsave('./Outputs/Cotton_Map_Cluster.pdf',
       full_cluster,
       width = 8,
       height = 5
)



