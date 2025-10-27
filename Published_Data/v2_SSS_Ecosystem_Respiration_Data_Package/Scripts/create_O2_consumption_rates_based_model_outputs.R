# Initilize parameter 
rm(list=ls())

library(sf)
library(tidyverse)

## calculating the O2 consumption rate (gO2/m2/day) from river corridor model outputs

## the model's raw model outputs are O2 consumption rate (mole) and NO2 and NO3 consumption rate
## aerobic respiration's reaction network f1=0.65*1/3
## anaerobic respiration reaction network f2=0.65 (NO3->NO2), f3=0.99 (NO2-> N2)
## O2 consumption rate (O2 mole)/f1*32g/365day/m2(surface area)

# Author: Kyongho Son (kyongho.son@pnnl.gov)

#published data downloaded from https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1962818


############  set wd 

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./..")

########## stream surface area 

nexss<-read_csv('./Published_Data/Son_et_al_2022_Respiration_Data_Package/model_inputs/nexss_inputs.csv')
nexss$area_m2<-nexss$length_m*10^(nexss$logw_m)

############ annual/seasonal substrate data year data-----------

#vertical direction # # actual N and DOC for vertical #
annual_10vert<-read_table('./Published_Data/Son_et_al_2022_Respiration_Data_Package/model_outputs/vert_ann.dat')

annual_10lat<-read_table('./Published_Data/Son_et_al_2022_Respiration_Data_Package/model_outputs/lat_annual.dat')

# extract 2nd and 3rd year cumulative respiration amounts (moles)
# model stabilizes after 2-3 year
tmp=subset(annual_10vert,annual_10vert$day==730)
tmp2=subset(annual_10vert,annual_10vert$day==1095)

# compute net annual  cumulative respiration---------

#subtract year 2 from year 3 to get respiration for just year 3
tmp3<-as.data.frame(cbind(tmp$COMID,tmp2$ver_o2_cons_mol-tmp$ver_o2_cons_mol,tmp2$ver_no3_prod_mol-tmp$ver_no3_prod_mol,tmp2$ver_no2_prod_mol-tmp$ver_no2_prod_mol))
colnames(tmp3)=c("COMID","ver_o2_cons_mole","ver_no3_cons_mole","ver_no2_cons_mole")
tmp3<-tmp3[,c(1,2)]
ver_annual_hr_CR<-tmp3

# repeat for lat
# extract 2nd and 3rd year cumulative respiration amounts (moles)
tmp=subset(annual_10lat,annual_10lat$day==730)
tmp2=subset(annual_10lat,annual_10lat$day==1095)

tmp3<-as.data.frame(cbind(tmp$COMID,tmp2$lat_o2_cons_mol-tmp$lat_o2_cons_mol,tmp2$lat_no3_prod_mol-tmp$lat_no3_prod_mol,tmp2$lat_no2_prod_mol-tmp$lat_no2_prod_mol))
colnames(tmp3)=c("COMID","lat_o2_cons_mole","lat_no3_cons_mole","lat_no2_cons_mole")
tmp3<-tmp3[,c(1,2)]
lat_annual_hr_CR<-tmp3

#### merge lateral and vertical O2 consumption (moles)
resp_hr_CR_annual_o2_consump_mole<-merge(ver_annual_hr_CR,lat_annual_hr_CR,by.x="COMID",by.y="COMID")

nhd_CR_stream_annual_o2_consum<-merge(nhd_CR_poly,resp_hr_CR_annual_o2_consump_mole, by.x = "COMID",by.y="COMID")
nhd_CR_stream_annual_o2_consum$tot_o2_cons_mole<-nhd_CR_stream_annual_o2_consum$ver_o2_cons_mole+nhd_CR_stream_annual_o2_consum$lat_o2_cons_mole
nhd_CR_stream_annual_o2_consum$tot_o2_cons_mole_day<-nhd_CR_stream_annual_o2_consum$tot_o2_cons_mole/365 #convert form annual to per day

nhd_CR_stream_annual_o2_consum<-merge(nhd_CR_stream_annual_o2_consum,nexss[,c("comid_nhd","area_m2")],by.x="COMID",by.y="comid_nhd") # merge with Nexss data to get surface area

nhd_CR_stream_annual_o2_consum$tot_o2_cons_g_m2_day<-nhd_CR_stream_annual_o2_consum$tot_o2_cons_mole_day*32/nhd_CR_stream_annual_o2_consum$area_m2

## Yakima river basin
nhd_CR_stream_annual_o2_consum_df<-data.frame(nhd_CR_stream_annual_o2_consum)
nhd_CR_stream_annual_o2_consum_df$geometry<-NULL


# the output of this script can be found within v2_SSS_ER_d50_TotalOxygenConsumed.csv as "Total_Oxygen_Consumed_g_per_m2_per_day"

