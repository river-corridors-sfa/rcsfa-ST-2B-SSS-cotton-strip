# ==============================================================================

# Create input file for stream metabolizer using data from SSS data package. 
# Published at https://data.ess-dive.lbl.gov/datasets/doi:10.15485/1969566.

# Status: Complete

# ==============================================================================

# Author: Brieanne Forbes (brieanne.forbes@pnnl.gov)
# 24 October 2024

# remove all files
rm(list = ls(all = TRUE))

# ==============================================================================

library(tidyverse)
library(plotly)
library(scales)
library(tsrobprep)
library(crayon)

# ================================= set wd ================================

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./..")

# ============================ find files =============================

DO_files <- list.files('./Published_Data/v3_SSS_Data_Package/Sensor_Data/miniDOT/Data', '.csv', full.names = T)

metadata <- './Published_Data/v3_SSS_Data_Package/v2_SSS_Field_Metadata.csv' %>% 
 read_csv()
  
baro_files <- list.files('./Published_Data/v3_SSS_Data_Package/Sensor_Data/BarotrollAtm/Data', '.csv', full.names = T)

hobo_files <- list.files('./Published_Data/v3_SSS_Data_Package/Sensor_Data/DepthHOBO/Data', '.csv', full.names = T)

depth_summary <- './Published_Data/v3_SSS_Data_Package/v3_SSS_Water_Depth_Summary.csv' %>%
  read_csv(comment = '#')%>%
  mutate(Date = ymd(Date),
         datetime = as_datetime(paste(Date, Start_Time)))

# ============================ set plot theme =============================

theme_set(
  theme(
    text = element_text(family = 'serif', face = 'plain'),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    line = element_line(size = 0.05),
    axis.line = element_line(size = 0.5),
    panel.background = element_rect(color = 'white'),
    panel.border = element_rect(
      colour = 'black',
      fill = NA,
      size = 0.5
    ),
    plot.title = element_text(size = 20, face = 'bold'),
    axis.ticks.length = unit(.25, 'cm'),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, family = 'serif'),
    strip.placement = "outside",
    plot.subtitle = element_text(size = 14, margin = margin(b = 10)),
    plot.margin = margin(0, 1, 0, 0, "cm")
  )
)

# ============================ make empty tibble =============================

hobo_time_difference_combine <- tibble(Parent_ID = as.character(),
                                       transect_depth_datetime = ymd_hms(),
                                       hobo_reference_datetime = ymd_hms())


# ============================ loop through files =============================

for (file in DO_files) {
  
  data <- read_csv(file, comment = '#', na = c('', '-9999',-9999, NA, 'N/A')) %>%
    select(-Dissolved_Oxygen_Saturation, -Battery)
  
  parent_ID <- str_extract(file, "[A-Z]{3}\\d{3}")
  
  # subset to 15 mins
  subset <- data %>%
    filter(minute(DateTime) %% 15 == 0)
  
  ## ============================ tsrobprep cleaning =============================
  
  set.seed(7)
  
  auto_clean <- auto_data_cleaning(
    data = subset$Dissolved_Oxygen, # Dissolved Oxygen data
    S = 96, # 96 intervals for daily seasonality (96 rows of 15 min data = 1 day of data)
    tau = NULL, # Performs lasso to determine tau
    no.of.last.indices.to.fix = nrow(subset), # Fix all data points
    indices.to.fix = NULL, # Automatically fix indices
    detect.outliers.pars = list(method = c("IQR"),  # Use IQR for outlier detection
                                threshold = c(1.5)  # Common threshold for IQR
                                ))
    
    clean_DO <- subset %>%
      add_column(Cleaned_DO = auto_clean$clean.data[,1, drop = T])
  
    
    ## ============================ fix sample day ==============================
    # use metadata to revert last two and first two points before and after sample day back to original data
    
    site_metadata <- metadata %>%
      filter(Parent_ID == parent_ID)
    
    sample_day <- site_metadata$Sample_Date
    
    datetime_to_fix <- c(as_datetime(paste0(sample_day-(1), '23:30:00')), # second to last point before sample day
                         as_datetime(paste0(sample_day-(1), '23:45:00')), # last point before sample day
                         as_datetime(paste0(sample_day+(1), '00:00:00')), # first point after sample day
                         as_datetime(paste0(sample_day+(1), '00:15:00'))) # second point after sample day
    
    clean_DO <- clean_DO %>%
      mutate(Dissolved_Oxygen_final = case_when(DateTime %in% datetime_to_fix ~ Dissolved_Oxygen,
                                                 TRUE ~ Cleaned_DO))
    
    
    ## ============================== check data ================================
    # assuming the cleaning algorithm did not work when cleaned values are more than 20 mg/L
    # in this case, use original data  
    # 17 mg/L is the reported maximum in the DO summary file, rounding up to 20 for the threshold
    
    if(max(clean_DO$Cleaned_DO) > 20){
      
      clean_DO <- clean_DO %>%
        select(DateTime, Parent_ID, Site_ID, Temperature, Dissolved_Oxygen)
      
      
      cat(red$bold(str_c(parent_ID, ' has data >20 mg/L. Removing cleaned DO and using original DO.')), "\n")
      
    } else {
      
      clean_DO <- clean_DO %>%
        select(DateTime, Parent_ID, Site_ID, Temperature, Dissolved_Oxygen_final) %>%
        rename(Dissolved_Oxygen = Dissolved_Oxygen_final)
      
    }
    
    
    
    ## ============================ remove biofouling =============================
    # biofouling identified with visual inspection
    
    # SSS004 - remove everything after 8/21
    # SSS005 - remove everything after 8/15
    # SSS016 - remove everything after 8/23
    # SSS028 - remove 8/22-8/23 and 8/26-8/27
    # SSS039 - remove 8/5-8/7
    # SSS046 - remove everything after 8/26
    
    if(parent_ID == 'SSS004'){
      
      clean_DO <- clean_DO %>%
        filter(DateTime<as_datetime('2022-08-22 00:00:00'))
      
    }else if(parent_ID == 'SSS005'){
      
      clean_DO <- clean_DO %>%
        filter(DateTime<as_datetime('2022-08-16 00:00:00'))
      
    }else if(parent_ID == 'SSS016'){
      
      clean_DO <- clean_DO %>%
        filter(DateTime<as_datetime('2022-08-24 00:00:00'))
      
    }else if(parent_ID == 'SSS028'){

      clean_DO <- clean_DO %>%
        filter(date(DateTime) != '2022-08-22') %>%
        filter(date(DateTime) != '2022-08-23') %>%
        filter(date(DateTime) != '2022-08-26') %>%
        filter(date(DateTime) != '2022-08-27')
      
    }else if(parent_ID == 'SSS039'){
      
      clean_DO <- clean_DO %>%
        filter(date(DateTime) != '2022-08-05') %>%
        filter(date(DateTime) != '2022-08-06') %>%
        filter(date(DateTime) != '2022-08-07') 
      
    }else if(parent_ID == 'SSS046'){
      
      clean_DO <- clean_DO %>%
        filter(DateTime<as_datetime('2022-08-27 00:00:00'))
      
    }
    
    ## ============================ remove data with dam influence =============================
    # dam influence identified with visual inspection

    
    if(parent_ID == 'SSS006'){
      
      
      
      clean_DO <- clean_DO %>%
        filter(date(DateTime) != '2022-08-07') %>%
        filter(date(DateTime) != '2022-08-08') %>%
        filter(date(DateTime) != '2022-08-14') %>%
        filter(date(DateTime) != '2022-08-15') %>%
        filter(date(DateTime) != '2022-08-21') %>%
        filter(date(DateTime) != '2022-08-22') %>%
        filter(date(DateTime) != '2022-08-23')%>%
        filter(DateTime < as_datetime('2022-08-27 00:00:00'))
      
    }
    ## =============== manually remove outliers/interpolate SSS001 ====================
      
    # some outliers not fully resolved with outlier detection, linearly interpolating 
    
      clean_DO <- clean_DO %>%
        mutate(Dissolved_Oxygen = case_when((DateTime >= as_datetime('2022-07-26 18:44:00') & DateTime <= as_datetime('2022-07-26 20:00:00')) ~ NA,
                                            TRUE ~ Dissolved_Oxygen),
               Dissolved_Oxygen = case_when((DateTime >= as_datetime('2022-08-02 05:14:00') & DateTime <= as_datetime('2022-08-02 05:45:00')) ~ NA,
                                            TRUE ~ Dissolved_Oxygen),
               Dissolved_Oxygen = round(zoo::na.approx(Dissolved_Oxygen, na.rm = FALSE), 3))
    
    ## =============== manually remove data ====================
    
    # weird DO sat on these days resulting in bad model fit, removing the days
    
    if(parent_ID == 'SSS008'){
      
      clean_DO <- clean_DO %>%
      filter(date(DateTime) != '2022-08-11') %>%
      filter(date(DateTime) != '2022-08-12')
      
    }else if(parent_ID == 'SSS017'){
      
      clean_DO <- clean_DO %>%
        filter(date(DateTime) != '2022-08-18') %>%
        filter(date(DateTime) != '2022-08-21') %>%
        filter(date(DateTime) != '2022-08-30')
      
    }

  ## ============================ create input file =============================
    
    ### ============================ merge DO with baro and hobo ===================
    
    baro_data <- baro_files[grepl(parent_ID, baro_files)] %>% # find baro file for same site
      read_csv(comment = '#', na = c('', '-9999',-9999, NA, 'N/A')) %>%
      select(-any_of("Air_Temperature")) %>%
      rename(BaroTROLL_Barometric_Pressure_mBar = Pressure)
    
    hobo_data <- hobo_files[grepl(parent_ID, hobo_files)] %>% # find hobo file for same site
      read_csv(comment = '#', na = c('', '-9999',-9999, NA, 'N/A')) %>%
      rename(HOBO_Temperature_degC = Temperature,
             HOBO_Absolute_Pressure_mbar = Absolute_Pressure)
    
    all_data <- hobo_data %>%
      full_join(baro_data, by = c('DateTime', 'Parent_ID', 'Site_ID')) %>%
      full_join(clean_DO, by = c('DateTime', 'Parent_ID', 'Site_ID'))
      
    ### ============================ calculate depth ===================
    
    site_avg_summary <- depth_summary %>%
      filter(Parent_ID == parent_ID)
    
    site_avg_depth_cm <- site_avg_summary %>%
      pull(Average_Depth)
    
    if(parent_ID == 'SSS024'){ 
      # (W20) hobo data file was lost, missing first half (before 8/13), 
      # using nearby site (W10) pressure to interpolate pressure and calculate depth 
      
      W10_hobo_data <- hobo_files[grepl('SSS036', hobo_files)] %>%
        read_csv(comment = '#', na = '-9999') %>%
        select(DateTime, Absolute_Pressure, Temperature) %>%
        rename(Pressure_W10 = Absolute_Pressure,
               Temp_W10 = Temperature)%>%
        filter(minute(DateTime) %% 15 == 0)
      
      join_W10_W20_full <- W10_hobo_data %>%
        full_join(all_data) %>%
        rename(Pressure_W20 = HOBO_Absolute_Pressure_mbar,
               Temp_W20 = HOBO_Temperature_degC)
      
      join_W10_W20_second_half <- join_W10_W20_full %>% # filter to time series where W20 has data to make the regression
        filter(!is.na(Pressure_W20))
      
      # linear model 
      pressure_lm <- lm(formula = join_W10_W20_second_half$Pressure_W20 ~ join_W10_W20_second_half$Pressure_W10)
      pressure_lm_summary <- summary(pressure_lm)
      
      temp_lm <- lm(formula = join_W10_W20_second_half$Temp_W20 ~ join_W10_W20_second_half$Temp_W10)
      temp_lm_summary <- summary(temp_lm)
      
      ggplot(join_W10_W20_full, aes(x = Pressure_W10, y = Pressure_W20))+
        geom_point()+
        geom_smooth(method = 'lm')
      
      ggplot(join_W10_W20_full, aes(x = Temp_W10, y = Temp_W20))+
        geom_point()+
        geom_smooth(method = 'lm')
      
      # Extract coefficients 
      pressure_intercept <- pressure_lm_summary$coefficients[1, 1]
      pressure_slope <- pressure_lm_summary$coefficients[2, 1]
      
      temp_intercept <- temp_lm_summary$coefficients[1, 1]
      temp_slope <- temp_lm_summary$coefficients[2, 1]
      
      # Perform the interpolation where W20 is missing data
      interp <- join_W10_W20_full %>%
        mutate(W20_Pressure_interp = case_when(
          is.na(Pressure_W20) & date(DateTime) < '2022-08-13' ~ (pressure_slope * Pressure_W10) + pressure_intercept,
          TRUE ~ Pressure_W20
        ),
        W20_Pressure_interp = round(W20_Pressure_interp, 1),
        W20_Temp_interp = case_when(
          is.na(Temp_W20) & date(DateTime) < '2022-08-13' ~ (temp_slope * Temp_W10) + temp_intercept,
          TRUE ~ Temp_W20
        ),
        W20_Temp_interp = round(W20_Temp_interp, 3)) %>% 
        select(DateTime, W20_Pressure_interp, W20_Temp_interp)
      
      all_data <- all_data %>%
        left_join(interp) %>%
        select(-HOBO_Absolute_Pressure_mbar, -HOBO_Temperature_degC) %>%
        rename(HOBO_Absolute_Pressure_mbar = W20_Pressure_interp,
               HOBO_Temperature_degC = W20_Temp_interp) # replace hobo data with interpolated data
    } 
    
    all_data_depth <- all_data %>%
      mutate(compensated_hobo_water_pressure_mbar = HOBO_Absolute_Pressure_mbar - BaroTROLL_Barometric_Pressure_mBar,
             density_kg_per_m3 = (999.84847 + (0.06337563 * HOBO_Temperature_degC) - (0.008523829 * HOBO_Temperature_degC^2) + (0.0000694324 * HOBO_Temperature_degC^3) - (0.0000003821216 * HOBO_Temperature_degC^4)),
             depth_from_pressure_m = (compensated_hobo_water_pressure_mbar*100)/(9.80  * density_kg_per_m3))
    
    if(parent_ID == 'SSS003'|parent_ID == 'SSS005'|parent_ID == 'SSS014'|parent_ID == 'SSS015'|parent_ID == 'SSS016'|parent_ID == 'SSS017'|parent_ID == 'SSS024'|parent_ID == 'SSS011'){ 
      
      # No data at time of transect depth, finding time closest,
      # these are the ones with the closest time AFTER the transect was taken
      
      hobo_data_reference_depth_m <- all_data_depth %>%
        filter(DateTime >= site_avg_summary$datetime) %>%
        head(1) %>%
        pull(depth_from_pressure_m)
      
      hobo_time_difference_combine <- hobo_time_difference_combine %>%
        add_row(Parent_ID = parent_ID,
        transect_depth_datetime =  site_avg_summary$datetime,
        hobo_reference_datetime = all_data_depth %>%
          filter(DateTime >= site_avg_summary$datetime) %>%
          head(1) %>%
          pull(DateTime))
      
    } else if(parent_ID == 'SSS010'|parent_ID == 'SSS023'|parent_ID == 'SSS036'){ 
      
      # No data at time of transect depth, finding time closest
      # these are the ones with the closest time BEFORE the transect was taken
      
      hobo_data_reference_depth_m <- all_data_depth %>%
        filter(DateTime <= site_avg_summary$datetime,
               !is.na(depth_from_pressure_m)) %>%
        tail(1) %>%
        pull(depth_from_pressure_m)
      
      hobo_time_difference_combine <- hobo_time_difference_combine %>%
        add_row(Parent_ID = parent_ID,
                transect_depth_datetime =  site_avg_summary$datetime,
                hobo_reference_datetime = all_data_depth %>%
                  filter(DateTime <= site_avg_summary$datetime) %>%
                  tail(1) %>%
                  pull(DateTime))
      
    }else{
      
      hobo_data_reference_depth_m <- all_data_depth %>%
        filter(DateTime <= site_avg_summary$datetime + 450 & DateTime >= site_avg_summary$datetime - 450,
               !is.na(depth_from_pressure_m)) %>%
        pull(depth_from_pressure_m)
      
    }
    
    if(parent_ID == 'SSS013'){
      # SSS013 (site T07) is missing hobo data for the beginning of the time series
      # resulting in the avg depth and the first hobo measurement being ~14 days apart. 
      # There is a very close USGS gage so we are using USGS discharge and a rating 
      # curve calculated from USGS historic data to estimate the depth 
      
      discharge <- './Stream_Metabolizer/Inputs/USGS_Kiona_Discharge.csv' %>%
        read_csv(comment = '#') %>%
        filter(agency_cd != '5s') %>%
        mutate(date = mdy(datetime),
               time = paste0(tz_cd, ":00"),       
               DateTime = as_datetime(paste(date, time)) - hours(1)) %>%
        rename(Discharge = '151886_00060_cd') %>%
        select(DateTime, Discharge)
      
      depth <- discharge %>%
        mutate(depth_ft = 0.0214 * (Discharge^0.6238),
               time_series_average_depth_cm = depth_ft*30.48) %>%
        select(DateTime, time_series_average_depth_cm)
      
      all_data_depth <- all_data_depth %>%
        left_join(depth)
      
      rm(discharge)
      rm(depth)
      
      
    } else {

    all_data_depth <- all_data_depth %>%
      mutate(offset_cm = (hobo_data_reference_depth_m * 100) - site_avg_depth_cm,
             time_series_average_depth_cm = (depth_from_pressure_m * 100)  - offset_cm) %>%
      filter(!is.na(Dissolved_Oxygen))
    
    }
    
    if(parent_ID == 'SSS024'){ # W10 was deployed two days later than W20 so filling in depth with avg of first day
      
      avg_depth_july29 <- all_data_depth %>%
        filter(date(DateTime) == '2022-07-29') %>%
        summarise(avg_depth = mean(time_series_average_depth_cm)) %>%
        pull()
      
      all_data_depth <- all_data_depth  %>%
        arrange(DateTime)%>%
        mutate(time_series_average_depth_cm = replace_na(time_series_average_depth_cm, avg_depth_july29))
    }

    
    ### ============================ output cleaned input file ===================
    
    input_file <- all_data_depth %>%
      rename(Depth = time_series_average_depth_cm,
             Pressure = BaroTROLL_Barometric_Pressure_mBar) %>%
      add_column(Latitude = metadata %>% filter(Parent_ID == parent_ID) %>% pull(Latitude),
                 Longitude = metadata %>% filter(Parent_ID == parent_ID) %>% pull(Longitude)) %>%
      mutate(Depth = round(Depth/100, 2)) %>%
      select(DateTime, Parent_ID, Site_ID, Latitude, Longitude, Temperature, Dissolved_Oxygen, Pressure, Depth) %>%
      arrange(DateTime)
    
    input_file_name <- str_c('./Stream_Metabolizer/Inputs/Sensor_Files/v2_', parent_ID, '_Temp_DO_Press_Depth.csv')
    
    # after the data are outputted, header rows with metadata are added
    write_csv(input_file, input_file_name, na = '-9999')
    
    ### ============================ plot all input data =============================
    
    input_file_plot1 <- ggplot(input_file, aes(x = DateTime, y = Dissolved_Oxygen)) +
      geom_point() + 
      labs(title = str_c("Parent ID: ", parent_ID, "                 Stream Metabolizer Input Data"),
           x = "DateTime", y = "Dissolved Oxygen (mg/L)", color = NULL) + 
      scale_x_datetime(labels = date_format("%Y-%m-%d %H:%M"))
    
    input_file_plot2 <- ggplot(input_file, aes(x = DateTime, y = Temperature)) +
      geom_point() + 
      labs(title = str_c("Parent ID: ", parent_ID, "                 Stream Metabolizer Input Data"),
           x = "DateTime", y = "Temperature (deg C)", color = NULL) + 
      scale_x_datetime(labels = date_format("%Y-%m-%d %H:%M"))
    
    input_file_plot3 <- ggplot(input_file, aes(x = DateTime, y = Pressure)) +
      geom_point() + 
      labs(title = str_c("Parent ID: ", parent_ID, "                 Stream Metabolizer Input Data"),
           x = "DateTime", y = "Barometric Pressure (mBar)", color = NULL) + 
      scale_x_datetime(labels = date_format("%Y-%m-%d %H:%M"))
    
    input_file_plot4 <- ggplot(input_file, aes(x = DateTime, y = Depth)) +
      geom_point() + 
      labs(title = str_c("Parent ID: ", parent_ID, "                 Stream Metabolizer Input Data"),
           x = "DateTime", y = "Depth (m)", color = NULL) + 
      scale_x_datetime(labels = date_format("%Y-%m-%d %H:%M"))
    
    
    plotly <- subplot(input_file_plot1, 
                      input_file_plot2, 
                      input_file_plot3,
                      input_file_plot4, 
                      nrows = 4, 
                      shareY = TRUE,
                      shareX = TRUE)
    
    #change wd to plot folder so it only outputs html and not additional files
    setwd('./Stream_Metabolizer/Inputs/Sensor_Files/Plots/')
    
    plotly_outname <-  str_c(parent_ID, '_Temp_DO_Press_Depth_Plot.html')

    htmlwidgets::saveWidget(as_widget(plotly), plotly_outname, selfcontained = T)
    
    #set back to original dir
    setwd(dirname(current_path))
    
}

hobo_time_difference_combine <- hobo_time_difference_combine %>%
  mutate(datetime_difference = transect_depth_datetime - hobo_reference_datetime)
