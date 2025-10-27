15 Nov 2024
Summary of the script from AI Incubator with edits from Brieanne Forbes


# Stream Metabolizer Input File Creation Script

## Script Description

1. **Setup and Initialization:**
   - Removes all existing objects in the environment.
   - Loads necessary libraries.
   - Sets the working directory to the location of the current script.

2. **File Discovery:**
   - Identifies relevant data files (Dissolved Oxygen, metadata, barometric pressure, depth) from the SSS data package directories.

3. **Plot Theme Setup:**
   - Defines a common theme for all generated plots.

4. **Empty Tibble Creation:**
   - Creates an empty tibble to store time differences between transect depths and sensor measurement times.

5. **Data Processing Loop:**
   - Iterates through Dissolved Oxygen files, performing the following steps for each file:
     - Reads the data file and extracts relevant columns.
     - Subsets the data to 15-minute intervals.
     - Cleans the data using the `tsrobprep` library to remove outliers and fix missing values.
     - Adjusts data for specific sample days based on metadata.
     - Ensures cleaned values are within expected ranges, reverting to original data if necessary.

6. **Biofouling Removal:**
   - Removes data points affected by biofouling based on pre-identified dates for specific sites:
      - **`SSS004`**: Remove everything after 8/21
      - **`SSS005`**: Remove everything after 8/15
      - **`SSS016`**: Remove everything after 8/23
      - **`SSS028`**: Remove 8/22-8/23 and 8/26-8/27
      - **`SSS039`**: Remove 8/5-8/7
      - **`SSS046`**: Remove everything after 8/26

7. **Dam Influence Removal:**
   - Removes data points affected by dam influence for Parent ID `SSS006`:
      - Remove data on 8/7, 8/8, 8/14, 8/15, 8/21, 8/22, and 8/23
      - Remove everything after 8/26

8. **Manual Outlier Removal and Interpolation for `SSS001`:**
   - Manually removes and linearly interpolates outliers on 7/26 and 8/2.

9. **Special Handling for `SSS008` and `SSS017`**
   - Removes data where there is poor model fit
   - SSS008: Removes data on 8/11 and 8/12 
   - SSS017: Removes data on 8/18, 8/21, and 8/30

10. **Input File Creation:**
   - Merges Dissolved Oxygen data with barometric and hobo (depth) data.
   - Calculates water depth from pressure using site-specific adjustments and linear models for missing data.
   - Handles missing depth data for specific sites as follows:
     - **`SSS024` (W20)**: Interpolates pressure and temperature from nearby site `SSS036` (W10) due to missing hobo data.
     - **`SSS013` (T07)**: Uses USGS discharge and a rating curve to estimate depth due to missing hobo data for the beginning of the time series.
   - Computes average depth and time series depth with site-specific offsets for other Parent IDs with special handling of calculation timing:
     - **`SSS003`, `SSS005`, `SSS014`, `SSS015`, `SSS016`, `SSS017`, `SSS024`, `SSS011`**: Finds the closest time AFTER the transect depth was taken.
     - **`SSS010`, `SSS023`, `SSS036`**: Finds the closest time BEFORE the transect depth was taken.

11. **Output Files:**
   - Creates cleaned input files with necessary columns for the Stream Metabolizer tool.
   - Saves the input files to `./Stream_Metabolizer/Inputs/Sensor_Files/`.

12. **Plot Generation:**
   - Generates interactive plots for Dissolved Oxygen, Temperature, Pressure, and Depth.
   - Saves the plots as HTML files in `./Stream_Metabolizer/Inputs/Sensor_Files/Plots/`.

13. **Time Difference Calculation:**
   - Computes the time differences between transect depth measurements and hobo data for specific sites.

## Output

- **Cleaned Input Files:** Located in `./Stream_Metabolizer/Inputs/Sensor_Files/`.
- **Interactive Plots:** Saved in `./Stream_Metabolizer/Inputs/Sensor_Files/Plots/`.

## Usage

Run the script in an R environment with the above libraries installed. Ensure the working directory is set to the location of this script, and the required data files are located in the specified paths within the SSS data package directory.