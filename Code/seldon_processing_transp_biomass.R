# Title     : Harvest data analysis
# Created by: Terence Seldon Kwafo
# Created on: 11/09/2022 .... Modified 07/25/2024


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

#load libraries 
list.of.packages <- c("dplyr","here","plotrix",
                      "gridExtra", "lattice", "DataCombine", "latticeExtra",
                      "reshape2", "squash", "Hmisc", "signal", "rgl", "corrplot",
                      "lme4", "lmerTest", "tidyr", "pracma", "TREX", "ggplot2", "ggpubr", 
                      "RcppRoll","PhotoGEA", "BioCroField", "BioCro")
invisible(lapply(list.of.packages, library, character.only = TRUE))

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
# Notes:
# biomass
# measurement               unit    
# length_row_biomass_part - metres (m)
# length_total_abvg_biomass - metres (m)
# row_spacing             - metres (m)
# grain_yield             - grams (g) - 2m row
# leaf_area               - square centimetres (cm2)
# leaf_mass               - grams (g)
# stem_mass               - grams (g)
# root_mass               - grams (g)
# pod_mass                - grams (g)
# total_abvg_biomass      - grams (g)
# leaf_litter_mass         - grams (g)


# nutrients
# Element - Unit
# N	        %	
# S         %	
# P         %	
# K         %	
# MG        %	
# Ca        %	
# Na        %
# B	        ppm = mg/kg
# Zn	      ppm = mg/kg
# Mn	      ppm = mg/kg
# Fe	      ppm = mg/kg
# Cu	      ppm = mg/kg
# Al        ppm uom


#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

################################################################
############  Data import  ####################################

getwd()
setwd("C:/Users/seldon/Box/Transpiration&Nutrient/biomass&nutrient/")
setwd("C:/Users/seldon/Box/Transpiration&Nutrient/biomass&nutrient/scipts/")
source('settings.R')




# Set the main directory
main_dir <- "C:/Users/seldon/Box/Transpiration&Nutrient/biomass&nutrient/collected data"
output_dir <-"C:/Users/seldon/Box/Transpiration&Nutrient/biomass_output"

# List all subdirectories (one for each year)
years <- list.dirs(main_dir, recursive = FALSE)

# Loop through each year and assign data frames
for (year_dir in years) {
  # Extract the year from the directory name
  year <- basename(year_dir)
  
  # Construct file paths for biomass and nutrient data

  nutrient_file <- file.path(year_dir, paste0(year, "_soybean_nutrient_data", ".csv"))
  biomass_file <- file.path(year_dir, paste0(year, "_soybean_biomass_data", ".csv"))
  
  assign(paste0("nutrient_data_", year), read.csv(nutrient_file))
  assign(paste0("biomass_data_", year), read.csv(biomass_file))
}

# Check one of the assigned data frames
print(nutrient_data_2022)


######################################################################
############  Conversions  ################################

## Convert nutrient percentages to ppm for columns N to Na and save to a new dataframe

# Function to convert nutrient percentages to ppm
# Nutrient N to N are given in percent by Waypoint
# To convert to ppm or mg/kg multiply by 10,0000 ppm

convert_to_ppm <- function(data, columns) {
  data[columns] <- lapply(data[columns], function(x) (x / 100) * 10000)
  return(data)
}

# List of nutrient columns to convert
nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Na")

# Convert nutrient percentages to ppm for each year's data
nutrient_data_list <- list(nutrient_data_2022, nutrient_data_2023) # Add nutrient_data_2024 if applicable
converted_nutrient_data_list <- lapply(nutrient_data_list, convert_to_ppm, columns = nutrient_columns)

# Extract converted data for each year
nutrient_data_2022_ppm <- converted_nutrient_data_list[[1]]
nutrient_data_2023_ppm <- converted_nutrient_data_list[[2]]
# nutrient_data_2024_ppm <- converted_nutrient_data_list[[3]] # Uncomment if applicable

# Display the first few rows of the converted data
cat("Nutrient Data in ppm for 2022:\n")
print(head(nutrient_data_2022_ppm))

cat("Nutrient Data in ppm for 2023:\n")
print(head(nutrient_data_2023_ppm))

# cat("Nutrient Data in ppm for 2024:\n")
# print(head(nutrient_data_2024_ppm)) 


######################################################################
### Biomass data

#Note: All units are in grams for data

# # Ensure the necessary columns exist in both datasets
# required_columns <- c("treatment", "variety", "tissue_type", "development_stage", "year")
# if (!all(required_columns %in% colnames(biomass_data_2022)) || !all(required_columns %in% colnames(biomass_data_2023))) {
#   stop("Both datasets must contain 'treatment', 'variety', 'tissue_type', 'development_stage', and 'year' columns.")
# }


# Define biomass columns
biomass_columns <- c("leaf_mass", "stem_mass", "root_mass", "pod_mass","length_total_biomass", 
                     "total_abvg_biomass", "grain_yield", "litter_mass", "leaf_area")

# Combine the raw data
combined_biomass_data <- bind_rows(
  biomass_data_2022 %>% mutate(year = 2022),
  biomass_data_2023 %>% mutate(year = 2023),
  biomass_data_2024 %>% mutate(year = 2024)
)

#Note: Here we want to find the tissue yield by calculating the 
# 1 - Above ground area (AGB) (g/m2) = Total AGB (1 m/2 m)/ row spacing (m) then
# 2 - Tissue above ground ratio (dimensionless) = Tissue mass (Partitioning) (g) /AGB (Partitioning = leaf+stem+pods) (g)
# Therefore Tissue yield (g/m2)= Tissue AGB ratio * Above-Ground Biomass (AGB) (g/m2)
# Convertions
# - 1 g / m^2 * (1 Mg / 1e6 g) * (1e4 m^2 / 1 ha) = 1e-2 Mg / ha
# - 1 m^2 / g * (1e6 g / 1 Mg) * (1 ha / 1e4 m^2) = 1e2 ha / Mg
# - 1 g / cm^2 * (100 cm / m)^2 = 1e4 g / m^2


# Define the row spacing 
row_spacing <- 0.762  # in meters

# Calculate AGB (Above-Ground Biomass) for each observation considering variable row length
combined_biomass_df <- combined_biomass_data

combined_biomass_df <- combined_biomass_df %>%
  mutate(AGB_per_m2 = total_abvg_biomass * (1 / length_total_biomass) * (1 / row_spacing))

# Finding tissue yield
combined_biomass_df <- combined_biomass_df %>%
  # Calculate the Tissue AGB ratio using the sum of leaf_mass, stem_mass, and pod_mass
  mutate(
    leaf_AGB_ratio = leaf_mass / (coalesce(leaf_mass, 0) + coalesce(stem_mass, 0) + coalesce(pod_mass, 0)),
    stem_AGB_ratio = stem_mass / (coalesce(leaf_mass, 0) + coalesce(stem_mass, 0) + coalesce(pod_mass, 0)),
    root_AGB_ratio = root_mass / (coalesce(leaf_mass, 0) + coalesce(stem_mass, 0) + coalesce(pod_mass, 0)),
    pod_AGB_ratio = pod_mass / (coalesce(leaf_mass, 0) + coalesce(stem_mass, 0) + coalesce(pod_mass, 0))
  )%>%
  mutate(
    leaves_tissue_yield = leaf_AGB_ratio * AGB_per_m2,
    stems_tissue_yield = stem_AGB_ratio * AGB_per_m2,
    roots_tissue_yield = root_AGB_ratio * AGB_per_m2,
    pods_tissue_yield  = pod_AGB_ratio * AGB_per_m2,
    grain_tissue_yield = grain_yield * (1 / length_total_biomass) * (1 / row_spacing)
  )


# Convert the wide format data to long format
combined_biomass_long <- combined_biomass_df%>%
  pivot_longer(
    cols = starts_with("leaves_tissue_Yield"):starts_with("grain_tissue_yield"),  # Select all Tissue_Yield columns
    names_to = "tissue_type",  # New column name for tissue types
    values_to = "tissue_yield"  # New column name for tissue yield values
  ) %>%
  # Clean up the tissue_type column to remove "_Tissue_Yield" and just keep "leaves", "stems", "roots", "pods"
  mutate(tissue_type = sub("_tissue_yield", "", tissue_type))


#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

# Define biomass columns
biomass_columns_2 <- c("leaf_mass", "stem_mass", "root_mass", "pod_mass","leaf_area", 
                     "total_abvg_biomass", "grain_yield", "litter_mass", 
                     "leaves_tissue_yield", "stems_tissue_yield", "roots_tissue_yield",
                     "pods_tissue_yield", "grain_tissue_yield")

# Treatment means for biomass columns 
mean_biomass_ci <- combined_biomass_df %>%
  group_by(treatment, variety, development_stage, year) %>%
  summarise(across(all_of(biomass_columns_2), 
                   list(
                     mean = ~mean(., na.rm = TRUE), 
                     sd = ~sd(., na.rm = TRUE), 
                     n = ~sum(!is.na(.)),
                     se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),  
                     ci_lower = ~mean(., na.rm = TRUE) - 1.96 * (sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))),
                     ci_upper = ~mean(., na.rm = TRUE) + 1.96 * (sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))))
                   ), 
                   .names = "{col}_{fn}"))

#########################################################################################################

# Chart function to plot graphs for means of biomass components

biomass_barchart <- function(
    big_df,
    development_stage,
    biomass_column,
    y_axis_title = "Mean Value", 
    output_folder = NULL
) {
  # Filter data
  rows_to_keep <- big_df$development_stage %in% development_stage
  big_df_subset <- big_df[rows_to_keep, ]
  
  # construct column names
  mean_col <- biomass_column
  ci_lower_col <- gsub("_mean", "_ci_lower", biomass_column)
  ci_upper_col <- gsub("_mean", "_ci_upper", biomass_column)
  
  # Check that the columns exist
  if (!all(c(mean_col, ci_lower_col, ci_upper_col) %in% colnames(big_df_subset))) {
    stop(paste("Missing columns for", biomass_column))
  }
  
  # Prepare plot data
  plot_data <- big_df_subset %>%
    select(treatment, variety, year, 
           !!sym(mean_col), 
           !!sym(ci_lower_col), 
           !!sym(ci_upper_col)) %>%
    rename(
      mean = !!sym(mean_col),
      ci_lower = !!sym(ci_lower_col),
      ci_upper = !!sym(ci_upper_col)
    )
  
  # Update variety labels to uppercase
  plot_data$variety <- factor(plot_data$variety, 
                              levels = c("hs", "loda"), 
                              labels = c("HS", "Loda"))
  
  # Create the bar chart
  p <- ggplot(plot_data, aes(x = treatment, y = mean, fill = variety)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), 
             width = 0.7, color = "black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, position = position_dodge(0.8), 
                  color = "black", size = 1.2) + 
    scale_fill_manual(values = c("#009999", "#FF6633"), name = "Variety") +
    labs(
      x = "Treatment",
      y = y_axis_title, 
      fill = "Variety"
    ) +
    theme_classic2() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      strip.background = element_blank(),
      axis.text = element_text(size = 20, colour = "black"),
      axis.title = element_text(size = 24),
      axis.title.x = element_text(vjust = -2),
      axis.title.y = element_text(vjust = 3),
      plot.title = element_text(size = 36),
      plot.subtitle = element_text(size = 24),
      strip.text = element_text(size = 24),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 24),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.spacing = unit(1.5, "lines")
    ) +
    facet_wrap(~year) +
    scale_x_discrete(labels = c("Ambient", "Elevated"))

  print(p)
}


library(stringr) #For str_to_title
# Define the biomass columns of interest
biomass_columns <- c("leaf_area_mean", "root_mass_mean", "grain_yield_mean", 
                     "total_abvg_biomass_mean", "leaves_tissue_yield_mean", "roots_tissue_yield_mean")

# Define parameters
year <- c(2022, 2023, 2024)  # Specify the years of interest


# Call the biomass_barchart function

# Leaf area (cm2) - check R5 and R6
development_stage <- c("R6")
biomass_barchart(
  big_df = mean_biomass_ci,
  development_stage,
  biomass_column = "leaf_area_mean",
  y_axis_title = "Leaf Area (cm2) "
)

# Root mass (g)
development_stage <- c("R8")
biomass_barchart(
  big_df = mean_biomass_ci,
  development_stage,
  biomass_column = "root_mass_mean",
  y_axis_title = "Root Mass (g) "
)

# Grain yield (g/m2)
development_stage <- c("R8")
biomass_barchart(
  big_df = mean_biomass_ci,
  development_stage,
  biomass_column = "grain_yield_mean",
  y_axis_title = "Grain Yield (g/m2) "
)

# Above ground biomass (g)
development_stage <- c("R8")
biomass_barchart(
  big_df = mean_biomass_ci,
  development_stage,
  biomass_column = "total_abvg_biomass_mean",
  y_axis_title = "Total Aboveground Biomass (g)"
)

# Leaf yield (g/m2)
development_stage <- c("R6")
biomass_barchart(
  big_df = mean_biomass_ci,
  development_stage,
  biomass_column = "leaves_tissue_yield_mean",
  y_axis_title = "Leaf Yield (g/m2) "
)

# Root yield (g/m2)
development_stage <- c("R8")
biomass_barchart(
  big_df = mean_biomass_ci,
  development_stage,
  biomass_column = "roots_tissue_yield_mean",
  y_axis_title = "Root Yield (g/m2)"
)

#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
######## Creating biomass relative change graph
# Define a function to calculate relative change and plot for biomass data
biomass_barchart<- function(
    big_df,
    development_stage,
    year,
    caption,
    biomass_columns
)
{
  # Filter the data based on development stage, year
  rows_to_keep <- 
    big_df$development_stage %in% development_stage & 
    big_df$year %in% year
  
  big_df_subset <- big_df[rows_to_keep, ]
  
  # Calculate relative change and confidence intervals for biomass components
  diff_list <- by(
    big_df_subset,
    list(big_df_subset$development_stage, big_df_subset$year, big_df_subset$variety),
    function(x) {
      ambient <- x[x$treatment == 'ambient', , drop = FALSE]  
      elevated <- x[x$treatment == 'elevated', , drop = FALSE]  
      
      if (nrow(ambient) == 0 | nrow(elevated) == 0) return(NULL)  # Ensure both are non-empty data frames
      
      relative_change <- lapply(biomass_columns_2, function(biomass) {
        biomass_mean_ambient <- paste0(biomass, "_mean")
        biomass_mean_elevated <- paste0(biomass, "_mean")
        biomass_sd_ambient <- paste0(biomass, "_sd")
        biomass_sd_elevated <- paste0(biomass, "_sd")
        biomass_n_ambient <- paste0(biomass, "_n")
        biomass_n_elevated <- paste0(biomass, "_n")
        
        ambient_mean <- ambient[[biomass_mean_ambient]]
        elevated_mean <- elevated[[biomass_mean_elevated]]
        ambient_sd <- ambient[[biomass_sd_ambient]]
        elevated_sd <- elevated[[biomass_sd_elevated]]
        n_ambient <- ambient[[biomass_n_ambient]]
        n_elevated <- elevated[[biomass_n_elevated]]
        
        # Ensure there are non-NA values
        if (is.na(ambient_mean) | is.na(elevated_mean) | is.null(ambient_mean) | is.null(elevated_mean)) {
          return(list(relative_change = NA, ci_lower = NA, ci_upper = NA))
        }
        
        # Calculate standard errors
        sem_ambient <- ambient_sd / sqrt(n_ambient)
        sem_elevated <- elevated_sd / sqrt(n_elevated)
        
        # Calculate the relative change and the CI
        relative_change <- (elevated_mean - ambient_mean) / ambient_mean * 100
        
        # 95% confidence intervals (using 1.96 as a t-score )
        ci_lower <- relative_change - 1.96 * sem_elevated
        ci_upper <- relative_change + 1.96 * sem_elevated
        
        list(
          relative_change = relative_change,
          ci_lower = ci_lower,
          ci_upper = ci_upper
        )
      })
      
      # Combine results into a data frame
      relative_changes <- sapply(relative_change, function(rc) rc$relative_change)
      ci_lowers <- sapply(relative_change, function(rc) rc$ci_lower)
      ci_uppers <- sapply(relative_change, function(rc) rc$ci_upper)
      
      data.frame(
        biomass = biomass_columns_2,
        relative_change = relative_changes,
        ci_lower = ci_lowers,
        ci_upper = ci_uppers,
        development_stage = x$development_stage[1],
        year = x$year[1],
        variety = x$variety[1]
      )
    }
  )
  
  diff_df <- do.call(rbind, diff_list)
  
  if (is.null(diff_df)) {
    stop("No valid data for the selected combination.")
  }
  
  diff_df$biomass <- factor(diff_df$biomass, levels = biomass_columns_2)
  
  # Plot

  p <- ggplot(diff_df, aes(x = biomass, y = relative_change, fill = variety)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), 
             width = 0.7, 
             color="black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, 
                  position = position_dodge(0.8), 
                  color = "black", 
                  size = 1.2) + 
    scale_fill_brewer(palette = "Set3") +  
    #facet_grid(~ variety + development_stage) +  # Facet by variety and development stage
    ggtitle(paste("Mean biomass",development_stage, year)) +
    xlab("Biomass components") +
    ylab("Relative change in Biomass (%)") +
    theme_classic() +
    theme(        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Change axis text size
                  axis.text.y = element_text(size = 14),  # Change y-axis text size
                  axis.title.x = element_text(size = 16),  # Change x-axis label size
                  axis.title.y = element_text(size = 16),  # Change y-axis label size
                  plot.title = element_text(size = 18, face = "bold"),  # Change plot title size
                  legend.text = element_text(size = 14),  # Change legend text size
                  legend.title = element_text(size = 16)  # Change legend title size
    )

  
  print(p)
}


#############################################################################################

biomass_barchart <- function(
    big_df,
    development_stage,
    year,
    caption,
    output_folder = NULL
)
{
  # Filter the data based on development stage and year
  rows_to_keep <- 
    big_df$development_stage %in% development_stage & 
    big_df$year %in% year
  
  big_df_subset <- big_df[rows_to_keep, ]
  
  # Calculate relative change and confidence intervals for biomass components
  diff_list <- by(
    big_df_subset,
    list(big_df_subset$development_stage, big_df_subset$year, big_df_subset$variety),
    function(x) {
      ambient <- x[x$treatment == 'ambient', , drop = FALSE]  
      elevated <- x[x$treatment == 'elevated', , drop = FALSE]  
      
      if (nrow(ambient) == 0 | nrow(elevated) == 0) return(NULL)  # Ensure both are non-empty data frames
      
      relative_change <- lapply(biomass_columns_2, function(biomass) {
        biomass_mean_ambient <- paste0(biomass, "_mean")
        biomass_mean_elevated <- paste0(biomass, "_mean")
        biomass_sd_ambient <- paste0(biomass, "_sd")
        biomass_sd_elevated <- paste0(biomass, "_sd")
        biomass_n_ambient <- paste0(biomass, "_n")
        biomass_n_elevated <- paste0(biomass, "_n")
        
        ambient_mean <- ambient[[biomass_mean_ambient]]
        elevated_mean <- elevated[[biomass_mean_elevated]]
        ambient_sd <- ambient[[biomass_sd_ambient]]
        elevated_sd <- elevated[[biomass_sd_elevated]]
        n_ambient <- ambient[[biomass_n_ambient]]
        n_elevated <- elevated[[biomass_n_elevated]]
        
        # Ensure there are non-NA values
        if (is.na(ambient_mean) | is.na(elevated_mean) | is.null(ambient_mean) | is.null(elevated_mean)) {
          return(list(relative_change = NA, ci_lower = NA, ci_upper = NA))
        }
        
        # Calculate standard errors
        sem_ambient <- ambient_sd / sqrt(n_ambient)
        sem_elevated <- elevated_sd / sqrt(n_elevated)
        
        # Calculate the relative change
        relative_change <- (elevated_mean - ambient_mean) / ambient_mean * 100
        
        # Error propagation for 95% confidence intervals
        ratio_se <- relative_change * sqrt((sem_elevated / elevated_mean)^2 + (sem_ambient / ambient_mean)^2)
        ci_lower <- relative_change - 1.96 * ratio_se
        ci_upper <- relative_change + 1.96 * ratio_se
        
        list(
          relative_change = relative_change,
          ci_lower = ci_lower,
          ci_upper = ci_upper
        )
      })
      
      # Combine results into a data frame
      relative_changes <- sapply(relative_change, function(rc) rc$relative_change)
      ci_lowers <- sapply(relative_change, function(rc) rc$ci_lower)
      ci_uppers <- sapply(relative_change, function(rc) rc$ci_upper)
      
      data.frame(
        biomass = biomass_columns_2,
        relative_change = relative_changes,
        ci_lower = ci_lowers,
        ci_upper = ci_uppers,
        development_stage = x$development_stage[1],
        year = x$year[1],
        variety = x$variety[1]
      )
    }
  )
  
<<<<<<< Updated upstream
  diff_df <- do.call(rbind, diff_list)
=======
  do.call(rbind, diff_list)
}

# Step 1: Filter the dataset for specific development stages and years
development_stage <- c("R8")  # Specify the development stages of interest
year <- c(2022, 2023, 2024)  # Specify the years of interest

# Filter the data using the `filter_data` function
big_df_subset <- filter_data(mean_biomass_ci, development_stage, year)

# Step 2: Define the biomass columns to analyze
biomass_columns_2 <- c("grain_yield")

# Step 3: Call the `calculate_relative_change` function
relative_changes_df <- calculate_relative_change(big_df_subset, biomass_columns_2)

# View the results
print(relative_changes_df)


####### Bar chart plot
biomass_barchart <- function(
    big_df,
    development_stage,
    year,
    biomass_columns,
    output_folder = NULL
) {
  # Filter data
  rows_to_keep <- 
    big_df$development_stage %in% development_stage & 
    big_df$year %in% year
  
  big_df_subset <- big_df[rows_to_keep, ]
  
  # Calculate relative change
  diff_df <- calculate_relative_change(big_df_subset, biomass_columns)
>>>>>>> Stashed changes
  
  if (is.null(diff_df)) {
    stop("No valid data for the selected combination.")
  }
  
  diff_df$biomass <- factor(diff_df$biomass, levels = biomass_columns_2)
  
  # Plot
  p <- ggplot(diff_df, aes(x = biomass, y = relative_change, fill = variety)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), 
             width = 0.7, 
             color = "black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, 
                  position = position_dodge(0.8), 
                  color = "black", 
                  size = 1.2) + 
    scale_fill_brewer(palette = "Set3") +  
    ggtitle(paste("Relative Change in Biomass", development_stage, year)) +
    xlab("Biomass Components") +
    ylab("Relative Change (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Change axis text size
      axis.text.y = element_text(size = 14),                         # Change y-axis text size
      axis.title.x = element_text(size = 16),                        # Change x-axis label size
      axis.title.y = element_text(size = 16),                        # Change y-axis label size
      plot.title = element_text(size = 18, face = "bold"),           # Change plot title size
      legend.text = element_text(size = 14),                         # Change legend text size
      legend.title = element_text(size = 16)                         # Change legend title size
    )
  
  print(p)
}

<<<<<<< Updated upstream
biomass_barchart(mean_biomass_ci, 'V7_R1', '2022', 'R1 Biomass 2022', biomass_columns_2)
biomass_barchart(mean_biomass_ci, 'R5', '2022','R5 Biomass 2022', biomass_columns_2)
biomass_barchart(mean_biomass_ci, 'R6', '2022', 'R6 Biomass 2022', biomass_columns_2)
biomass_barchart(mean_biomass_ci, 'R8', '2022', 'R8 Biomass 2022', biomass_columns_2)
=======
# Define your parameters
development_stage <- c("R8")  # Specify the development stages of interest
year <- c(2022, 2023, 2024)  # Specify the years of interest
biomass_columns_2 <- c("grain_yield")  # Example biomass components
>>>>>>> Stashed changes

### 2023
biomass_barchart(mean_biomass_ci, 'V7_R1', '2023', 'R1 Biomass 2023', biomass_columns_2)
biomass_barchart(mean_biomass_ci, 'R5', '2023','R5 Biomass 2023', biomass_columns_2)
biomass_barchart(mean_biomass_ci, 'R6', '2023', 'R6 Biomass 2022', biomass_columns_2)
biomass_barchart(mean_biomass_ci, 'R8', '2023', 'R8 Biomass 2023', biomass_columns_2)

###################################################################################################
### Nutrient data 

#Note: All units are in ppm or mg/kg for raw data

# # Ensure the necessary columns exist in both datasets
# required_columns <- c("treatment", "variety", "tissue_type", "development_stage", "year")
# if (!all(required_columns %in% colnames(nutrient_data_2022_ppm)) || !all(required_columns %in% colnames(nutrient_data_2023_ppm))) {
#   stop("Both datasets must contain 'treatment', 'variety', 'tissue_type', 'development_stage', and 'year' columns.")
# }

# Define nutrient columns
nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")

# Combine the raw data
combined_nutrient_data <- bind_rows(
  nutrient_data_2022_ppm %>% mutate(year = 2022),
  nutrient_data_2023_ppm %>% mutate(year = 2023)
)

# Data frame to use for nutrient yield calculations 
combined_nutrient_df <- combined_nutrient_data

# Convert nutrient concentrations from ppm (mg/kg) to g/g
#1 ppm = 1 mg/kg = 1×10−6 g/g

nutrient_columns <- colnames(combined_nutrient_df)[which(colnames(combined_nutrient_df) == "N"):which(colnames(combined_nutrient_df) == "Al")]
combined_nutrient_df[nutrient_columns] <- combined_nutrient_df[nutrient_columns] * 1e-6

# Merge nutrient data with biomass data on common columns
combined_data_df <- merge(combined_nutrient_df , combined_biomass_long, 
                          by = c("treatment", "variety", "block_id", "development_stage", "year", "ring_id", "tissue_type"))


#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#



#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

### Find response ratio using nutrient concentrations
# Note: The response ratio, r = E/A, where A(ambient co2) and E(elevated co2) are the mean concentrations of an element
# (E - A) / A) * 100
# 1 - Stems and Leaves - R1 & R5
# 2 - Pods - R8
# Test for 95% confidence limits

# Calculate Mean and Standard Error for Nutrient Concentration
nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")

mean_nutrient_concentration_ci <- combined_data_df %>%
  group_by(treatment, variety, tissue_type, development_stage, year) %>%
  summarise(across(all_of(nutrient_columns), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        n = ~sum(!is.na(.))), 
                   .names = "{col}_{fn}"))

<<<<<<< Updated upstream
# Define a function to calculate confidence intervals
element_barchart_df_ggplot <- function(
    big_df,
    development_stage,
    year,
    tissue,
    caption
)
{
  # Filter the data based on development stage, year, and tissue type
  rows_to_keep <- 
    big_df$development_stage %in% development_stage & 
    big_df$year %in% year & 
    big_df$tissue_type %in% tissue
  
  big_df_subset <- big_df[rows_to_keep, ]
  
  element_names <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")
  
  # Calculate relative change and confidence intervals
=======
# Write the treatment_cumulative_flow data frame to a CSV file
write.csv(mean_nutrient_concentration_ci, "mean_nutrient_concentration.csv", row.names = FALSE)

# Display confirmation message
cat("Data successfully written to mean_nutrient_concentration.csv\n")

calculate_relative_change_nutrient <- function(big_df_subset, element_names) {
>>>>>>> Stashed changes
  diff_list <- by(
    big_df_subset,
    list(big_df_subset$development_stage, big_df_subset$year, big_df_subset$tissue_type, big_df_subset$variety),
    function(x) {
      ambient <- x[x$treatment == 'ambient', , drop = FALSE] 
      elevated <- x[x$treatment == 'elevated', , drop = FALSE]
      
      if (nrow(ambient) == 0 | nrow(elevated) == 0) return(NULL)  # Ensure both are non-empty data frames
      
      relative_change <- lapply(element_names, function(en) {
        en_mean_ambient <- paste0(en, "_mean")
        en_mean_elevated <- paste0(en, "_mean")
        en_sd_ambient <- paste0(en, "_sd")
        en_sd_elevated <- paste0(en, "_sd")
        en_n_ambient <- paste0(en, "_n")
        en_n_elevated <- paste0(en, "_n")
        
        ambient_mean <- ambient[[en_mean_ambient]]
        elevated_mean <- elevated[[en_mean_elevated]]
        ambient_sd <- ambient[[en_sd_ambient]]
        elevated_sd <- elevated[[en_sd_elevated]]
        n_ambient <- ambient[[en_n_ambient]]
        n_elevated <- elevated[[en_n_elevated]]
        
        # Ensure there are non-NA values
        if (is.na(ambient_mean) | is.na(elevated_mean) | is.null(ambient_mean) | is.null(elevated_mean)) {
          return(list(relative_change = NA, ci_lower = NA, ci_upper = NA))
        }
        
        # Calculate standard errors
        sem_ambient <- ambient_sd / sqrt(n_ambient)
        sem_elevated <- elevated_sd / sqrt(n_elevated)
        
        # Calculate the relative change and error propagation for CI
        relative_change <- (elevated_mean - ambient_mean) / ambient_mean * 100
        
        # Error propagation for 95% confidence intervals
        ratio_se <- relative_change * sqrt((sem_elevated / elevated_mean)^2 + (sem_ambient / ambient_mean)^2)
        ci_lower <- relative_change - 1.96 * ratio_se
        ci_upper <- relative_change + 1.96 * ratio_se
        
        list(
          relative_change = relative_change,
          ci_lower = ci_lower,
          ci_upper = ci_upper
        )
      })
      
      # Combine results into a data frame
      relative_changes <- sapply(relative_change, function(rc) rc$relative_change)
      ci_lowers <- sapply(relative_change, function(rc) rc$ci_lower)
      ci_uppers <- sapply(relative_change, function(rc) rc$ci_upper)
      
      data.frame(
        nutrient = element_names,
        relative_change = relative_changes,
        ci_lower = ci_lowers,
        ci_upper = ci_uppers,
        development_stage = x$development_stage[1],
        year = x$year[1],
        tissue_type = x$tissue_type[1],
        variety = x$variety[1]
      )
    }
  )
  
  diff_df <- do.call(rbind, diff_list)
  
  if (is.null(diff_df)) {
    stop("No valid data for the selected combination.")
  }
  
  diff_df$nutrient <- factor(diff_df$nutrient, levels = element_names)
  
  # Plot
  p <- ggplot(diff_df, aes(x = nutrient, y = relative_change, fill = interaction(variety, tissue_type))) +
    geom_bar(stat = "identity", position = position_dodge(0.8), 
             width = 0.7, 
             color="black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, 
                  position = position_dodge(0.8), 
                  color = "black", 
                  size = 1.2) + 
    scale_fill_brewer(palette = "Set3", name = "Tissue") +  
    ggtitle(caption) +
    xlab("Nutrient") +
    ylab("Relative change in concentration (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Change axis text size
      axis.text.y = element_text(size = 14),                         # Change y-axis text size
      axis.title.x = element_text(size = 16),                        # Change x-axis label size
      axis.title.y = element_text(size = 16),                        # Change y-axis label size
      plot.title = element_text(size = 18, face = "bold"),           # Change plot title size
      legend.text = element_text(size = 14),                         # Change legend text size
      legend.title = element_text(size = 16)                         # Change legend title size
    )
  
  print(p)
}


<<<<<<< Updated upstream
### 2022
element_barchart_df_ggplot(mean_nutrient_concentration_ci, 'V7_R1', '2022', c('leaves', 'stems'), 'R1 leaves & stems')
element_barchart_df_ggplot(mean_nutrient_concentration_ci, 'R5', '2022', c('leaves', 'stems'), 'R5 leaves & stems')
element_barchart_df_ggplot(mean_nutrient_concentration_ci, 'R8', '2022', 'pods', 'R8 pods')

### 2023
element_barchart_df_ggplot(mean_nutrient_concentration_ci, 'V7_R1', '2023', c('leaves', 'stems'), 'R1 leaves & stems')
element_barchart_df_ggplot(mean_nutrient_concentration_ci, 'R5', '2023', c('leaves', 'stems'), 'R5 leaves & stems')
element_barchart_df_ggplot(mean_nutrient_concentration_ci, 'R8', '2023', 'pods', 'R8 pods')
=======
#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

### Finding nutrient yield

# Note: Here we want to find the nutrient yield by calculating the 
# 1 - Tissue nutrient yield (g/m2) = Tissue nutrient concentration (g/g) * Tissue yield (g/m2)
# 2 - Find mean for each nutrient concentration and nutrient yield by treatment, variety, tissue_type, development_stage, year, block_id


nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")

# Loop through each nutrient and calculate nutrient yield
for (nutrient in nutrient_columns) {
  combined_data_df[[paste0(nutrient, "_yield")]] <- combined_data_df[[nutrient]] * combined_data_df$tissue_yield
}

# Calculate Mean for Nutrient Yield
nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")
yield_columns <- paste0(nutrient_columns, "_yield")

mean_nutrient_yield <- combined_data_df %>%
  group_by(treatment, variety, tissue_type, development_stage, year) %>%
  summarise(across(all_of(yield_columns), mean, na.rm = TRUE))


#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
##### Nutrient yield Relative Change
nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")
yield_columns <- paste0(nutrient_columns, "_yield")

mean_nutrient_yield_ci <- combined_data_df %>%
  group_by(treatment, variety, tissue_type, development_stage, year) %>%
  summarise(across(all_of(yield_columns), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE), 
                        n = ~sum(!is.na(.))), 
                   .names = "{col}_{fn}"))

# Write the treatment_cumulative_flow data frame to a CSV file
write.csv(mean_nutrient_yield_ci, "mean_nutrient_yield.csv", row.names = FALSE)

# Display confirmation message
cat("Data successfully written to mean_nutrient_yield.csv\n")

# Define a function to calculate relative change with confidence intervals
calculate_relative_change_yield <- function(big_df_subset, yield_columns) {
  diff_list <- by(
    big_df_subset,
    list(big_df_subset$development_stage, big_df_subset$year, big_df_subset$tissue_type, big_df_subset$variety),
    function(x) {
      ambient <- x[x$treatment == 'ambient', , drop = FALSE]
      elevated <- x[x$treatment == 'elevated', , drop = FALSE]
      
      if (nrow(ambient) == 0 | nrow(elevated) == 0) return(NULL)
      
      relative_change <- lapply(yield_columns, function(yc) {
        yield_mean_ambient <- paste0(yc, "_mean")
        yield_mean_elevated <- paste0(yc, "_mean")
        yield_sd_ambient <- paste0(yc, "_sd")
        yield_sd_elevated <- paste0(yc, "_sd")
        
        ambient_mean <- ambient[[yield_mean_ambient]]
        elevated_mean <- elevated[[yield_mean_elevated]]
        ambient_sd <- ambient[[yield_sd_ambient]]
        elevated_sd <- elevated[[yield_sd_elevated]]
        
        if (is.na(ambient_mean) | is.na(elevated_mean) | is.null(ambient_mean) | is.null(elevated_mean)) {
          return(list(relative_change = NA, ci_lower = NA, ci_upper = NA))
        }
        
        # Calculate relative change
        relative_change <- (elevated_mean - ambient_mean) / ambient_mean * 100
        
        # Error propagation for 95% confidence intervals
        ratio_se <- relative_change * sqrt((elevated_sd / elevated_mean)^2 + (ambient_sd / ambient_mean)^2)
        ci_lower <- relative_change - 1.96 * ratio_se
        ci_upper <- relative_change + 1.96 * ratio_se
        
        list(
          relative_change = relative_change,
          ci_lower = ci_lower,
          ci_upper = ci_upper
        )
      })
      
      # Combine results into a data frame
      relative_changes <- sapply(relative_change, function(rc) rc$relative_change)
      ci_lowers <- sapply(relative_change, function(rc) rc$ci_lower)
      ci_uppers <- sapply(relative_change, function(rc) rc$ci_upper)
      
      data.frame(
        nutrient = yield_columns,
        relative_change = relative_changes,
        ci_lower = ci_lowers,
        ci_upper = ci_uppers,
        development_stage = x$development_stage[1],
        year = x$year[1],
        tissue_type = x$tissue_type[1],
        variety = x$variety[1]
      )
    }
  )
  
  do.call(rbind, diff_list)
}

##### plot for nutrient yield
nutrient_yield_barchart <- function(
    big_df,
    development_stage,
    year,
    tissue,
    yield_columns,
    title = NULL,
    output_folder = NULL
) {
  # Filter data
  rows_to_keep <- 
    big_df$development_stage %in% development_stage & 
    big_df$year %in% year & 
    big_df$tissue_type %in% tissue
  
  big_df_subset <- big_df[rows_to_keep, ]
  
  # Calculate relative change
  diff_df <- calculate_relative_change_yield(big_df_subset, yield_columns)
  
  if (is.null(diff_df)) {
    stop("No valid data for the selected combination.")
  }
  
  diff_df$nutrient <- factor(diff_df$nutrient, levels = yield_columns)
  
  # Update variety labels to uppercase
  diff_df$variety <- factor(diff_df$variety, 
                            levels = c("hs", "loda"), 
                            labels = c("HS", "Loda"))
  
  # Create the bar chart
  p <- ggplot(diff_df, aes(x = nutrient, y = relative_change, fill = variety)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), 
             width = 0.7, 
             color = "black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, 
                  position = position_dodge(0.8), 
                  color = "black", 
                  size = 1.2) + 
    scale_fill_manual(values = c("#009999", "#FF6633"), name = "Variety") +  # Custom palette
    labs(
      title = title,
      x = "Nutrients",
      y = "Relative Change in Nutrient Yield (%)",
      fill = "Variety"  # Capitalize legend title
    ) +
    theme_classic2() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate x-axis text
      axis.ticks.x = element_blank(),  # Remove x-axis tick marks
      axis.text.y = element_text(size = 20, colour = "black"),
      axis.title = element_text(size = 24),
      axis.title.x = element_text(vjust = -2),
      axis.title.y = element_text(vjust = 3),
      plot.title = element_text(size = 36),
      strip.text = element_text(size = 24),
      legend.position = "right",
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 24),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.spacing = unit(1.5, "lines")
    ) +
    facet_wrap(~year)  # Group by year
  
  print(p)
}


# Define parameters
yield_columns <- paste0(c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu"), "_yield")
development_stage <- c("R8")  # Specify the development stages of interest
tissue_type <- c("pods")  # Example tissue type
year <- c(2022, 2023)  # Specify the years of interest

# Call the nutrient yield bar chart function
nutrient_yield_barchart(
  big_df = mean_nutrient_yield_ci,
  development_stage = development_stage,
  year = year,
  tissue = tissue_type,
  yield_columns = yield_columns,
  #title = "Nutrient Yield Relative Change"
)
>>>>>>> Stashed changes

#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
###
#---------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#


### Statistics 

### Nutrient concentration



### NUtrient yield
# Define the nutrient columns
nutrient_columns <- c("N", "P", "K", "Ca", "Mg", "S", "Zn", "Fe", "Mn", "Cu")
yield_columns <- paste0(nutrient_columns, "_yield")
results_list <- list()  # To store results for each nutrient

# Loop over each nutrient yield column
for (nutrient in yield_columns) {
  # Formula for the mixed model
  formula <- as.formula(paste(nutrient, "~ treatment * variety * tissue_type + (1|block_id) + (1|year)"))
  
  lmerModel <- lmer(as.formula(paste(el, "~ treatment + variety + (1|year) + (1|block_id) + 
                                     (1|treatment:variety) + (1|block_id:treatment)")),
                    data = tissue_data)
  # Fit the mixed-effects model
  lmer_model <- lmer(formula, data = combined_data_df)
  
  # Extract the summary of the model- p-values
  model_summary <- summary(lmer_model)
  
  # Extract the fixed effects (including interactions) and p-values
  effect_names <- rownames(model_summary$coefficients)
  fixed_effects <- model_summary$coefficients[, "Estimate"]
  p_values <- model_summary$coefficients[, "Pr(>|t|)"]
  
  # Perform ANOVA on the mixed-effects model to get F-values and p-values for the fixed effects
  anova_results <- anova(lmer_model)  # ANOVA for the mixed-effects model
  
  # Combine fixed effects and ANOVA results
  results_df <- data.frame(
    nutrient = nutrient,
    effect = effect_names,
    fixed_effects = fixed_effects,
    fixed_effects_p_values = p_values
  )
  
  # Add ANOVA results (F-values and p-values) to the results dataframe
  anova_effects <- rownames(anova_results)
  anova_f_values <- anova_results[, "F value"]
  anova_p_values <- anova_results[, "Pr(>F)"]
  
  # Create a data frame for ANOVA results
  anova_df <- data.frame(
    nutrient = nutrient,
    effect = anova_effects,
    anova_f_values = anova_f_values,
    anova_p_values = anova_p_values
  )
  
  # Merge fixed effects and ANOVA results into a single data frame for the current nutrient
  combined_results <- merge(results_df, anova_df, by = c("nutrient", "effect"), all = TRUE)
  
  # Add the combined results to the list
  results_list[[nutrient]] <- combined_results
}

# Combine all the results into a single data frame
final_results_df <- do.call(rbind, results_list)

# View the results
print(final_results_df)

# Prompt the user to select a folder to save the results
folder <- choose.dir(default = "", caption = "Select folder to save results")

# Check if folder selection was successful
if (!is.na(folder) && folder != "") {
  # Create a file path for the results CSV file
  file_path <- file.path(folder, "mixed_effects_nutrient_results.csv")
  
  # Save the results to the CSV file
  write.csv(final_results_df, file_path, row.names = FALSE)
  message("Results saved successfully in: ", file_path)
} else {
  message("Folder selection was canceled. Results not saved.")
}


# Model without year as a fixed effect
model_without_year <- lmer(as.formula(paste(nutrient, "~ treatment * variety * tissue_type + (1|block_id) + (1|year)")), data = combined_data_df)

# Model with year as a fixed effect
model_with_year <- lmer(as.formula(paste(nutrient, "~ treatment * variety * tissue_type + year + (1|block_id)")), data = combined_data_df)

# Perform a likelihood ratio test (LRT) to compare the two models
anova(model_without_year, model_with_year)
### Biomass


# Convert block_id, year, and tissue_type to factors
combined_mean_nutrients$block_id <- as.factor(combined_mean_nutrients$block_id)
combined_mean_nutrients$year <- as.factor(combined_mean_nutrients$year)
combined_mean_nutrients$tissue_type <- as.factor(combined_mean_nutrients$tissue_type)
combined_mean_nutrients$development_stage <- as.factor(combined_mean_nutrients$development_stage)


#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

### Final graphs

lmerout_j <- lmer(J_at_25_mean ~ treatment + variety + treatment:variety + (1|block_id) +
                    (1|block_id:treatment), mean_aci_parameters)
anova(lmerout_j)
rand(lmerout_j)
summary(lmerout_j)

formula <- as.formula(paste(nutrient, "~ treatment * variety + (1|block_id) + (1|year)"))

lmerModel <- lmer(as.formula(paste(el, "~ treatment + variety + (1|year) + (1|block_id) + 
                                     (1|treatment:variety) + (1|block_id:treatment)")),
                  data = tissue_data)