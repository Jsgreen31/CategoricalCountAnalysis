##
# Script to Analyze categorical count data using contingency tables
# Author: J. Green
# Date: Apr 26 2025
# NOTE: Script was modified from an original script by Wendy Monk

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(tibble)
library(car)
library(emmeans) 
library(rcompanion) 
library(epitools)
library(DescTools)
library(vcd)
library(gmodels)
library(lubridate)

# Step 1: Wrangle the data and organize it into Contingency Tables

# Load species x site matrix
full_species_matrix <- read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/Rarefaction Analysis/speciesXsitematrix.csv")  

# Convert to long format and filter detected species
species_long <- full_species_matrix %>%
  pivot_longer(cols = -location, names_to = "species", values_to = "presence") %>%
  filter(presence > 0)

# Load sampling units for each design
grid_sites <- read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/Sample Networks/Systematic Grid Plots.csv")
targeted_sites <- read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/Sample Networks/Targeted Plots.csv")
random_strat_sites <- read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/Sample Networks/random_stratified_samples.csv")

# Compile valid sites in my study area (exclude those outside)
valid_sites <- bind_rows(grid_sites, targeted_sites, random_strat_sites) %>% distinct(location)

#Load detection data
wildtrax_detections <- bind_rows(
  read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/HawkEars Validation:Analysis/wildtrax_validation_2021.csv") %>%
    mutate(location = as.character(location)),
  read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/HawkEars Validation:Analysis/wildtrax_validation_2022.csv") %>%
    mutate(location = as.character(location)),
  read.csv("/Users/Jorda/Desktop/Data/Cleaned Data/HawkEars Validation:Analysis/wildtrax_validation_2024.csv") %>%
    mutate(location = as.character(location))
)

# Filter detections to valid study sites 
filtered_detections <- wildtrax_detections %>%
  filter(location %in% valid_sites$location)


# Convert location columns to character vectors
grid_locations <- as.character(grid_sites$location)
targeted_locations <- as.character(targeted_sites$location)
random_strat_locations <- as.character(random_strat_sites$location)

# Assign sampling designs and count number of unique sites detected per species
species_site_counts <- filtered_detections %>%
  mutate(
    grid = location %in% grid_locations,
    targeted = location %in% targeted_locations,
    random_strat = location %in% random_strat_locations
  ) %>%
  pivot_longer(cols = c(grid, targeted, random_strat), names_to = "sampling_design", values_to = "included") %>%
  filter(included) %>%
  group_by(species_code, sampling_design) %>%
  summarise(num_sites_detected = n_distinct(location), .groups = "drop")

# Count total detections per species per sampling design
species_detection_counts <- filtered_detections %>%
  mutate(
    grid = location %in% grid_locations,
    targeted = location %in% targeted_locations,
    random_strat = location %in% random_strat_locations
  ) %>%
  pivot_longer(cols = c(grid, targeted, random_strat), names_to = "sampling_design", values_to = "included") %>%
  filter(included) %>%
  group_by(species_code, sampling_design) %>%
  summarise(num_detections = n(), .groups = "drop")

# Define total sampling units for each design (if nescessary)
total_grid_sites <- 118
total_targeted_sites <- 53
total_random_strat_sites <- 176

# Filter for selected species
selected_species <- c("CONI", "EWPW", "EAWP", "GCKI", "OVEN", "CSWA", "HETH" , "WTSP" , "SOSP")

filtered_species_site_counts <- species_site_counts %>% filter(species_code %in% selected_species)
filtered_species_detection_counts <- species_detection_counts %>% filter(species_code %in% selected_species)

# Add numSampleUnits to the dataset based on sampling design
filtered_species_site_counts <- filtered_species_site_counts %>%
  mutate(numSampleUnits = case_when(
    sampling_design == "Systematic Grid" ~ total_grid_sites,
    sampling_design == "Targeted" ~ total_targeted_sites,
    sampling_design == "Stratified" ~ total_random_strat_sites
  ))

filtered_species_detection_counts <- filtered_species_detection_counts %>%
  mutate(numSampleUnits = case_when(
    sampling_design == "Grid" ~ total_grid_sites,
    sampling_design == "Targeted" ~ total_targeted_sites,
    sampling_design == "Stratified" ~ total_random_strat_sites
  ))

# Arrange data into contingency table format (num sites)
contTableDataSites <- xtabs(num_sites_detected ~ species_code + sampling_design, data = filtered_species_site_counts)

# Create contingency table (sites)
CrossTable(contTableDataSites,
           digits = 3,
           expected = TRUE,
           prop.r = TRUE,
           prop.c = TRUE,
           chisq = TRUE,
           sresid = TRUE,
           format = "SPSS")

# Arrange data into contingency table format (num detections)
contTableDataDetections <- xtabs(num_detections ~ species_code + sampling_design, 
                                 data = filtered_species_detection_counts)

# Create contingency table (detections)
CrossTable(contTableDataDetections,
           digits = 3,
           expected = TRUE,
           prop.r = TRUE,
           prop.c = TRUE,
           chisq = TRUE,
           sresid = TRUE,
           format = "SPSS")

# It seems that there is a statistically significant difference in the sampling design 
# used relative to the # of sites and detections for each species! 

# Can explore the patterns of standardised residuals to look at which cross-classifications deviate from the expected
# values and drive the lack of independence.
# The formula is (observed - expected)/sqrt(expected).

# The units are in standard deviations, so a residual greater than 2 or less than -2
# represents a departure significant at the 95% level.

# We can plot the data using a mosaic plot i.e. a rectangle is divided into several rectangles
# whose areas are proportional to the counts (one rectangle for each cell)

# Note that it will shade any cells where you have significantly greater (blue) or significantly
# fewer (red) observations than expected

## Expected (num sites)
vcd::strucplot(contTableDataSites, shade = TRUE,
               labeling_args=list(set_varnames=c(species_code="Species",
                                                sampling_design = "Sample Design")),
               type = "expected")

## Observed (num sites)
vcd::strucplot(contTableDataSites, shade = TRUE,
               labeling_args=list(set_varnames=c(species_code ="Species",
                                                 sampling_design = "Sample Design")),
               type = "observed")

# Expected (num detections)
vcd::strucplot(contTableDataDetections, shade = TRUE,
               labeling_args=list(set_varnames=c(species_code="Species",
                                                 sampling_design = "Sample Design")),
               type = "expected")

## Observed (num detections)
vcd::strucplot(contTableDataDetections, shade = TRUE,
               labeling_args=list(set_varnames=c(species_code="Species",
                                                 sampling_design = "Sample Design")),
               type = "observed")

## Or via an association plot i.e. where the observed and expected frequencies so the
# height of a rectangle is proportional to the residual and the width is
# proportional to the square root of the expected counts.

# The rectangles for each row in the table are positioned relative to a baseline
# representing independence shown by a dotted line. Cells with observed > expected
# frequency rise above the line; cells that contain less than the expected frequency fall below it

# num sites plot
vcd::assoc(contTableDataSites,
           shade = TRUE,
           compress = FALSE,
           labeling_args=list(set_varnames=c(species_code="Species",
                                             sampling_design = "Sample Design")))

# num detections plot
vcd::assoc(contTableDataDetections,
           shade = TRUE,
           compress = FALSE,
           labeling_args=list(set_varnames=c(species_code="Species",
                                             sampling_design = "Sample Design")))

# Great now we have the results plotted in an intuitive way to interpret
# But we should compare to see if scaling locations with less sampling events influences results

# Determine number of samples (restrict to a single year of data if sampled multiple years. Max 7)
# First let's determine how many recordings were processed at each site and how many were sampled in > 1 year


# Calculate number of total visits per site (unique recording days)

# Extract year from recording_datetime
filtered_detections <- filtered_detections %>%
  mutate(year = year(recording_datetime))

# Count distinct visits per site per year
visits_per_site <- filtered_detections %>%
  distinct(location, recording_datetime, year) %>%
  count(location, year, name = "visits")

# Identify the year with the most visits per site
top_year_per_site <- visits_per_site%>%
  group_by(location) %>%
  slice_max(visits, n = 1, with_ties = FALSE) %>%
  ungroup()

# Filter original data to only keep entries from the top year for each site
cleaned_detections <- filtered_detections %>%
  inner_join(top_year_per_site, by = c("location", "year"))

# Summarize the total number of visits after cleaning
total_visits <- cleaned_detections %>% 
  distinct(location, recording_datetime) %>% 
  count(location, name = "visits")

# now with the number of visits sorted, let's clean the detection data to only include unique detections (1SPT)

# Modify detection data to ensure each species is counted once per recording
unique_detections <- filtered_detections %>%
  distinct(location, recording_datetime, species_code)

# Now calculate the total number of detections for each species at each site
scaled_species_detection_counts <- unique_detections %>%
  group_by(species_code, location) %>%
  summarise(num_detections = n(), .groups = "drop")

# Calculate the proportion of visits with detections
detection_proportion <- scaled_species_detection_counts %>%
  left_join(total_visits, by = "location") %>%
  mutate(proportion_detections = num_detections / visits)

# with the detection proportions sorted, let's choose a common surey effort by determining the mean # visits

mean(total_visits$visits)

# Since we have a mean of 6.19, let's use 6 samples as our reference for survey effort

# Merge total visits into proportion deetection
scaled_df <- detection_proportion %>%
  mutate(
    scaled_detection = num_detections / visits *  6
  )

# Now that we have our detections scaled, let's add them up based on a per species/sample design basis
# and re-run the contingency analysis to see if there are differences

# first let's sort out the sample designs and aggregate the detections accordingly

grid_sites_clean <- grid_sites %>%
  mutate(location = as.character(location))

targeted_sites_clean <- targeted_sites %>%
  mutate(location = as.character(location))

random_strat_sites_clean <- random_strat_sites %>%
  mutate(location = as.character(location), design = "Random Stratified") %>%
  select(location, design)

design_df <- bind_rows(grid_sites_clean, targeted_sites_clean, random_strat_sites_clean)







# Combine site design info with scaled detection data
scaled_df_design <- scaled_df %>%
  mutate(
    sampling_design = case_when(
      location %in% grid_sites_clean$location ~ "Grid",
      location %in% targeted_sites_clean$location ~ "Targeted",
      location %in% random_strat_sites_clean$location ~ "Stratified",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(species_code %in% selected_species, !is.na(sampling_design))

# Aggregate scaled detections by species and sampling design
aggregated_scaled_detections <- scaled_df_design %>%
  group_by(species_code, sampling_design) %>%
  summarise(total_scaled_detections = sum(scaled_detection), .groups = "drop")

# Create contingency table 
contTableDataScaled <- xtabs(total_scaled_detections ~ species_code + sampling_design, 
                             data = aggregated_scaled_detections)

# View CrossTable of scaled detections
CrossTable(contTableDataScaled,
           digits = 3,
           expected = TRUE,
           prop.r = TRUE,
           prop.c = TRUE,
           chisq = TRUE,
           sresid = TRUE,
           format = "SPSS")


# Mosaic plot (scaled detections - observed)
vcd::strucplot(contTableDataScaled,
               shade = TRUE,
               labeling_args = list(set_varnames = c(species_code = "Species",
                                                     sampling_design = "Sample Design")),
               type = "observed")

# Association plot (scaled detections)
vcd::assoc(contTableDataScaled,
           shade = TRUE,
           compress = FALSE,
           labeling_args = list(set_varnames = c(species_code = "Species",
                                                 sampling_design = "Sample Design")))

# Great now we can compare the scaled to unscaled results to determine how bias in survey effort
# is handled!


























