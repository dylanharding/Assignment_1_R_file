###########################################
######## BINF6210 Assignment 1 ############
########    Dylan Harding      ############
###########################################

# Thank you to.....
# Sameh Mohamed
# Mirza Ahmadi
# For their support/guidance on this project

## Introduction --------

# Let's start by clearing the global environment
rm(list = ls())

## Please set working directory to source file location ##

# Below are the packages I used for this analysis. Please install as needed.
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("viridis")
# install.packages("ggplot2")
# install.packages("maps")
# install.packages("mapdata")
# install.packages("sf")
# install.packages("stringr")
library(terra)
library(tidyverse)
library(ggplot2)
library(vegan)
library(viridis)
library(mapdata)
library(maps)
library(sf)
library(stringr)
library(styler)
theme_set(theme_light())
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")

## Using BOLD API to download specimen data for the families Salmonidae and Caligidae
# dfBOLD_salmonidae = read_tsv("http://www.boldsystems.org/index.php/API_Public/specimen?taxon=Salmonidae&format=tsv")
# dfBOLD_caligidae = read_tsv("http://www.boldsystems.org/index.php/API_Public/specimen?taxon=caligidae&format=tsv")

## Writing TSV files to harddisk so original data is saved
# write_tsv(dfBOLD_salmonidae, "salmonidae_BOLD_data.tsv")
# write_tsv(dfBOLD_caligidae, "caligidae_BOLD_data.tsv")

# Reading Salmonidae, Caligidae, and salmon farm TSV files stored locally. (See code citations for salmon farm dataframe)

dfBOLD_salmonidae <- read_tsv("../data/salmonidae_BOLD_data.tsv")
dfBOLD_caligidae <- read_tsv("../data/caligidae_BOLD_data.tsv")
dfsalmon_farms <- read_tsv("../data/PUBLIC_CanadianSalmonFarms_Sheet1.tsv")


?read_csv()

# Checking that each data frame was imported successfully
class(dfBOLD_salmonidae)
summary(dfBOLD_salmonidae)
class(dfBOLD_caligidae)
summary(dfBOLD_caligidae)
class(dfsalmon_farms)
summary(dfsalmon_farms)

# Filter for only Canadian specimens in BOLD dataframes
dfsalmonidae_Can <- dfBOLD_salmonidae %>%
  filter(country == "Canada")
dfcaligidae_Can <- dfBOLD_caligidae %>%
  filter(country == "Canada")

# The code below reduces the redundancy and makes the code more cleaner.This apporach can help the data filtering step to be more clear and concise, ultimately making the code more efficient. 
dfsalmonidae_Can <- filter(dfBOLD_salmonidae, country == "Canada")
dfcaligidae_Can <- filter(dfBOLD_caligidae, country == "Canada")

# Checking column names from the salmon farm database
names(dfsalmon_farms)

# Removing all spaces and periods from column names and converting all to lower case for easier manipulation
names.original <- as.vector(names(dfsalmon_farms))
names.original
class(names.original)

names.edited <- str_replace_all(
  string = names.original,
  pattern = "[. ]", replacement = "_"
) %>%
  str_to_lower()

# Assigning edited names to column headers and checking
names(dfsalmon_farms) <- names.edited
names(dfsalmon_farms)

# Changing names of latitude/longitude headers for consistency with BOLD data
dfsalmon_farms <- dfsalmon_farms %>%
  rename(
    lat = latitude,
    lon = longitude,
  )

# Checking for unique farm types and numbers of entries
unique(dfsalmon_farms$farm_type)
length(dfsalmon_farms$farm_type)

# Refining farm types, removing unknowns
dfsalmon_farms <- dfsalmon_farms %>%
  filter(!farm_type %in% c("N/A", "Abandoned farm cages?", "TBD")) %>%
  mutate(farm_type = str_replace_all(farm_type, pattern = "Seafarm", replacement = "Sea Farm"))

# Rechecking cleaner farm type column
unique(dfsalmon_farms$farm_type)
length(dfsalmon_farms$farm_type)

# Remove values unneeded for next section
rm(names.edited, names.original)

## Accumulation Curve --------

# First, I'll make an accumulation curve using BIN ID and Region data for the family Salmonidae. This will give us an idea of the sampling completeness in Canada.

# Counting number of salmonidae bins/region
dfsalmonidae_bins_x_region <- dfsalmonidae_Can %>%
  filter(!is.na(bin_uri)) %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  count(bin_uri
)

# Restructuring data to be used in the specaccum() function
dfsalmonidae_spread <- pivot_wider(
  data = dfsalmonidae_bins_x_region,
  names_from = bin_uri,
  values_from = n
)

# Converting NAs --> 0s
dfsalmonidae_spread[is.na(dfsalmonidae_spread)] <- 0

# Moving region column to row names
dfsalmonidae_spread <- dfsalmonidae_spread %>%
  remove_rownames() %>%
  column_to_rownames(var = "region")

accum_curve <- specaccum(dfsalmonidae_spread)

## FIGURE 1 #######
# Plotting salmonidae accumulation curve
plot(accum_curve,
  xlab = "Regions Sampled",
  ylab = "BIN Richness",
  main = str_wrap("Fig 1. Accumulation Curve for Salmonidae in Canada",
    width = 45
  ),
  col = "blue"
)

# The plateau of the curve towards the right side of the X-axis indicates good sampling completeness, meaning we expect our data to be a relatively good representation of salmonidae diversity in Canada.

# Remove all objects not needed for next section
rm(accum_curve, dfsalmonidae_bins_x_region, dfsalmonidae_spread)

## Geographic Distributions-------

# The purpose of this section is to plot geographic data (in this case lat & lon), for both Salmonidae and Caligidae on the map. Although not statistical in nature, this should give us a visual idea of geographic correlation between the two families in Canada.

# Filtering BOLD dataframesfor relevant lat/lon data for each family, removing NAs.
dfsalmonidae_filtered <- dfsalmonidae_Can %>%
  select(family_name, lat, lon) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(family_name))

dfcaligidae_filtered <- dfcaligidae_Can %>%
  select(family_name, lat, lon) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(family_name))

# This code block below will reduce the redundancy the code above and is simplified fitlering.
clean_data <- function(data) {
  data %>%
    select(family_name, lat, lon) %>%
    filter(complete.cases(family_name, lat, lon))
}

dfsalmonidae_filtered <- clean_data(dfsalmonidae_Can)
dfcaligidae_filtered <- clean_data(dfcaligidae_Can)


# Filter salmon farm data for lat/lon, temporarily creating new column "family_name" for consistency with BOLD data
dfsalmon_farms_filtered <- dfsalmon_farms %>% 
  select(lat, lon) %>% 
  filter(!is.na(as.numeric(lat))) %>% 
  filter(!is.na(as.numeric(lon))) %>% 
  mutate(family_name = "Salmon Farm")

# Checking range of lat/lon data
# Lat range for Canada ~= 41.5 to 83.5
# Lon range for Canada ~= -141 to -52.5
range(dfsalmonidae_filtered$lat)
range(dfcaligidae_filtered$lat)
range(dfsalmon_farms_filtered$lat)
range(dfsalmonidae_filtered$lon)
range(dfcaligidae_filtered$lon)
range(dfsalmon_farms_filtered$lon)

# This form of code outputs well-formated results and is clearer to understand.The code would more maintainable and is not repetitive.
check_lat_lon_range <- function(data, name) {
  cat(paste("\n", name, " Latitude Range: ", range(data$lat), "\n"))
  cat(paste(name, " Longitude Range: ", range(data$lon), "\n"))
}

check_lat_lon_range(dfsalmonidae_filtered, "Salmonidae")
check_lat_lon_range(dfcaligidae_filtered, "Caligidae")
check_lat_lon_range(dfsalmon_farms_filtered, "Salmon Farms")


#The block of codes below achieve generating a scatter plot of latitide and longitude rang of the given species. The scatter plot will visualize and detect data biases and outliers in the data set. This plot can be used a a data exploration step to see if there are outliers or any overlapping data. This data exploration step can be used to give an idea of what visual data form could be used to answer the "research question" 

dfsalmonidae_filtered$lat <- as.numeric(dfsalmonidae_filtered$lat)
dfsalmonidae_filtered$lon <- as.numeric(dfsalmonidae_filtered$lon)
dfcaligidae_filtered$lat <- as.numeric(dfcaligidae_filtered$lat)
dfcaligidae_filtered$lon <- as.numeric(dfcaligidae_filtered$lon)
dfsalmon_farms_filtered$lat <- as.numeric(dfsalmon_farms_filtered$lat)
dfsalmon_farms_filtered$lon <- as.numeric(dfsalmon_farms_filtered$lon)


dfsalmonidae_filtered$Dataset <- "Salmonidae"
dfcaligidae_filtered$Dataset <- "Caligidae"
dfsalmon_farms_filtered$Dataset <- "Salmon Farms"


combined_data <- rbind(dfsalmonidae_filtered, dfcaligidae_filtered, dfsalmon_farms_filtered)


ggplot(combined_data, aes(x = lon, y = lat, color = Dataset)) +
  geom_point(alpha = 0.8, size = 6, shape = 16) +
  labs(title = "Geographic Distribution of Species and Farms in Canada",
       x = "Longitude",
       y = "Latitude",
       color = "Dataset") +
  xlim(-141, -52.5) +  
  ylim(41.5, 83.5) +   
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  
    legend.position = "right"                
  )
# Satisfied with lat/lon ranges, data can be combined and converted to sf object
dfcombined_lonlat_data <- rbind(dfsalmonidae_filtered, dfcaligidae_filtered, dfsalmon_farms_filtered)

map_data <- st_as_sf(dfcombined_lonlat_data, coords = c("lon", "lat"), crs = 4326)

# Checking for sf data
class(map_data)

## FIGURE 2 #######
# Plotting distribution map with prepared salmonidae and caligidae sf data
ggplot() +
  geom_sf(data = st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)), fill = "lightgrey") +
  geom_sf(data = map_data, aes(color = family_name), size = 1) +
  facet_wrap(~ family_name) +  # Create separate panels for each family
  scale_color_manual(values = c("Caligidae" = "blue", "Salmonidae" = "red", "Salmon Farm" = "green")) +
  labs(
    title = str_wrap("Fig 2. Geographic distributions of Caligidae, Salmonidae, and Salmon Farms in Canada", width = 60),
    color = NULL
  ) +
  coord_sf(xlim = c(-145, -55), ylim = c(40, 75)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", colour = "black")
  )

dfsalmon_farms_filtered <- dfsalmon_farms %>%
  
# The map demonstrates a geographic overlap between both families, particularly on the west coast of Canada. One stand-out is the notably fewer specimen samples for Caligidae compared to Salmonidae.

# Remove all objects not needed for next section
rm(map_data, dfsalmonidae_filtered, dfcaligidae_filtered, dfcombined_lonlat_data, dfsalmon_farms_filtered)

## Salmon Farms--------

# As a final point, and for discussion relevancy, in this section we will examine data concerning Canadian salmon farms, downloaded from a publicly available database (see citations). Specifically, we will plot a breakdown of salmon farm types using a stacked bar graph.

# Preparing data for graphing by grouping, total counts
dffarms_x_province <- dfsalmon_farms %>%
  group_by(region_or_province, farm_type) %>%
  summarize(total_farms = n(), .groups = "drop")

## FIGURE 3 #######
# Creating stacked bar graph
ggplot(dffarms_x_province, aes(x = region_or_province, y = total_farms, fill = farm_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = str_wrap("Fig 3. Total Number of Canadian Salmon Farms by Region/Province and Farm Type", width = 50),
    x = "Region/Province",
    y = "Total Number of Farms",
    fill = "Farm Types"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(fill = "white", colour = "black")
  ) +
  scale_fill_viridis_d(option = "E"
)

# Another potential way to convey the information presented in the stacked bar graph. The GLMM test below indicates that Sea Farm type is much more common in the given provinces. This test will account for the fact whether there is sampling bias as the Sea Farm type has shown disproportionately high count in comparison to other farm types. The statistical test also directs a new light to the data, by adding more context.
set.seed(123)  
dfsalmon_farms$farm_size <- runif(nrow(dfsalmon_farms), min = 1, max = 100)

dffarms_x_province <- dfsalmon_farms %>%
  group_by(region_or_province, farm_type) %>%
  summarize(
    total_farms = n(),
    avg_farm_size = mean(farm_size, na.rm = TRUE),
    .groups = "drop"
  )

library(lme4)

glmm_model <- glmer(total_farms ~ farm_type + (1 | region_or_province), data = dffarms_x_province, family = poisson)
summary(glmm_model)


# The resulting graph highlights the vast number of sea farms compared to the other types, a potentially important factor for salmon diversity, as we will discuss

###### END OF CODE ######