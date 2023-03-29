rm(list = ls())
library(tidyverse)
library(terra)
library(rnrfa)
library(lubridate)
library(measurements)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
source("source/deerSDM_misc_functions.R")

# ------------------------------------------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------------------------------------------
ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))
inaturalist <- read.csv('data-raw/inaturalist_deer_england.csv')
BDS <- read.csv('data-raw/BDS_distribution_data_2005_2016.csv')
cull <- read.csv('data-raw/cull_2017_present.csv')
RTA <- read.csv('data-raw/Highways Agency RTA data.xlsx - Sheet1.csv')
DVC <- read.csv('data-raw/Langbein_DVC_MAIN_Yr2010_2020_EngWal_v_Feb2023.csv') %>% 
     filter(CaseForestsAll_or_nonCaseFC != "LangbeinCaseForests")
GBC <- read.csv('data-raw/Noble_etal_2012_DeerDensityMaps2005-09.csv')
BBS <- read.csv('data-raw/BBS WHITE REF1674141171263692/BBS WHITE REF1674141171263692.csv') 
traffic_data_2020_2022 <- read.csv("data-raw/England_traffic_data_2020_2022.csv")
traffic_data_2010_2020 <- read.csv("data-raw/England_traffic_data_2010_2020.csv")
roads_distance <- rast("spatial_data/roads_distance_Eng60m.tif")
named_roads <- vect('spatial_data/England_named_roads.shp') %>% 
      terra::project(crs(km_projection))

# Species of interest
species_list <- c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi", "Cervus nippon", "Hydropotes inermis")

# ------------------------------------------------------------------------------------------------------------------------
# iNaturalist records
# ------------------------------------------------------------------------------------------------------------------------
# Filter and clean inaturalist
inaturalist_vect <- inaturalist %>% 
      filter(scientific_name %in% species_list,
             num_identification_disagreements == 0,
             # Remove observations with spatial resolution lower than 1km
             positional_accuracy < 1000) %>% 
      dplyr::select(scientific_name, latitude, longitude, observed_on, positional_accuracy, user_id) %>% 
      rename(species = scientific_name) %>% 
      vect(geom = c("longitude", "latitude"), crs = crs("epsg:4326")) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") %>% 
      vect(geom = c("x", "y"), crs =crs(km_projection), keepgeom = T) %>% 
      # Add effort, as distance to road
      extract(roads_distance, ., bind = T) 

grid <- rast(extent = ext(ROI), resolution = c(10, 10), crs = crs(km_projection)) %>% 
      crop(ROI)
values(grid) <- 1:nrow(values(grid))
names(grid) <- "grid_ID"

inaturalist_final_2004_2020 <- extract(grid, inaturalist_vect, bind = T) %>% 
      data.frame() %>% 
      group_by(grid_ID) %>% 
      mutate(n_users_100km2 = length(unique(user_id)),
             distance_to_road = round(distance_to_road)) %>% 
      ungroup() %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection)) %>% 
      as.data.frame(geom = "XY") %>% 
      group_by(species, x, y, n_users_100km2, distance_to_road) %>% 
      summarize(presence = 1) %>% 
      ungroup()

# ------------------------------------------------------------------------------------------------------------------------
# British Deer Survey records
# ------------------------------------------------------------------------------------------------------------------------
BDS_long <- BDS %>% 
      dplyr::select(1:9, "SUR_05", "SUR_2010", "SUR_2016", everything()) %>% 
      pivot_longer(cols = 13:30, names_to = "var", values_to = "presence") %>% 
      separate(var, sep = "_", into = c("species", "year")) %>% 
      mutate(year = as.integer(year),
             year = replace(year, year == 05, 2005),
             species = replace(species, species == "FAL", "Dama dama"),
             species = replace(species, species == "ROE", "Capreolus capreolus"),
             species = replace(species, species == "RED", "Cervus elaphus"),
             species = replace(species, species == "MUN", "Muntiacus reevesi"),
             species = replace(species, species == "SIK", "Cervus nippon"),
             species = replace(species, species == "CWD", "Hydropotes inermis")) %>% 
      filter(species %in% species_list) %>% 
      dplyr::select("OS_Tile", "species", "year", "presence", "X_BottLeft", "Y_BottLeft", "X_Cent", "Y_Cent")

BDS_survey_overview <- BDS %>% 
      dplyr::select("OS_Tile", "SUR_05", "SUR_2010", "SUR_2016") %>% 
      pivot_longer(cols = 2:4, names_to = "x", values_to = "surveyed") %>% 
      separate(x, sep = "_", into = c("x", "year")) %>% 
      mutate(year = as.integer(year),
             year = replace(year, year == 05, 2005),
             surveyed = ifelse(surveyed == "Yes", 1, 0)) %>% 
      dplyr::select(-x)

# Get some kind of effort for BDS data, as number of times square was visited
BDS_effort <- BDS_survey_overview %>% 
      group_by(OS_Tile) %>% 
      summarize(n_years_surveyed = sum(surveyed))

BDS_merged <- merge(BDS_long, BDS_survey_overview, all.x = TRUE, by = c("OS_Tile", "year")) %>% 
      filter(surveyed == 1) %>% 
      mutate(presence = ifelse(presence == "Yes", 1, 0)) %>% 
      dplyr::select(-surveyed)

BDS_final_2005_2016 <- merge(BDS_merged, BDS_effort, all.x = TRUE, by = c("OS_Tile")) %>% 
      dplyr::select(-X_BottLeft, -Y_BottLeft) %>% 
      vect(geom = c("X_Cent", "Y_Cent"), crs = crs("epsg:27700"), keepgeom = FALSE) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") %>% 
      group_by(species, x, y, n_years_surveyed) %>% 
      summarize(years_present = sum(presence)) %>% 
      mutate(presence = ifelse(years_present > 0, 1, 0))

# ------------------------------------------------------------------------------------------------------------------------
# Deer culling records
# ------------------------------------------------------------------------------------------------------------------------
cull_final_2017_2020 <- cull %>% 
      mutate(Grid.Ref = gsub(" ", "", Grid.Ref, fixed=T),
             easting = rnrfa::osg_parse(Grid.Ref)$easting,
             northing = rnrfa::osg_parse(Grid.Ref)$northing) %>% 
      group_by(Grid.Ref, easting, northing, Species) %>% 
      summarize(count = n()) %>% 
      drop_na() %>% 
      mutate(Species = replace(Species, Species == "Roe Deer", "Capreolus capreolus"),
             Species = replace(Species, Species == "Red Deer", "Cervus elaphus"),
             Species = replace(Species, Species == "Fallow Deer", "Dama dama"),
             Species = replace(Species, Species == "Muntjac Deer", "Muntiacus reevesi"),
             Species = replace(Species, Species == "Sika Deer", "Cervus nippon"),
             Species = replace(Species, Species == "Chinese Water Deer", "Hydropotes inermis")) %>% 
      filter(Species %in% species_list) %>% 
      vect(geom = c("easting", "northing"), crs = crs("epsg:27700"), keepgeom = FALSE) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") %>% 
      mutate(presence = 1) %>% 
      rename(species = Species)

write.csv(cull_final_2017_2020, 'data-processed/cull_final_2017_2020.csv')

# ------------------------------------------------------------------------------------------------------------------------
# Game bag census 10km data
# ------------------------------------------------------------------------------------------------------------------------
extract_mean <- function(x) {
      numbers <- as.numeric(unlist(strsplit(gsub("[^0-9]+", " ", x), " ")))
      mean(numbers)
}

GBC_final <- GBC %>% 
      rowwise() %>% 
      mutate(across(2:5, extract_mean),
             easting = rnrfa::osg_parse(os_grid)$easting,
             northing = rnrfa::osg_parse(os_grid)$northing) %>% 
      ungroup() %>% 
      pivot_longer(cols = 2:5, names_to = "species", values_to = "NGC") %>% 
      drop_na() %>% 
      mutate(species = replace(species, species == "roe_deer", "Capreolus capreolus"),
             species = replace(species, species == "red_deer", "Cervus elaphus"),
             species = replace(species, species == "fallow_deer", "Dama dama"),
             species = replace(species, species == "muntjac", "Muntiacus reevesi"))

write.csv(GBC_final, 'data-processed/GBC_final_1961_2009.csv')

# ------------------------------------------------------------------------------------------------------------------------
# Deer density records -  literature
# ------------------------------------------------------------------------------------------------------------------------
deer_density <- read.csv('data-raw/deer_density_database.csv') %>% 
      filter(Country != 'Scotland', Deer.species != "Sika") %>% 
      dplyr::select(Author.s., URL, Deer.species, Lat, Long, Area, Area.unit, density.value..specific., units, Location) %>% 
      mutate(Area = ifelse(Area.unit == "ha", Area*0.01, Area)) %>% 
      rename(density = density.value..specific.) %>% 
      mutate(density = as.numeric(density),
             density = ifelse(units == "n / ha", density*0.01, density), 
             units = 'n / km2',
             Area.unit = 'km2') %>% 
      filter(!is.na(density)) %>% 
      rowwise() %>% 
      mutate(Lat_chd = ifelse(str_detect(Lat, "[A-Za-z]"),
                          str_extract_all(Lat, "\\(?[0-9,.]+\\)?")[[1]][1], NA),
             Lat_chm = ifelse(str_detect(Lat, "[A-Za-z]"),
                               str_extract_all(Lat, "\\(?[0-9,.]+\\)?")[[1]][2], NA),
             Long_chd = ifelse(str_detect(Long, "[A-Za-z]"),
                               str_extract_all(Long, "\\(?[0-9,.]+\\)?")[[1]][1], NA),
             Long_chm = ifelse(str_detect(Long, "[A-Za-z]"),
                              str_extract_all(Long, "\\(?[0-9,.]+\\)?")[[1]][2], NA),
             lat_dec = conv_unit(paste0(Lat_chd, " ", Lat_chm), from = 'deg_dec_min', to = 'dec_deg'),
             long_dec = conv_unit(paste0(Long_chd, " ", Long_chm), from = 'deg_dec_min', to = 'dec_deg'),
             Lat = as.numeric(gsub(",", "", ifelse(!is.na(Lat) & !is.na(lat_dec), lat_dec, Lat))),
             Long = as.numeric(ifelse(!is.na(Long) & !is.na(long_dec), long_dec, Long))) %>% 
      dplyr::select(1:10) %>% 
      vect(geom = c("Long", "Lat"), crs = crs("epsg:4326")) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") 
write.csv(deer_density, 'data-processed/deer_density_literature.csv')

# ------------------------------------------------------------------------------------------------------------------------
# Deer Vehicle Collisions - Jochen Langbein
# ------------------------------------------------------------------------------------------------------------------------
DVC_edited <- DVC %>% 
      dplyr::select("INC_year", "OSEasting", "OSNorthing", "DeerSpeciesIfKnown", "RoadNo_if_allocated") %>% 
      rename(year = INC_year,
             easting = OSEasting,
             northing = OSNorthing,
             species = DeerSpeciesIfKnown,
             Road_name = RoadNo_if_allocated) %>% 
      mutate(species = replace(species, species == "Fallow", "Dama dama"),
             species = replace(species, species == "ROE deer", "Capreolus capreolus"),
             species = replace(species, species == "Red deer", "Cervus elaphus"),
             species = replace(species, species == "Muntjac", "Muntiacus reevesi"),
             species = replace(species, species == "Sika", "Cervus nippon"),
             species = replace(species, species == "Chinese Water", "Hydropotes inermis"),
             Road_name = toupper(Road_name),
             Road_name = replace(Road_name, Road_name == "X-Z-ROAD NO NOT ALLOCATED", NA)) %>%  
      mutate(Road_name = str_squish(sub(" .*", "", Road_name))) %>% 
      filter(species %in% species_list) %>% 
      vect(geom = c("easting", "northing"), crs = crs("epsg:27700"), keepgeom = FALSE) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") %>% 
      group_by(Road_name, x, y, species) %>% 
      summarise(n_years_obs = n()) %>% 
      ungroup() %>% 
      mutate(index = row_number())

nrow(DVC_edited[is.na(DVC_edited$Road_name), ])

# Check if any points have identical coordinates, but different road names
ambiguous_road_names <- DVC_edited %>%
      filter(!is.na(Road_name)) %>% 
      group_by(x, y) %>% 
      mutate(n = length(unique(Road_name))) %>% 
      filter(n > 1)

# There are a few hundred. Set these to NA.
DVC_edited$Road_name[DVC_edited$index %in% ambiguous_road_names$index] <- NA

# Add road names to points that are spatially close to a named road
DVC_NA_names <- DVC_edited %>% 
      filter(is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

named_roads_points <- as.points(named_roads, multi = TRUE)

added_closest_named_road <- get_nearest_value(DVC_NA_names, named_roads_points, max_dist = 300, value = "ref")

DVC_names_added <- merge(DVC_edited, added_closest_named_road, all.x = T, by = c("x", "y")) %>% 
      mutate(Road_name = ifelse(is.na(Road_name) & !is.na(ref), ref, Road_name)) %>% 
      dplyr::select(-ref) 

nrow(DVC_names_added[is.na(DVC_names_added$Road_name), ])

# Check if there are any new ambiguous road names 
ambiguous_road_names_2 <- DVC_names_added %>%
      filter(!is.na(Road_name)) %>% 
      group_by(x, y) %>% 
      mutate(n = length(unique(Road_name))) %>% 
      filter(n > 1)

# Set ambiguous road names to NA again
DVC_names_added$Road_name[DVC_names_added$index %in% ambiguous_road_names_2$index] <- NA

# Export to check result in QGIS
# writeVector(named_roads_points, "named_roads_points.shp", overwrite = T)
# writeVector(DVC_edited %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE), "spatial_data/all_DVC_points.shp", overwrite = T)
# writeVector(DVC_names_added %>% filter(!is.na(Road_name)) %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE), "spatial_data/DVC_names_added.shp", overwrite = T)
# writeVector(DVC_names_added %>% filter(is.na(Road_name)) %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE), "spatial_data/DVC_remaining_NA.shp", overwrite = T)

# Some points have NA road names, even though they are very close to points with the correct road name.
# Carry over the correct road names, if points are 250m or less apart. Try running more than once.
DVC_non_NA_roads <- DVC_names_added %>% 
      filter(!is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

DVC_NA_roads <- DVC_names_added %>% 
      filter(is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

added_closest_named_DVC <- get_nearest_value(DVC_NA_roads, DVC_non_NA_roads, max_dist = 300, value = "Road_name", add_unique_value_name = TRUE)

DVC_names_added_2 <- merge(DVC_names_added, added_closest_named_DVC, all.x = T, by = c("x", "y")) %>% 
      mutate(Road_name = ifelse(is.na(Road_name) & !is.na(Road_name_new), Road_name_new, Road_name)) %>% 
      dplyr::select(-Road_name_new) 

nrow(DVC_names_added_2[is.na(DVC_names_added_2$Road_name), ])

DVC_non_NA_roads_2 <- DVC_names_added_2 %>% 
      filter(!is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

DVC_NA_roads_2 <- DVC_names_added_2 %>% 
      filter(is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

added_closest_named_DVC_2 <- get_nearest_value(DVC_NA_roads_2, DVC_non_NA_roads_2, max_dist = 300, value = "Road_name", add_unique_value_name = TRUE)

DVC_names_added_3 <- merge(DVC_names_added_2, added_closest_named_DVC_2, all.x = T, by = c("x", "y")) %>% 
      mutate(Road_name = ifelse(is.na(Road_name) & !is.na(Road_name_new), Road_name_new, Road_name)) %>% 
      dplyr::select(-Road_name_new) 

nrow(DVC_names_added_3[is.na(DVC_names_added_3$Road_name), ])

DVC_non_NA_roads_3 <- DVC_names_added_3 %>% 
      filter(!is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

DVC_NA_roads_3 <- DVC_names_added_3 %>% 
      filter(is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

added_closest_named_DVC_3 <- get_nearest_value(DVC_NA_roads_3, DVC_non_NA_roads_3, max_dist = 300, value = "Road_name", add_unique_value_name = TRUE)

DVC_names_added_4 <- merge(DVC_names_added_3, added_closest_named_DVC_3, all.x = T, by = c("x", "y")) %>% 
      mutate(Road_name = ifelse(is.na(Road_name) & !is.na(Road_name_new), Road_name_new, Road_name)) %>% 
      dplyr::select(-Road_name_new) 

nrow(DVC_names_added_4[is.na(DVC_names_added_4$Road_name), ])

# Export to check result in QGIS
# writeVector(DVC_names_added_4 %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE), "spatial_data/DVC_names_added_4.shp", overwrite = T)

# Make data set with reliable road traffic data (same roads, nearby traffic records)
traffic_vect_2010_2020 <- traffic_data_2010_2020 %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

traffic_data_edited_2010_2020 <- traffic_data_2010_2020 %>% 
      mutate(Road_name = gsub("\\(|\\)", "", Road_name))

DVC_names_added_4_vect <- DVC_names_added_4 %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

added_traffic <- get_nearest_value(DVC_names_added_4_vect, traffic_vect_2010_2020, value = c("All_motor_vehicles", "Road_name"))

DVC_traffic_added <- merge(DVC_names_added_4, added_traffic, all.x = T, by = c("x", "y", "Road_name"))

nrow(DVC_traffic_added[is.na(DVC_traffic_added$All_motor_vehicles), ])

# Check if there are records that have a road name, but no effort. In those cases, get the traffic data from the nearest record with matching
# road type.
DVC_named_no_traffic <- DVC_traffic_added %>% 
      filter(is.na(All_motor_vehicles) & !is.na(Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

traffic_data_simple <- traffic_data_edited_2010_2020 %>% 
      mutate(Road_name_simple = sub("^([[:alpha:]]*).*", "\\1", Road_name)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

added_simple_traffic <- get_nearest_value(DVC_named_no_traffic, traffic_data_simple, max_dist = 5000, value = c("All_motor_vehicles", "Road_name_simple"))

DVC_traffic_added_2 <- DVC_traffic_added %>% 
      mutate(Road_name_simple = sub("^([[:alpha:]]*).*", "\\1", Road_name)) %>% 
      merge(., added_simple_traffic, all.x = T, by = c("x", "y", "Road_name_simple")) %>% 
      mutate(All_motor_vehicles = ifelse(is.na(All_motor_vehicles.x) & !is.na(All_motor_vehicles.y), All_motor_vehicles.y, All_motor_vehicles.x)) %>% 
      dplyr::select(-All_motor_vehicles.x, -All_motor_vehicles.y)

nrow(DVC_traffic_added_2[is.na(DVC_traffic_added_2$All_motor_vehicles), ])

# In some cases matching didn't work, even though DVC points were very close to traffic data points. Merge these, if distance less than 300m
DVC_no_traffic <- DVC_traffic_added_2 %>% 
      filter(is.na(All_motor_vehicles)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)

added_closest_traffic <- get_nearest_value(DVC_no_traffic, traffic_data_simple, max_dist = 300, value = "All_motor_vehicles")

DVC_traffic_added_3 <- DVC_traffic_added_2 %>% 
      merge(., added_simple_traffic, all.x = T, by = c("x", "y")) %>% 
      mutate(All_motor_vehicles = ifelse(is.na(All_motor_vehicles.x) & !is.na(All_motor_vehicles.y), All_motor_vehicles.y, All_motor_vehicles.x)) %>% 
      dplyr::select(-All_motor_vehicles.x, -All_motor_vehicles.y)

nrow(DVC_traffic_added_3[is.na(DVC_traffic_added_3$All_motor_vehicles), ])

DVC_final_2010_2020 <- DVC_traffic_added_3 %>% 
      filter(!is.na(All_motor_vehicles))

# ------------------------------------------------------------------------------------------------------------------------
# Breeding Bird Survey Records - BTO
# ------------------------------------------------------------------------------------------------------------------------
BBS_1995_2022_final <- BBS %>% 
      rename_all(tolower) %>% 
      rename(species = english_name) %>% 
      filter(!is.na(count)) %>% 
      mutate(easting = rnrfa::osg_parse(square)$easting,
             northing = rnrfa::osg_parse(square)$northing,
             species = replace(species, species == "Roe Deer", "Capreolus capreolus"),
             species = replace(species, species == "Red Deer", "Cervus elaphus"),
             species = replace(species, species == "Fallow Deer", "Dama dama"),
             species = replace(species, species == "Reeves's Muntjac", "Muntiacus reevesi")) %>% 
      # Replace NAs with zero in the time columns, so we can calculate effort
      mutate_at(vars(matches("start|end")), ~replace_na(., 0)) %>% 
      mutate(zero_count = rowSums(select(., colnames(.)[grepl("start|end", colnames(.))]) == 0)) %>% 
      mutate(t1_early_duration = minute_difference(t1endearly, t1startearly),
             t2_early_duration = minute_difference(t2endearly, t2startearly),
             t1_late_duration = minute_difference(t1endlate, t1startlate),
             t2_late_duration = minute_difference(t2endlate, t2startlate),
             number_of_visits = 4 - zero_count/2) %>% 
      # Only keep records where single transect duration didn't exceed 3 hours
      filter_at(vars(t1_early_duration, t2_early_duration, t1_late_duration, t2_late_duration), ~ . <= 180) %>% 
      # Remove rows where no visit occurred
      filter(number_of_visits > 0) %>% 
      mutate(total_visit_duration = t1_early_duration + t2_early_duration + t1_late_duration + t2_late_duration,
             average_visit_duration = total_visit_duration / number_of_visits) %>% 
      filter(total_visit_duration > 0) %>% 
      group_by(square, easting, northing, species, total_visit_duration, number_of_visits, average_visit_duration) %>% 
      summarize(count = sum(count)) %>% 
      vect(geom = c("easting", "northing"), crs = crs("epsg:27700"), keepgeom = FALSE) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") %>% 
      mutate(presence = ifelse(count > 0, 1, 0))

# ------------------------------------------------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------------------------------------------------
save(inaturalist_final_2004_2020, BDS_final_2005_2016, DVC_final_2010_2020, BBS_1995_2022_final, file ="data-processed/deer_records.RData")








