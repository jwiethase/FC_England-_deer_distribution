rm(list = ls())
library(PointedSDMs)
library(INLA)
library(tidyverse)
library(terra)
library(rgdal)
library(inlabru)
sapply(list.files(path = 'source', full.names = T), source, .GlobalEnv)

# ------------------------------------------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------------------------------------------
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

# Species, one of: "Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi"
species_choice <- "Cervus elaphus"

# First guess of spatial range was 1/3 of the study area
# Full N-S study area extent is 645km, estimated range should be smaller than that
# estimated_range = (ext(ROI)[4] - ext(ROI)[3])/3
# max.edge = estimated_range/12
covariates_used <- c('elevation_log_1.sortBase', 'elevation_log_2.sortBase',
                     'linear_woody_features_1.sortBase', 'linear_woody_features_2.sortBase',
                     'grassland_log_1.sortBase', 'grassland_log_2.sortBase',
                     'tree_cover_density_log_1.sortBase', 'tree_cover_density_log_2.sortBase',
                     'builtup_log_1.sortBase', 'builtup_log_2.sortBase',
                     'night_time_light_log_1.sortBase', 'night_time_light_log_2.sortBase')
# covariates_used <- c('arable_log_1.sortBase', 'arable_log_2.sortBase', 
#                      'tree_cover_density_log_1.sortBase', 'tree_cover_density_log_2.sortBase',
#                      'builtup_log_1.sortBase', 'builtup_log_2.sortBase', 
#                      'elevation_1.sortBase', 'elevation_2.sortBase')

# What makes sense?
# Cover of builtup areas is great umbrella: Correlates with light at night, residential population, road cover
# Tree cover density: Very likely important, and correlates well with broadleaf forest.
# Not so much coniferous forest, but there isn't much in England to begin with
# Arable: Very specific type of habitat, likely to influence deer presence. Grassland is also very 
# important, but correlates strongly with arable - less grass means more arable
# Soil moisture: Correlates strongly with elevation and human modification, but contains
# more variation than elevation for muntiac, and is hence more informative

# Try elevation, try 2 splines

offset_vars <- c('n_users_100km2', 'All_motor_vehicles')
# "broadleaf_1.sortBase", "broadleaf_2.sortBase",
# "coniferous_1.sortBase", "coniferous_2.sortBase",

# ------------------------------------------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------------------------------------------
covariates <- raster::stack(readRDS('data-processed/raster_covariates.RDS'))
covariates_sp <- as(covariates, 'SpatialPixelsDataFrame')
load('data-processed/deer_records.RData')

# Study area
ROI <- rgdal::readOGR('spatial_data/england_outline_27700_simple.shp') %>% 
      spTransform(., km_projection)

# ------------------------------------------------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------------------------------------------------
BDS_filtered <- BDS_final_2005_2016 %>%
      filter(species == species_choice) %>% 
      as.data.frame(.)

inat_filtered <- inaturalist_final_2004_2020 %>%
      filter(species == species_choice) %>% 
      as.data.frame(.) %>% 
      mutate(distance_to_road = scale(log(distance_to_road+1)),
             n_users_100km2 = scale(log(n_users_100km2)),
             presence = NULL)

DVC_filtered <- DVC_final_2010_2020 %>%
      filter(species == species_choice) %>% 
      as.data.frame(.) %>% 
      mutate(All_motor_vehicles = scale(log(All_motor_vehicles)),
             presence = NULL)

BBS_filtered <- BBS_1995_2022_final %>%
      filter(species == species_choice) %>% 
      as.data.frame(.) %>% 
      mutate(total_visit_duration = scale(log(total_visit_duration)),
             presence = ifelse(count > 0, 1, 0))

NROW(BDS_filtered %>% filter(presence == 1)) + NROW(inat_filtered) + NROW(DVC_filtered) + NROW(BBS_filtered %>% filter(presence == 1))

covariates_selection <- covariates_sp[, names(covariates_sp) %in% covariates_used]
names(covariates_selection)
covariates_selection_df <- data.frame(covariates_selection) %>% dplyr::select(-x, -y)
usdm::vif(covariates_selection_df)

# Only scale covariates that have not been split into splines (those have already been scaled)
# covariates_selection@data <- data.frame(scale(covariates_selection@data[!names(covariates_selection) %in% grep("\\.s", names(covariates_selection), value = TRUE)]), 
                                        # covariates_selection@data[names(covariates_selection) %in% grep("\\.s", names(covariates_selection), value = TRUE)])

# ------------------------------------------------------------------------------------------------------------------------
# Mesh
# ------------------------------------------------------------------------------------------------------------------------
# max.edge = diff(range(st_coordinates(intensity_sf)[, 1]))/(3*5)
all_obs_x <- c(BDS_filtered$x[BDS_filtered$presence == 1], inat_filtered$x, DVC_filtered$x, BBS_filtered$x[BBS_filtered$presence == 1])
all_obs_y <- c(BDS_filtered$y[BDS_filtered$presence == 1], inat_filtered$y, DVC_filtered$y, BBS_filtered$y[BBS_filtered$presence == 1])

max.edge = round(diff(range(all_obs_x))/19, digits = 1)
max.edge = 23
region.bdry <- inla.sp2segment(ROI)
mesh <- inla.mesh.2d(boundary = region.bdry, 
                     cutoff = max.edge/2, 
                     max.edge = c(max.edge, max.edge*4), 
                     offset = c(max.edge, max.edge*5))
mesh$crs <- km_projection

# Ref:
# "If you want a computational mesh, you have to modify the mesh after running inference. This is not an inconsistency in the 
# modeling approach, since it is only a computational issue. I recommend running inference at first with the settings above, 
# then running inference, and looking at inla.result$summary.hyperpar$mean and finding the posterior mean estimate for the range. 
# Then, go back and redo max.edge to be between 1/5 and 1/10 of this range. You can have a smaller max.edge but this is just 
# an unnecessary “waste of time”. Similarly, ensure the outer extension is close to the estimate for the range."

# ------------------------------------------------------------------------------------------------------------------------
# Priors
# ------------------------------------------------------------------------------------------------------------------------
prior.range = c(round(0.2*diff(range(all_obs_y)), digits = 1), 0.5)
prior.sigma = c(round(0.2*diff(range(all_obs_y))/25, digits = 1), 0.1)
# prior.range = c(256.7, 0.5)
# prior.sigma = c(10.3, 0.1)
inat.prior.range = c(round(0.2*diff(range(inat_filtered$y)), digits = 1), 0.5)
inat.prior.sigma = c(round(0.2*diff(range(inat_filtered$y))/25, digits = 1), 0.1)
# inat.prior.range = c(256.7, 0.5)
# inat.prior.sigma = c(10.3, 0.1)
DVC.prior.range = c(round(0.2*diff(range(DVC_filtered$y)), digits = 1), 0.5)
DVC.prior.sigma = c(round(0.2*diff(range(DVC_filtered$y))/25, digits = 1), 0.1)
# DVC.prior.range = c(256.7, 0.5)
# DVC.prior.sigma = c(10.3, 0.1)

# ------------------------------------------------------------------------------------------------------------------------
# Prepare model
# ------------------------------------------------------------------------------------------------------------------------
model_name <- paste0(sub(" ", "_", species_choice), "_E", max.edge, "_R", prior.range[1], "_", prior.range[2], "_S", prior.sigma[1], "_", prior.sigma[2], ".RData")

model_setup <- intModel(BDS_filtered, inat_filtered, DVC_filtered, BBS_filtered,
                        Mesh = mesh, Projection = CRS(km_projection), responseCounts = 'count', spatialCovariates = covariates_selection,
                        responsePA = 'presence', Coordinates = c('x', 'y'), 
                        pointCovariates = 'total_visit_duration', Offset = offset_vars)

model_setup$addBias('DVC_filtered')
model_setup$addBias('inat_filtered')

model_setup$specifySpatial(sharedSpatial = TRUE, prior.range = prior.range,
                           prior.sigma = prior.sigma)
model_setup$specifySpatial(Bias = 'inat_filtered', prior.range = inat.prior.range, prior.sigma = inat.prior.sigma)
model_setup$specifySpatial(Bias = 'DVC_filtered', prior.range = DVC.prior.range, prior.sigma = DVC.prior.sigma)

save.image(paste0("model-files/", model_name))

