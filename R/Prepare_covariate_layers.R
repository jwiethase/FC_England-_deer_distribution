rm(list = ls())
library(tidyverse)
library(terra)
library(corrplot)
library(INLA)
library(usdm)
source("source/deerSDM_misc_functions.R")
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

# ------------------------------------------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------------------------------------------
ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))
# Morera-Pujol et al. 2022:
swf_100m <- rast('spatial_data_raw/SWF_ENG2015_100m.tif') %>% 
      terra::project(crs(km_projection))
# Morera-Pujol et al. 2022:
tcd_100m <- rast('spatial_data_raw/TCD2018_100m_Eng.tif') %>% 
      terra::project(crs(km_projection))
# Morera-Pujol et al. 2022:
elevation_100m <- rast('spatial_data_raw/NASA_dem_Eng100m.tif') %>% 
      terra::project(crs(km_projection))

human_modification <- rast('spatial_data_raw/human_modification_2016_Eng1000m.tif') %>% 
      terra::project(crs(km_projection))
# Gill & Morgan 2010:
soil_moisture <- rast('spatial_data_raw/EIDC_soil_moisture_Eng1000m.tif') %>% 
      terra::project(crs(km_projection))
# Additional:
CEH_lcm_perc_aggr_1000m <- rast('spatial_data_raw/gb2021lcm1km_percentage_aggregate.tif') %>% 
      terra::project(crs(km_projection))
# Additional:
CEH_dom_aggr_1000m <- rast('spatial_data_raw/gb2021lcm1km_dominant_target.tif') %>% 
      terra::project(crs(km_projection))

linear_woody_features <- rast('spatial_data_raw/Eng_LWF_30m_2016.tif') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) 
# Ciach & Frohlich 2019:
NTL_VIIRS_500m <- rast('spatial_data_raw/NTL_Eng500m.tif') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) 
# Population density (as more explicit alternative to human modification layer)
pop_density <- rast('spatial_data_raw/UK_residential_population_2011_1_km.asc') %>% 
      terra::project(crs(km_projection)) %>% 
      crop(., ROI) %>% 
      mask(., ROI) 

roads_Eng60m <- rast('spatial_data_raw/allRoads_Eng60m.tif') %>% 
      terra::project(crs(km_projection)) %>% 
      crop(., ROI) %>% 
      mask(., ROI) 

# ------------------------------------------------------------------------------------------------------------------------
# Reduce raster resolution
# ------------------------------------------------------------------------------------------------------------------------
# Make dummy raster to project to. Use OS 1km grid.
OS_grid_1km <- vect('spatial_data_raw/OSGB_Grid_1km.shp') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) %>% 
      rast(., resolution = 1) 

# Summarise linear woody features. Length of feature appears to be problematic, as sums of lengths are noticeably 
# different for rectangular area in North. Use area proportion of linear woody features instead, based on
# 30m resolution linear features.
lwf_1000m <- ifel(is.na(linear_woody_features), 0, 1) %>% mask(., ROI) %>% terra::aggregate(., fact = 30, fun = function(vals) {
      sum(vals == 1, na.rm = TRUE) / length(vals)}) %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'max')  
lwf_1000m[is.na(lwf_1000m)] <- 0
lwf_1000m <- lwf_1000m %>% mask(., ROI)
    
tcd_1000m <- clamp(tcd_100m, upper = 100, value = FALSE) %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'bilinear')
tcd_1000m[is.nan(tcd_1000m)] <- NA

#Hemami et al. 2005: Muntjac density was higher in forest blocks with a greater ratio of open habitat perimeter to forest area
forest_patchiness_1000m <- clamp(tcd_100m, upper = 100, value = FALSE) %>% mask(., ROI) %>% aggregate(., fact = 10, fun = 'sd') %>% resample(., OS_grid_1km, method = 'bilinear')
forest_patchiness_1000m[is.nan(forest_patchiness_1000m)] <- NA
names(forest_patchiness_1000m) <- "forest_patchiness_1000m"

elevation_1000m <- elevation_100m %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'average')
elevation_1000m[is.nan(elevation_1000m)] <- NA

ghm_1000m <- human_modification %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'bilinear')
ghm_1000m[is.nan(ghm_1000m)] <- NA

soil_moisture <- soil_moisture %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'bilinear')
soil_moisture[is.nan(soil_moisture)] <- NA

CEH_lcm_perc_aggr_1000m <- CEH_lcm_perc_aggr_1000m %>% crop(., ROI) %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'bilinear')
names(CEH_lcm_perc_aggr_1000m) <- c("perc_broadleaf", "perc_coniferous", "perc_arable", "perc_improved_grass", "perc_semi_nat_grass",
                                    "perc_mount_heath_bog", "perc_saltwater", "perc_freshwater", "perc_coastal", "perc_builtup_gardens")
CEH_lcm_perc_aggr_1000m[is.nan(CEH_lcm_perc_aggr_1000m)] <- NA

CEH_perc_grass_stack <- c(CEH_lcm_perc_aggr_1000m$perc_improved_grass, CEH_lcm_perc_aggr_1000m$perc_semi_nat_grass)

# Hemami et al. 2004, Zini et al. 2021
CEH_perc_grassland_1000m <- app(CEH_perc_grass_stack, sum)
names(CEH_perc_grassland_1000m) <- "all_grass"
CEH_perc_grassland_1000m[is.nan(CEH_perc_grassland_1000m)] <- NA

CEH_rural_1000m_sub <- terra::ifel(CEH_dom_aggr_1000m < 20, 1, 0) %>% 
      crop(., ROI) %>% mask(., ROI) 
CEH_suburban_1000m_sub <- terra::ifel(CEH_dom_aggr_1000m == 21, 1, 0) %>% 
      crop(., ROI) %>% mask(., ROI) 
CEH_urban_1000m_sub <- terra::ifel(CEH_dom_aggr_1000m == 20, 1, 0) %>% 
      crop(., ROI) %>% mask(., ROI) 

NTL_VIIRS_1000m <- NTL_VIIRS_500m %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'bilinear')
NTL_VIIRS_1000m[is.nan(NTL_VIIRS_1000m)] <- NA

# Population density: Many NA values where it should be zero. Replace with zero.
pop_density[is.nan(pop_density)] <- 0
pop_density <- pop_density %>% mask(., ROI)
names(pop_density) <- "residential_population"

# Road density
roads_1000m <- ifel(is.na(roads_Eng60m), 0, 1) %>% mask(., ROI) %>% terra::aggregate(., fact = 10, fun = function(vals) {
      sum(vals == 1, na.rm = TRUE) / length(vals)}) %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'bilinear')  
roads_1000m[is.na(roads_1000m)] <- 0
roads_1000m <- roads_1000m %>% mask(., ROI)
names(roads_1000m) <- "road_density"

# ------------------------------------------------------------------------------------------------------------------------
# Covariate selection
# ------------------------------------------------------------------------------------------------------------------------
# Check correlation between variables
full_stack <- c(lwf_1000m, tcd_1000m, forest_patchiness_1000m, elevation_1000m, 
                ghm_1000m, soil_moisture, CEH_lcm_perc_aggr_1000m, CEH_perc_grassland_1000m,
                NTL_VIIRS_1000m, pop_density, roads_1000m)
full_df <- data.frame(full_stack) %>% drop_na()

res <- cor(full_df)
png(filename = 'figures/corr_plot.png', width = 45, height = 30, units = "cm", res = 300)
corrplot(res, type = "upper", order = "hclust", 
         method = "color", addCoef.col="black", number.cex=0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)
dev.off()

# Make selection
final_stack <- c(lwf_1000m, CEH_perc_grassland_1000m, 
                 CEH_lcm_perc_aggr_1000m$perc_coniferous, 
                 elevation_1000m, CEH_lcm_perc_aggr_1000m$perc_broadleaf, CEH_lcm_perc_aggr_1000m$perc_arable,
                 soil_moisture, tcd_1000m, 
                 NTL_VIIRS_1000m, pop_density, CEH_lcm_perc_aggr_1000m$perc_builtup_gardens, roads_1000m,
                 CEH_lcm_perc_aggr_1000m$perc_mount_heath_bog)
names(final_stack) <- c('linear_woody_features', 'grassland', 
                        'coniferous', 'elevation', 'broadleaf', 'arable', 
                        'soil_moisture', 'tree_cover_density', 
                        'night_time_light', 'residential_population', 'builtup',
                        'road_cover', 'mount_heath_bog')

final_stack <- raster::stack(final_stack)

# Log-transform covariate where applicable
final_stack$elevation[final_stack$elevation < 0] <- 0
final_stack$elevation_log <- log(final_stack$elevation+1)
final_stack$soil_moisture_log <- log(final_stack$soil_moisture+1)
final_stack$grassland_log <- log(final_stack$grassland+1)
final_stack$linear_woody_features_log <- log(final_stack$linear_woody_features+1)
final_stack$arable_log <- log(final_stack$arable+1)
final_stack$tree_cover_density_log <- log(final_stack$tree_cover_density+1)
final_stack$night_time_light_log <- log(final_stack$night_time_light+1)
final_stack$residential_population_log <- log(final_stack$residential_population+1)
final_stack$builtup_log <- log(final_stack$builtup+1)
final_stack$road_cover_log <- log(final_stack$road_cover+1)
final_stack$mount_heath_bog_log <- log(final_stack$mount_heath_bog+1)
final_stack$coniferous_log <- log(final_stack$coniferous+1)

# ------------------------------------------------------------------------------------------------------------------------
# Make splines
# ------------------------------------------------------------------------------------------------------------------------
final_df <- as.data.frame(final_stack, xy = TRUE)  %>%
      filter(!if_all(-c(x, y), is.na))

final_df_cortest <- data.frame(final_df) %>% 
      drop_na() %>%
      dplyr::select(matches("log") | matches("^(x|y)$"))

res2 <- cor(final_df_cortest)
png(filename = 'figures/corr_plot_final.png', width = 45, height = 30, units = "cm", res = 300)
corrplot(res2, type = "upper", order = "hclust", 
         method = "color", addCoef.col="black", number.cex=0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)
dev.off()

# final_df <- read.csv("data-processed/final_covariates_df.csv") %>% dplyr::select(-X)
final_df <- prepareSplineData(final_df, final_df$elevation_log, method = "gmm", nsplines = 2) 
final_df <- prepareSplineData(final_df, final_df$elevation, method = "gmm", nsplines = 2) 
final_df <- prepareSplineData(final_df, final_df$linear_woody_features, method = "gmm", nsplines = 2) 
final_df <- prepareSplineData(final_df, final_df$arable, method = "gmm", nsplines = 2) 
final_df <- prepareSplineData(final_df, final_df$grassland, method = "gmm", nsplines = 2) 
final_df <- prepareSplineData(final_df, final_df$grassland_log, method = "gmm", nsplines = 2) 
final_df <- prepareSplineData(final_df, final_df$tree_cover_density_log, method = "gmm", nsplines = 2)
final_df <- prepareSplineData(final_df, final_df$tree_cover_density, method = "gmm", nsplines = 2)
final_df <- prepareSplineData(final_df, final_df$builtup, method = "gmm", nsplines = 2)
final_df <- prepareSplineData(final_df, final_df$builtup_log, method = "gmm", nsplines = 2)
final_df <- prepareSplineData(final_df, final_df$night_time_light_log, method = "gmm", nsplines = 2)
final_df <- prepareSplineData(final_df, final_df$night_time_light, method = "gmm", nsplines = 2)
final_df <- prepareSplineData(final_df, final_df$residential_population_log, method = "quantile", nsplines = 2, user_cp_quantiles = c(0.01, 0.105))

# Combine the unscaled prediction vectors
rm(all.seq)
all.seq <- mget(ls(pattern = "*.seq"))
# all.seq <- ls(arable_log.seq, tree_cover_density_log.seq, builtup_log.seq, soil_moisture.seq)
covariates <- rast(final_df, type = 'xyz', crs = crs(km_projection))

# Make linear combinations for effect plots
arable_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                             arable_1.sortBase = arable_1.sortBase,
                                             arable_2.sortBase = arable_2.sortBase); names(arable_sortBase_lc) <- paste0("arable_sortBase_lc", 1:100)
grassland_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                         grassland_log_1.sortBase = grassland_log_1.sortBase,
                                         grassland_log_2.sortBase = grassland_log_2.sortBase); names(grassland_log_sortBase_lc) <- paste0("grassland_log_sortBase_lc", 1:100)
tree_cover_density_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                     tree_cover_density_log_1.sortBase = tree_cover_density_log_1.sortBase,
                                                     tree_cover_density_log_2.sortBase = tree_cover_density_log_2.sortBase); names(tree_cover_density_log_sortBase_lc) <- paste0("tree_cover_density_log_sortBase_lc", 1:100)
tree_cover_density_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                         tree_cover_density_1.sortBase = tree_cover_density_1.sortBase,
                                                         tree_cover_density_2.sortBase = tree_cover_density_2.sortBase); names(tree_cover_density_sortBase_lc) <- paste0("tree_cover_density_sortBase_lc", 1:100)
builtup_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                          builtup_log_1.sortBase = builtup_log_1.sortBase,
                                          builtup_log_2.sortBase = builtup_log_2.sortBase); names(builtup_log_sortBase_lc) <- paste0("builtup_log_sortBase_lc", 1:100)
builtup_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                              builtup_1.sortBase = builtup_1.sortBase,
                                              builtup_2.sortBase = builtup_2.sortBase); names(builtup_sortBase_lc) <- paste0("builtup_sortBase_lc", 1:100)
linear_woody_features_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                linear_woody_features_1.sortBase = linear_woody_features_1.sortBase,
                                                linear_woody_features_2.sortBase = linear_woody_features_2.sortBase); names(linear_woody_features_sortBase_lc) <- paste0("linear_woody_features_sortBase_lc", 1:100)
elevation_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                elevation_1.sortBase = elevation_1.sortBase,
                                                elevation_2.sortBase = elevation_2.sortBase); names(elevation_sortBase_lc) <- paste0("elevation_sortBase_lc", 1:100)
elevation_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                            elevation_log_1.sortBase = elevation_log_1.sortBase,
                                            elevation_log_2.sortBase = elevation_log_2.sortBase); names(elevation_log_sortBase_lc) <- paste0("elevation_log_sortBase_lc", 1:100)
night_time_light_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                            night_time_light_log_1.sortBase = night_time_light_log_1.sortBase,
                                            night_time_light_log_2.sortBase = night_time_light_log_2.sortBase); names(night_time_light_log_sortBase_lc) <- paste0("night_time_light_log_sortBase_lc", 1:100)
all_lc <- c(grassland_log_sortBase_lc, arable_sortBase_lc, tree_cover_density_log_sortBase_lc, tree_cover_density_sortBase_lc, builtup_log_sortBase_lc, builtup_sortBase_lc, 
            elevation_sortBase_lc, elevation_log_sortBase_lc, linear_woody_features_sortBase_lc, night_time_light_log_sortBase_lc)

pdf(file = "figures/covars_plot_1.pdf")
plot(covariates[[1:16]])
dev.off()

pdf(file = "figures/covars_plot_2.pdf")
plot(covariates[[17:32]])
dev.off()

pdf(file = "figures/covars_plot_3.pdf")
plot(covariates[[33:length(names(covariates))]])
dev.off()

saveRDS(object = covariates, file = 'data-processed/raster_covariates.RDS')
save(all.seq, all_lc, file = 'data-processed/lincomb_files.RData')




