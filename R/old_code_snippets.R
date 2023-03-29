# Morera-Pujol et al. 2022:
swf_100m <- rast('spatial_data_raw/SWF_ENG2015_100m.tif') %>% 
      terra::project(crs(km_projection))
# Morera-Pujol et al. 2022:
tcd_100m <- rast('spatial_data_raw/TCD2018_100m_Eng.tif') %>% 
      terra::project(crs(km_projection))


# Morera-Pujol et al. 2022:
slope_100m <- rast('spatial_data_raw/NASA_slope_Eng100m.tif') %>% 
      terra::project(crs(km_projection))


swf_1000m <- clamp(swf_100m, upper = 100, value = FALSE) %>% resample(., OS_grid_1km, method = 'bilinear') 
swf_1000m[is.nan(swf_1000m)] <- NA
# Indicator layer for swf data gaps
indicator_swf <- swf_1000m
values(indicator_swf)[!is.na(values(indicator_swf))] <- 1
values(indicator_swf)[is.na(values(indicator_swf))] <- 0
indicator_swf <- indicator_swf %>% mask(., ROI)
swf_1000m[is.na(swf_1000m)] <- 0
swf_1000m <- swf_1000m %>% mask(., ROI)


tcd_1000m <- clamp(tcd_100m, upper = 100, value = FALSE) %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'sum')
tcd_1000m[is.nan(tcd_1000m)] <- NA
# traffic_data <- data.table::fread("~/Downloads/dft_traffic_counts_raw_counts.csv")
# traffic_2020_2022 <- traffic_data %>%
#       filter(Year > 2019) %>%
#       dplyr::select(Road_name, Easting, Northing, All_motor_vehicles) %>%
#       group_by(Road_name, Easting, Northing) %>%
#       summarize(All_motor_vehicles = sum(All_motor_vehicles)) %>%
#       mutate(Easting = as.numeric(Easting), Northing = as.numeric(Northing)) %>%
#       vect(geom = c("Easting", "Northing"), crs = crs("epsg:27700")) %>%
#       terra::project(crs(km_projection)) %>%
#       terra::crop(ROI) %>%
#       as.data.frame(geom = "XY")
# traffic_2010_2020 <- traffic_data %>%
#       filter(Year >= 2010) %>%
#       dplyr::select(Road_name, Easting, Northing, All_motor_vehicles) %>%
#       group_by(Road_name, Easting, Northing) %>%
#       summarize(All_motor_vehicles = sum(All_motor_vehicles)) %>%
#       mutate(Easting = as.numeric(Easting), Northing = as.numeric(Northing)) %>%
#       vect(geom = c("Easting", "Northing"), crs = crs("epsg:27700")) %>%
#       terra::project(crs(km_projection)) %>%
#       terra::crop(ROI) %>%
#       as.data.frame(geom = "XY")
# write.csv(traffic_2020_2022, "data-raw/England_traffic_data_2020_2022.csv")
# write.csv(traffic_2010_2020, "data-raw/England_traffic_data_2010_2020.csv")

# ------------------------------------------------------------------------------------------------------------------------
# Deer collisions records - Highway Agency
# ------------------------------------------------------------------------------------------------------------------------
RTA_2020_2022 <- RTA %>% 
      dplyr::select(Road_Name, Easting, Northing) %>% 
      group_by(Road_Name, Easting, Northing) %>% 
      summarize(n_deer = n()) %>% 
      filter(Easting !=0) %>% 
      vect(geom = c("Easting", "Northing"), crs = crs("epsg:27700"), keepgeom = FALSE) %>% 
      terra::project(crs(km_projection)) %>% 
      terra::crop(ROI) %>% 
      as.data.frame(geom = "XY") %>% 
      mutate(index = row_number())

RTA_vect <- RTA_2020_2022 %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)
traffic_vect_2020_2022 <- traffic_data_2020_2022 %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)
plot(traffic_vect_2020_2022, col = "red")
plot(RTA_vect, add = T)

traffic_data_edited_2020_2022 <- traffic_data_2020_2022 %>% 
      rename(traffic_x = x, traffic_y = y) %>% 
      mutate(Road_name = gsub("\\(|\\)", "", Road_name))

# Data set with reliable road traffic data (same roads, nearby traffic records)
joined <- merge(RTA_2020_2022, traffic_data_edited_2020_2022, all.x = T, all.y = T, by.x = "Road_Name", by.y = "Road_name")
RTA_with_effort_2020_2022 <- joined %>% 
      group_by(index) %>% 
      mutate(dist_x = abs(x-traffic_x), dist_y = abs(y-traffic_y)) %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(dist_total = sum(dist_x, dist_y)) %>% 
      group_by(index) %>% 
      slice_min(order_by = dist_total) %>% 
      ungroup() %>% 
      dplyr::select(n_deer, All_motor_vehicles, x, y)


# iNaturalist effort:
# roads_distance <- rast('spatial_data_raw/distance_to_road.tif') %>% 
#       terra::project(crs(km_projection)) %>% 
#       mask(., ROI) 


# # DVC effort:
# traffic_data_2010_2020 <- read.csv("data-raw/England_traffic_data_2010_2020.csv") %>% 
#       vect(geom = c("x", "y"), crs =crs(km_projection), keepgeom = T)
# 
# # iNaturalist effort
# inaturalist_users <- vect("spatial_data/inaturalist_users.shp")


# all_roads <- rast('spatial_data_raw/allRoads_Eng60m.tif')
# all_roads_clamped <- all_roads %>% clamp(., lower = 0, upper = 0, value = FALSE)
# all_roads_rast <- subst(all_roads_clamped, 0, 1)
# all_roads_rast <- subst(all_roads_rast, NA, 0)
# writeRaster(all_roads_rast, filename = "spatial_data_raw/roads_sub_Eng60m.tif", overwrite = TRUE)

# roads_distance_1000m <- roads_distance %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'average')
# roads_distance_1000m[is.nan(roads_distance_1000m)] <- NA

# traffic_2010_2020_1000m <- rasterize(traffic_data_2010_2020, OS_grid_1km, field = "All_motor_vehicles", fun=sum) %>% mask(., ROI) 
# traffic_2010_2020_1000m[is.nan(traffic_2010_2020_1000m)] <- NA

# rast_10km <- OS_grid_1km
# res(rast_10km) <- 10
# inat_users_10km <- rasterize(inaturalist_users, rast_10km, fun=function(i){length(i)}, background=0) %>% mask(., ROI) 
# inat_users_10000m <- inat_users_10km %>% resample(., OS_grid_1km, method = 'max')

# ------------------------------------------------------------------------------------------------------------------------
# Export raster layers
# ------------------------------------------------------------------------------------------------------------------------
# writeRaster(swf_1000m, filename = "spatial_data/SmallWoodyFeatures_2015_Eng1000m.tif", overwrite=TRUE)
# writeRaster(tcd_1000m, filename = "spatial_data/TreeCoverDensity_2018_Eng1000m.tif", overwrite=TRUE)
# writeRaster(forest_patchiness_1000m, filename = "spatial_data/forest_patchiness_2018_Eng1000m.tif", overwrite=TRUE)
# writeRaster(elevation_1000m, filename = "spatial_data/NASA_elevation_Eng1000m.tif", overwrite=TRUE)
# writeRaster(slope_1000m, filename = "spatial_data/NASA_slope_Eng1000m.tif", overwrite=TRUE)
# writeRaster(pop_density, filename = "spatial_data/Eng_res_pop_1km_2011.tif", overwrite=TRUE)
# writeRaster(ghm_1000m, filename = "spatial_data/HumanModification_Eng1000m.tif", overwrite=TRUE)
# writeRaster(soil_moisture, filename = "spatial_data/soil_moisture_Eng1000m.tif", overwrite=TRUE)
# writeRaster(CEH_lcm_perc_aggr_1000m, filename = "spatial_data/CEH_lcm_perc_aggr_1000m.tif", overwrite=TRUE)
# writeRaster(CEH_perc_grassland_1000m, filename = "spatial_data/CEH_perc_grassland_1000m.tif", overwrite=TRUE)
# writeRaster(ancient_woodland_1000m, filename = "spatial_data/ancient_woodland_1000m.tif", overwrite=TRUE)
# writeRaster(lwf_1000m, filename = "spatial_data/lwf_1000m.tif", overwrite=TRUE)
# writeRaster(roads_distance_1000m, filename = "spatial_data/roads_distance_1000m.tif", overwrite=TRUE)
# writeRaster(traffic_2010_2020_1000m, filename = "spatial_data/traffic_2010_2020_1000m.tif", overwrite=TRUE)
# writeRaster(inat_users_10000m, filename = "spatial_data/inat_users_10000m.tif", overwrite=TRUE)



# OS_grid <- vect('spatial_data_raw/OSGB_Grid_1km.shp') %>% 
#       terra::project(crs(km_projection)) %>% 
#       mask(., ROI_vect) %>% 
#       rast(., resolution = res_predictions*3) 
#       
#       
linear_woody_features_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                               linear_woody_features_1.natural_splines = linear_woody_features_1.natural_splines,
                                                               linear_woody_features_2.natural_splines = linear_woody_features_2.natural_splines,
                                                               linear_woody_features_3.natural_splines = linear_woody_features_3.natural_splines); names(linear_woody_features_natural_splines_lc) <- paste0("linear_woody_features_natural_splines_lc", 1:100)
elevation_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                   elevation_1.natural_splines = elevation_1.natural_splines,
                                                   elevation_2.natural_splines = elevation_2.natural_splines,
                                                   elevation_3.natural_splines = elevation_3.natural_splines); names(elevation_natural_splines_lc) <- paste0("elevation_natural_splines_lc", 1:100)
grassland_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                   grassland_1.natural_splines = grassland_1.natural_splines,
                                                   grassland_2.natural_splines = grassland_2.natural_splines,
                                                   grassland_3.natural_splines = grassland_3.natural_splines); names(grassland_natural_splines_lc) <- paste0("grassland_natural_splines_lc", 1:100)
tree_cover_density_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                            tree_cover_density_1.natural_splines = tree_cover_density_1.natural_splines,
                                                            tree_cover_density_2.natural_splines = tree_cover_density_2.natural_splines,
                                                            tree_cover_density_3.natural_splines = tree_cover_density_3.natural_splines); names(tree_cover_density_natural_splines_lc) <- paste0("tree_cover_density_natural_splines_lc", 1:100)
soil_moisture_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                       soil_moisture_1.natural_splines = soil_moisture_1.natural_splines,
                                                       soil_moisture_2.natural_splines = soil_moisture_2.natural_splines,
                                                       soil_moisture_3.natural_splines = soil_moisture_3.natural_splines); names(soil_moisture_natural_splines_lc) <- paste0("soil_moisture_natural_splines_lc", 1:100)
broadleaf_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                   broadleaf_1.natural_splines = broadleaf_1.natural_splines,
                                                   broadleaf_2.natural_splines = broadleaf_2.natural_splines,
                                                   broadleaf_3.natural_splines = broadleaf_3.natural_splines); names(broadleaf_natural_splines_lc) <- paste0("broadleaf_natural_splines_lc", 1:100)
coniferous_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                    coniferous_1.natural_splines = coniferous_1.natural_splines,
                                                    coniferous_2.natural_splines = coniferous_2.natural_splines,
                                                    coniferous_3.natural_splines = coniferous_3.natural_splines); names(coniferous_natural_splines_lc) <- paste0("coniferous_natural_splines_lc", 1:100)
residential_population_log_natural_splines_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                                    residential_population_log_1.natural_splines = residential_population_log_1.natural_splines,
                                                                    residential_population_log_2.natural_splines = residential_population_log_2.natural_splines,
                                                                    residential_population_log_3.natural_splines = residential_population_log_3.natural_splines); names(residential_population_log_natural_splines_lc) <- paste0("residential_population_log_natural_splines_lc", 1:100)
all_lc_natural_splines <- c(linear_woody_features_natural_splines_lc, elevation_natural_splines_lc, grassland_natural_splines_lc, tree_cover_density_natural_splines_lc, soil_moisture_natural_splines_lc,
                            broadleaf_natural_splines_lc, coniferous_natural_splines_lc, residential_population_log_natural_splines_lc)


final_df <- prepareSplineData(final_df, final_df$grassland, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$coniferous, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$ancient_woodland, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$broadleaf, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$elevation, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$linear_woody_features, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$soil_moisture, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$arable, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$residential_population_log, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$tree_cover_density, method = "natural_splines", nsplines = 3)
final_df <- prepareSplineData(final_df, final_df$human_modification, method = "natural_splines", nsplines = 3)








# ------------------------------------------------------------------------------------------------------------------------
# Plot model output - smoothed spatial estimates
# ------------------------------------------------------------------------------------------------------------------------

inat <- predict(model, data, predictor = TRUE, fun = 'linear', n.samples = 1000, formula = ~ inat_filtered_biasField)
inat_raster <- rast(inat[["predictions"]])
DVC <- predict(model, data, predictor = TRUE, fun = 'linear', n.samples = 1000, formula = ~ DVC_filtered_biasField)
DVC_raster <- rast(DVC[["predictions"]])
shared <- predict(model, data, predictor = TRUE, fun = 'linear', n.samples = 1000, formula = ~ shared_spatial)
shared_raster <- rast(shared[["predictions"]])

species_scientific <- c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi")
species_common <- c("Red deer", "Roe deer", "Fallow deer", "Reeves's Muntjac")
spec_df <- data.frame(species_scientific, species_common)
spec_name <- spec_df$species_common[spec_df$species_scientific == species_choice]

pred_plot_median <- gplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = prediction_raster$adjusted_abundance_median) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("log(Estimated relative abundance)") + 
      scale_fill_viridis(name = "Median", na.value="transparent",
                         option = "H"); pred_plot_median

pred_plot_sd <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = prediction_raster$sd) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("log(Estimated standard error)") + 
      scale_fill_viridis(name = "Median", na.value="transparent",
                         option = "H"); pred_plot_sd










plot(log(extraction$density_1km2), extraction$median)



pop_size <- round(sum(values(prediction_raster$adjusted_abundance_median), na.rm = TRUE)*4)



# save.image(file = paste0('model-output/', model_name))


pred_plot_median <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = prediction_raster$adjusted_abundance_median) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle(paste0("Scaled abundance, population: ", pop_size)) + 
      scale_fill_viridis(name = "n/km\u00B2", na.value="transparent",
                         option = "H"); pred_plot_median

pred_plot_sd <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = prediction_raster$sd) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("log(Estimated standard error)") + 
      scale_fill_viridis(name = "Median", na.value="transparent",
                         option = "H")

corrplot <- ggplot(extraction, aes(x = median, y = log(density_1km2))) +
      geom_point() +
      geom_abline(intercept = coef(lm_fit)[1], slope = coef(lm_fit)[2], 
                  col = "blue", alpha = 0.5, linewidth = 1.2) +
      theme_bw() +
      ggtitle('Correlation: Estimate - observed') +
      xlab("log(Estimated relative abundance)") +
      ylab("log(Culling counts)") +
      annotate("text", x = min(extraction$median) + 0.2,
               y = max(log(extraction$density_1km2)) - 0.5,
               label = paste0("Intercept = ", round(intercept, 2), 
                              "\nSlope = ", round(slope, 2),
                              "\nR\u00B2 = ", r2), 
               size = 3, hjust = 0); corrplot

pred_inat_plot_median <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = inat_raster$median) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("iNaturalist spatial field") + 
      scale_fill_viridis(name = "Median", na.value="transparent")

pred_DVC_plot_median <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = DVC_raster$median) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("DVC spatial field") + 
      scale_fill_viridis(name = "Median", na.value="transparent")

pred_shared_plot_median <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = shared_raster$median) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Shared spatial field") + 
      scale_fill_viridis(name = "Median", na.value="transparent")



pred_plot <- egg::ggarrange(pred_plot_median, pred_plot_sd, corrplot, pred_shared_plot_median, pred_inat_plot_median, pred_DVC_plot_median, nrow = 2)
pred_plot_annotated <- ggpubr::annotate_figure(pred_plot, top = ggpubr::text_grob(spec_name, face = "bold", size = 14))

pdf(paste0("maps_", sub(".RData", ".pdf", model_name)), width = 18, height = 9)
pred_plot_annotated
dev.off()
# ggsave(filename = paste0("figures/", sub(".RData", ".png", model_name)), plot = pred_plot)
model_fixed$ID <- factor(model_fixed$ID, levels = c("broadleaf", "n_years_surveyed", "residential_population", "linear_woody_features",  "slope", "grassland", "coniferous",            
                                                    "ancient_woodland", "soil_moisture", "elevation",      
                                                    "total_visit_duration"))

# ------------------------------------------------------------------------------------------------------------------------
# Plot model output - effect plots
# ------------------------------------------------------------------------------------------------------------------------
# Forest plot, only useful for linear terms
model_fixed <- data.frame(model$summary.fixed) %>% 
      mutate(ID = rownames(.))
model_fixed$significant <- ifelse((model_fixed$X0.025quant > 0 & model_fixed$X0.975quant > 0)|(model_fixed$X0.025quant < 0 & model_fixed$X0.975quant < 0), "yes", "no")
model_fixed <- model_fixed[!grepl("intercept", model_fixed$ID), ]

forest_plot <- ggplot() + 
      geom_point(data = model_fixed, aes(y = ID, x = X0.5quant, col = significant)) +
      geom_errorbar(data = model_fixed, aes(y = ID, xmin = X0.025quant, xmax = X0.975quant, col = significant), width = 0.1) +
      geom_vline(aes(xintercept = 0), lty = 2, alpha = .3) +
      ylab("Posterior estimates") +
      xlab(element_blank()) +
      theme_minimal() +
      ggtitle(spec_name) +
      scale_colour_manual(values = c("black", "#D55E00"), guide = "none"); forest_plot

# Effect plots, for regression spline fits
original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
      mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
             sequence = as.numeric(gsub("\\D", "", rownames(.))))
original_values$covariate <- gsub('all.seq.all.seq.all.seq.', '', original_values$covariate)
original_values$orig_values[original_values$covariate == "residential_population_log"] <- exp(original_values$orig_values[original_values$covariate == "residential_population_log"])

scale_params <- median(model$summary.lincomb.derived$`0.5quant`)

effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                           sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                           quant_05 = model$summary.lincomb.derived$`0.5quant`,
                           quant_0025 = model$summary.lincomb.derived$`0.025quant`,
                           quant_0975 = model$summary.lincomb.derived$`0.975quant`)
effect_combs$covariate <- gsub('_sortBase', '', effect_combs$covariate)

effect_combs_m <- merge(original_values, effect_combs)


effects_plot <- ggplot(effect_combs_m) +
      geom_line(aes(x = orig_values, y = quant_05)) +
      geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
      geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
      ggthemes::theme_few() +
      facet_wrap(~ covariate, scale = 'free_x') +
      ggtitle(spec_name) +
      xlab("Covariate value") +
      ylab("Relative abundance"); effects_plot


ancient_woodland_1000m <- ifel(is.na(ancient_woodland), 0, 1) %>% mask(., ROI) %>% terra::aggregate(., fact = 10, fun = function(vals) {
      sum(vals == 1, na.rm = TRUE) / length(vals)    # Proportion of ancient woodland
}) %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'max')
ancient_woodland_1000m[is.nan(ancient_woodland_1000m)] <- NA

# Additional:
ancient_woodland <- rast('spatial_data_raw/ancient_woodland.tif') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) 
ancient_woodland[is.nan(ancient_woodland)] <- NA


# Morera-Pujol et al. 2022:
slope_100m <- rast('spatial_data_raw/NASA_slope_Eng100m.tif') %>% 
      terra::project(crs(km_projection))


slope_1000m <- slope_100m %>% mask(., ROI) %>% resample(., OS_grid_1km, method = 'average')
slope_1000m[is.nan(slope_1000m)] <- NA

# Population density (as more explicit alternative to human modification layer)
pop_density <- rast('spatial_data_raw/UK_residential_population_2011_1_km.asc') %>% 
      terra::project(crs(km_projection)) %>% 
      crop(., ROI) %>% 
      mask(., ROI) 

# Population density: Many NA values where it should be zero. Replace with zero.
pop_density[is.nan(pop_density)] <- 0
pop_density <- pop_density %>% mask(., ROI)
names(pop_density) <- "residential_population"


linear_woody_features_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                        linear_woody_features_1.sortBase = linear_woody_features_1.sortBase,
                                                        linear_woody_features_2.sortBase = linear_woody_features_2.sortBase); names(linear_woody_features_sortBase_lc) <- paste0("linear_woody_features_sortBase_lc", 1:100)
elevation_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                            elevation_1.sortBase = elevation_1.sortBase,
                                            elevation_2.sortBase = elevation_2.sortBase,
                                            elevation_3.sortBase = elevation_3.sortBase); names(elevation_sortBase_lc) <- paste0("elevation_sortBase_lc", 1:100)

soil_moisture_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                soil_moisture_1.sortBase = soil_moisture_1.sortBase,
                                                soil_moisture_2.sortBase = soil_moisture_2.sortBase,
                                                soil_moisture_3.sortBase = soil_moisture_3.sortBase); names(soil_moisture_sortBase_lc) <- paste0("soil_moisture_sortBase_lc", 1:100)

arable_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                             arable_log_1.sortBase = arable_log_1.sortBase,
                                             arable_log_2.sortBase = arable_log_2.sortBase); names(arable_log_sortBase_lc) <- paste0("arable_log_sortBase_lc", 1:100)
grassland_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                            grassland_1.sortBase = grassland_1.sortBase,
                                            grassland_2.sortBase = grassland_2.sortBase); names(grassland_sortBase_lc) <- paste0("grassland_sortBase_lc", 1:100)
tree_cover_density_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                         tree_cover_density_log_1.sortBase = tree_cover_density_log_1.sortBase,
                                                         tree_cover_density_log_2.sortBase = tree_cover_density_log_2.sortBase); names(tree_cover_density_log_sortBase_lc) <- paste0("tree_cover_density_log_sortBase_lc", 1:100)
night_time_light_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                       night_time_light_log_1.sortBase = night_time_light_log_1.sortBase,
                                                       night_time_light_log_2.sortBase = night_time_light_log_2.sortBase); names(night_time_light_log_sortBase_lc) <- paste0("night_time_light_log_sortBase_lc", 1:100)
road_cover_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                                 road_cover_log_1.sortBase = road_cover_log_1.sortBase,
                                                 road_cover_log_2.sortBase = road_cover_log_2.sortBase); names(road_cover_log_sortBase_lc) <- paste0("road_cover_log_sortBase_lc", 1:100)
builtup_log_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                              builtup_log_1.sortBase = builtup_log_1.sortBase,
                                              builtup_log_2.sortBase = builtup_log_2.sortBase); names(builtup_log_sortBase_lc) <- paste0("builtup_log_sortBase_lc", 1:100)
builtup_sortBase_lc <- inla.make.lincombs(BDS_filtered_intercept = rep(1, 100),
                                          builtup_1.sortBase = builtup_1.sortBase,
                                          builtup_2.sortBase = builtup_2.sortBase,
                                          builtup_3.sortBase = builtup_3.sortBase); names(builtup_sortBase_lc) <- paste0("builtup_sortBase_lc", 1:100)

all_lc_sortBase <- c(linear_woody_features_sortBase_lc, elevation_log_sortBase_lc, arable_log_sortBase_lc, 
                     grassland_sortBase_lc, tree_cover_density_log_sortBase_lc, night_time_light_log_sortBase_lc,
                     road_cover_log_sortBase_lc, builtup_log_sortBase_lc)
covariates_used <- c('linear_woody_features_1.sortBase', 'linear_woody_features_2.sortBase',
                     'elevation_log_1.sortBase', 'elevation_log_2.sortBase',
                     'arable_log_1.sortBase', 'arable_log_2.sortBase',
                     'tree_cover_density_log_1.sortBase', 'tree_cover_density_log_2.sortBase',
                     'night_time_light_log_1.sortBase', 'night_time_light_log_2.sortBase',
                     'road_cover_log_1.sortBase', 'road_cover_log_2.sortBase',
                     'x_coord', 'y_coord')

# Forest plot, only useful for linear terms
model_fixed <- data.frame(model$summary.fixed) %>% 
      mutate(ID = rownames(.))
model_fixed$significant <- ifelse((model_fixed$X0.025quant > 0 & model_fixed$X0.975quant > 0)|(model_fixed$X0.025quant < 0 & model_fixed$X0.975quant < 0), "yes", "no")
model_fixed <- model_fixed[!grepl("intercept", model_fixed$ID), ]

forest_plot <- ggplot() + 
      geom_point(data = model_fixed, aes(y = ID, x = X0.5quant, col = significant)) +
      geom_errorbar(data = model_fixed, aes(y = ID, xmin = X0.025quant, xmax = X0.975quant, col = significant), width = 0.1) +
      geom_vline(aes(xintercept = 0), lty = 2, alpha = .3) +
      ylab("Posterior estimates") +
      xlab(element_blank()) +
      theme_minimal() +
      ggtitle(species_choice) +
      scale_colour_manual(values = c("black", "#D55E00"), guide = "none"); forest_plot