args = commandArgs(trailingOnly = TRUE)
library(PointedSDMs)
library(inlabru)
suppressPackageStartupMessages(library(INLA, quietly=TRUE))
library(terra)
library(tidyterra)
library(tidyverse)
library(egg)
library(ggpubr)
library(rgeos)
library(rgdal)
library(gridExtra)
library(viridis)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
setwd("/users/jhw538/scratch/FC_deer/data")

i <- as.numeric(args[1])

covariates_used <- c('residential_population_log', 'arable', 'linear_woody_features_1.s', 'linear_woody_features_2.s', 'grassland_1.s', 'grassland_2.s', 'coniferous_log_1.s', 'coniferous_log_2.s', 'ancient_woodland_1.s', 'ancient_woodland_2.s', 'soil_moisture_log_1.s', 'soil_moisture_log_2.s', 'elevation_log_1.s', 'elevation_log_2.s', 'broadleaf_log_1.s', 'broadleaf_log_2.s')

#covariates_used <- c('residential_population_log', 'arable', 'linear_woody_features', 'grassland', 'coniferous', 'ancient_woodland', 'soil_moisture', 'elevation', 'broadleaf')

offset_vars <- c('distance_to_road', 'All_motor_vehicles')

multiplier <- c(0.2, 0.3, 0.4)
edges <- seq(12, 28, by = 1)
specs <- c("Capreolus capreolus", "Cervus elaphus", "Muntiacus reevesi", "Dama dama")

mult_combs <- crossing(multiplier, edges, specs)
comb_values <- mult_combs[i, ]
print(comb_values)
multiplier <- comb_values$multiplier

max.edge <- comb_values$edges
species_choice <- comb_values$specs

res_predictions <- 2

# ------------------------------------------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------------------------------------------
covariates <- raster::stack(readRDS('raster_covariates.RDS'))
covariates_sp <- as(covariates, 'SpatialPixelsDataFrame')
print(names(covariates_sp))
load('deer_records.RData')

validation_cull <- read.csv('cull_final_2017_2020.csv') %>% 
      filter(species == species_choice) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE) 

# Study area
ROI <- rgdal::readOGR('england_outline_27700_simple.shp') %>% 
      spTransform(., km_projection)

ROI_vect <- vect(ROI) %>% 
      terra::project(crs(km_projection))

OS_grid_100m <- vect('OSGB_Grid_1km.shp') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI_vect) %>% 
      rast(., resolution = 0.1) 

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
             n_users_25km2 = scale(log(n_users_25km2)),
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

covariates_selection <- covariates_sp[, names(covariates_sp) %in% covariates_used]
covariates_selection@data <- data.frame(scale(covariates_selection@data[!names(covariates_selection) %in% grep("\\.s", names(covariates_selection), value = TRUE)]), covariates_selection@data[names(covariates_selection) %in% grep("\\.s", names(covariates_selection), value = TRUE)])
print("Data prep done")

# ------------------------------------------------------------------------------------------------------------------------
# Mesh
# ------------------------------------------------------------------------------------------------------------------------
all_obs_x <- c(BDS_filtered$x[BDS_filtered$presence == 1], inat_filtered$x, DVC_filtered$x, BBS_filtered$x[BBS_filtered$presence == 1])
all_obs_y <- c(BDS_filtered$y[BDS_filtered$presence == 1], inat_filtered$y, DVC_filtered$y, BBS_filtered$y[BBS_filtered$presence == 1])

region.bdry <- inla.sp2segment(ROI)
mesh <- inla.mesh.2d(loc = cbind(all_obs_x, all_obs_y),
                     boundary = region.bdry, 
                     cutoff = max.edge/2, 
                     max.edge = c(max.edge, max.edge*4), 
                     offset = c(max.edge, max.edge*5))
mesh$crs <- km_projection

# ------------------------------------------------------------------------------------------------------------------------
# Priors
# ------------------------------------------------------------------------------------------------------------------------
prior.range = c(round(multiplier*diff(range(all_obs_y)), digits = 1), 0.5)
prior.sigma = c(round(multiplier*diff(range(all_obs_y))/25, digits = 1), 0.1)

inat.prior.range = c(round(multiplier*diff(range(inat_filtered$y)), digits = 1), 0.5)
inat.prior.sigma = c(round(multiplier*diff(range(inat_filtered$y))/25, digits = 1), 0.1)

DVC.prior.range = c(round(multiplier*diff(range(DVC_filtered$y)), digits = 1), 0.5)
DVC.prior.sigma = c(round(multiplier*diff(range(DVC_filtered$y))/25, digits = 1), 0.1)

print("# ------------------------------------------------------------------------------------------------------------------------")
print(paste0(c("iNaturalist range: ", "P(estimated R < Range): "), inat.prior.range)) 
print(paste0(c("iNaturalist sigma: ", "P(estimated sigma > sigma): "), inat.prior.sigma)) 
print(paste0(c("DVC range: ", "P(estimated R < Range): "), DVC.prior.range)) 
print(paste0(c("DVC sigma: ", "P(estimated sigma > sigma): "), DVC.prior.sigma)) 
print("# ------------------------------------------------------------------------------------------------------------------------")

# ------------------------------------------------------------------------------------------------------------------------
# Prepare model
# ------------------------------------------------------------------------------------------------------------------------
model_name <- paste0(sub(" ", "_", species_choice), "_E", max.edge, "_R", prior.range[1], "_", prior.range[2], "_S", prior.sigma[1], "_", prior.sigma[2], ".RData")
print(model_name)

print(species_choice)
print(head(BDS_filtered))
print(head(inat_filtered))
print(names(covariates_selection))
model_setup <- intModel(BDS_filtered, inat_filtered, DVC_filtered, BBS_filtered,
                        Mesh = mesh, Projection = CRS(km_projection), responseCounts = 'count', spatialCovariates = covariates_selection,
                        responsePA = 'presence', Coordinates = c('x', 'y'), 
                        pointCovariates = 'total_visit_duration', Offset = offset_vars)
print(model_setup)

model_setup$addBias('DVC_filtered')
model_setup$addBias('inat_filtered')
model_setup$specifySpatial(sharedSpatial = TRUE, prior.range = prior.range,
                           prior.sigma = prior.sigma)
model_setup$specifySpatial(Bias = 'inat_filtered', prior.range = inat.prior.range, prior.sigma = inat.prior.sigma)
model_setup$specifySpatial(Bias = 'DVC_filtered', prior.range = DVC.prior.range, prior.sigma = DVC.prior.sigma)

print("Model setup done")

# ------------------------------------------------------------------------------------------------------------------------
# Run model
# ------------------------------------------------------------------------------------------------------------------------
start.time <- Sys.time()
print(start.time)
model <- PointedSDMs::fitISDM(model_setup, options = list(num.threads = 16))
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
print(paste0("Model estimates: ", model$summary.hyperpar))
print(summary(model))

setwd("/users/jhw538/scratch/FC_deer/model_output_splines")

hyper_df <- data.frame(round(model[["summary.hyperpar"]], digits = 2))

if (any(c(hyper_df$mean[1], hyper_df$mean[3], hyper_df$mean[5]) > 600)){
stop("Estimated range larger than study area extent")}

if (any(c((hyper_df$sd[1]/hyper_df$mean[1]), (hyper_df$sd[3]/hyper_df$mean[3]), (hyper_df$sd[5]/hyper_df$mean[5])) < 0.1)){
stop("Standard deviation of estimated range suspiciously small")}

pdf(paste0("hyperpars_", sub(".RData", ".pdf", model_name)))
grid.table(hyper_df)
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# Plot model output
# ------------------------------------------------------------------------------------------------------------------------
required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / res_predictions)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / res_predictions)

data <- pixels(mesh = mesh, mask = ROI, nx = required_nx, ny = required_ny)

data <- cprod(data, data.frame(total_visit_duration = quantile(BBS_filtered$total_visit_duration, probs = 0.9)[[1]],
                               All_motor_vehicles = quantile(DVC_filtered$All_motor_vehicles, probs = 0.9)[[1]],
                               n_users_25km2 = quantile(inat_filtered$n_users_25km2, probs = 0.9)[[1]]))

formula = paste0("~ -1 + shared_spatial + ",  paste(offset_vars, collapse = " + "), " + ", paste(model[["names.fixed"]], collapse = " + "))
preds <- predict(model, data, fun = 'linear', n.samples = 1000, formula = as.formula(formula))

inat <- predict(model, data, predictor = TRUE, fun = 'linear', n.samples = 1000, formula = ~ inat_filtered_biasField)
DVC <- predict(model, data, predictor = TRUE, fun = 'linear', n.samples = 1000, formula = ~ DVC_filtered_biasField)
shared <- predict(model, data, predictor = TRUE, fun = 'linear', n.samples = 1000, formula = ~ shared_spatial)

pred_data <- data.frame(median = preds$predictions$median,
                        sd = preds$predictions$sd,
                        x = preds$predictions@coords[, 1],
                        y = preds$predictions@coords[, 2]) 

prediction_raster <- pred_data %>% 
      dplyr::select(x, y, median, sd) %>% 
      rast(., type = "xyz", crs = crs(km_projection))

cull_raster <- rasterize(validation_cull, prediction_raster, field = "count", fun = "sum") %>% 
      mask(., ROI_vect)

rast_stack <- c(prediction_raster$median, cull_raster)

rast_df <- rast_stack %>% 
      drop_na() %>% 
      data.frame()

rsq <- function (x, y) cor(x, y) ^ 2
r2 <- round(rsq(rast_df$median, log(rast_df$count_sum)), digits = 2)

lm_fit <- lm(log(count_sum) ~ median, data = rast_df)
intercept <- as.numeric(coef(lm_fit)[1])
slope <- as.numeric(coef(lm_fit)[2])

pred_inat <- data.frame(median = inat$predictions$median,
                        sd = inat$predictions$sd,
                        x = inat$predictions@coords[, 1],
                        y = inat$predictions@coords[, 2])

pred_DVC <- data.frame(median = DVC$predictions$median,
                       sd = DVC$predictions$sd,
                       x = DVC$predictions@coords[, 1],
                       y = DVC$predictions@coords[, 2])

pred_shared <- data.frame(median = shared$predictions$median,
                          sd = shared$predictions$sd,
                          x = shared$predictions@coords[, 1],
                          y = shared$predictions@coords[, 2])

species_scientific <- c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi")
species_common <- c("Red deer", "Roe deer", "Fallow deer", "Reeves's Muntjac")
spec_df <- data.frame(species_scientific, species_common)
spec_name <- spec_df$species_common[spec_df$species_scientific == species_choice]

pred_plot_median <- ggplot() +
      tidyterra::geom_spatvector(data = ROI_vect, fill = "white", colour = "black") +
      geom_spatraster(data = prediction_raster$median) +
      tidyterra::geom_spatvector(data = ROI_vect, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("log(Estimated relative abundance)") + 
      scale_fill_viridis(name = "Median", na.value="transparent",
                         option = "H") 

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

pred_inat_plot_median <- ggplot(pred_inat) +
      geom_tile(aes(x = x, y = y, fill = median)) +
      tidyterra::geom_spatvector(data = vect(ROI), fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("iNaturalist spatial field") +
      scale_fill_viridis(name = "Median")

pred_DVC_plot_median <- ggplot(pred_DVC) +
      geom_tile(aes(x = x, y = y, fill = median)) +
      tidyterra::geom_spatvector(data = vect(ROI), fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("DVC spatial field") +
      scale_fill_viridis(name = "Median")

pred_shared_plot_median <- ggplot(pred_shared) +
      geom_tile(aes(x = x, y = y, fill = median)) +
      tidyterra::geom_spatvector(data = vect(ROI), fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Shared spatial field") +
      scale_fill_viridis(name = "Median")

corrplot <- ggplot(rast_df, aes(x = mean, y = log(count_sum))) +
      geom_point() +
      geom_abline(intercept = intercept, slope = slope, 
                  col = "blue", alpha = 0.5, linewidth = 1.2) +
      theme_bw() +
      ggtitle('Correlation: Estimate - observed') +
      xlab("log(Estimated relative abundance)") +
      ylab("log(Culling counts)") +
      annotate("text", x = min(rast_df$mean) + 0.2,
               y = max(log(rast_df$count_sum)) - 0.5,
               label = paste0("Intercept = ", round(intercept, 2), 
                              "\nSlope = ", round(slope, 2),
                              "\nR\u00B2 = ", r2), 
               size = 3, hjust = 0)

# pred_plot <-  (pred_plot_mean + pred_plot_sd + pred_plot_lower + pred_plot_upper) / (pred_shared_plot_median + pred_inat_plot_median + pred_DVC_plot_median) + 
#       plot_annotation(title = paste0(spec_name, " (", sub(".RData", "", model_name), ")"),
#                       theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))+ 
#       plot_layout(nrow = 2)

pred_plot <- egg::ggarrange(pred_plot_median, pred_plot_sd, corrplot, pred_shared_plot_median, pred_inat_plot_median, pred_DVC_plot_median, nrow = 2)
pred_plot_annotated <- ggpubr::annotate_figure(pred_plot, top = ggpubr::text_grob(spec_name, face = "bold", size = 14))

pdf(paste0("maps_", sub(".RData", ".pdf", model_name)), width = 18, height = 9)
pred_plot_annotated
dev.off()

model_fixed <- data.frame(model$summary.fixed) %>% 
      mutate(ID = rownames(.))
model_fixed$significant <- ifelse((model_fixed$X0.025quant > 0 & model_fixed$X0.975quant > 0)|(model_fixed$X0.025quant < 0 & model_fixed$X0.975quant < 0), "yes", "no")
model_fixed <- model_fixed[!grepl("intercept", model_fixed$ID), ]

pdf(paste0("forest_", sub(".RData", ".pdf", model_name)))
ggplot() + 
      geom_point(data = model_fixed, aes(y = ID, x = X0.5quant, col = significant)) +
      geom_errorbar(data = model_fixed, aes(y = ID, xmin = X0.025quant, xmax = X0.975quant, col = significant), width = 0.1) +
      geom_vline(aes(xintercept = 0), lty = 2, alpha = .3) +
      ylab("Posterior estimates") +
      xlab(element_blank()) +
      theme_minimal() +
      scale_colour_manual(values = c("black", "#D55E00"), guide = "none")
dev.off()

