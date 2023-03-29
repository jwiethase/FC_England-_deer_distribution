rm(list = ls())
library(viridis)
library(tidyverse)
library(tidyterra)
library(terra)
library(patchwork)
library(grid)
library(gridExtra)
library(geodata)

# ------------------------------------------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------------------------------------------
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
rsq <- function (x, y) cor(x, y) ^ 2

ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))

species_scientific <- c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi")
species_common_list <- c("Red deer", "Roe deer", "Fallow deer", "Muntjac")
poly_file <- c("red_deer_density_polygons.shp", "roe_deer_density_polygons.shp", "fallow_deer_density_polygons.shp", "muntjac_density_polygons.shp")

gb <- gadm(country = "GB", level = 1, resolution = 2,
           path = "spatial_data/")

effects_list <- list()

for(i in 1:length(species_scientific)){
      species_choice <- species_scientific[4]
      species_common <- species_common_list[species_scientific == species_choice]
      
      # ------------------------------------------------------------------------------------------------------------------------
      # Load data
      # ------------------------------------------------------------------------------------------------------------------------
      prediction_raster <- rast(list.files("model-output", full.names = T)[grep(paste0("predRaster_", gsub(" ", "_", species_choice)), list.files("model-output", full.names = T))])
      med_rast <- prediction_raster$median
      range01 <- function(x){(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))}
      values(med_rast) <- range01(values(med_rast))
      # writeRaster(med_rast, filename = paste0(gsub(" ", "_", species_choice), "_rast.tif"))

      # load(list.files("model-output", full.names = T)[grep(paste0("modelOut_", gsub(" ", "_", species_choice)), list.files("model-output", full.names = T))])
      effects_data <- read.csv(list.files("model-output", full.names = T)[grep(paste0("linCombs_", gsub(" ", "_", species_choice)), list.files("model-output", full.names = T))])
      effects_list[[i]] <- effects_data
      
      validation_cull <- read.csv('data-processed/cull_final_2017_2020.csv') %>% 
            filter(species == species_choice) %>% 
            vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE) 
      
      validation_GBC <- read.csv('data-processed/GBC_final_2005_2009.csv') %>% 
            filter(species == species_choice) %>% 
            vect(geom = c("easting", "northing"), crs = crs("epsg:27700"), keepgeom = FALSE) %>% 
            terra::project(crs(km_projection))
      
      # Density data to adjust relative abundance estimates
      deer_density_polygons <- vect(paste0('spatial_data/', poly_file[species_scientific == species_choice])) %>% 
            terra::project(crs(km_projection))
      
      # Make a background layer for spatial context
      bbox <- extend(ext(ROI), c(+50, +50, +50, +50))
      
      gb_ROI <- gb %>% 
            terra::project(crs(km_projection)) %>% 
            crop(bbox)
      
      # ------------------------------------------------------------------------------------------------------------------------
      # Model validation - Forestry England culling data
      # ------------------------------------------------------------------------------------------------------------------------
      cull_raster <- rasterize(validation_cull, prediction_raster$median, field = "count", fun = "sum") %>% 
            mask(., ROI) 
      
      rast_stack <- c(prediction_raster$median, cull_raster)
      
      rast_df <- rast_stack %>% 
            drop_na() %>% 
            data.frame()
      print(paste0(species_choice, "_cull_N: ", NROW(rast_df)))
      r2 <- rsq(rast_df$median, log(rast_df$count_sum)); r2
      
      # Fit linear regression model
      lm_fit <- lm(log(count_sum) ~ median, data = rast_df)
      intercept <- as.numeric(coef(lm_fit)[1])
      slope <- as.numeric(coef(lm_fit)[2]); slope
      pearson <- cor(rast_df$median, log(rast_df$count_sum)); pearson
      
      corr_plot <- ggplot(rast_df, aes(x = median, y = log(count_sum))) +
            geom_point() +
            geom_abline(intercept = intercept, slope = slope, 
                        col = "blue", alpha = 0.5, linewidth = 1.2) +
            theme_bw() +
            ggtitle('Validation - Forestry England culling') +
            xlab("log(Estimated relative abundance)") +
            ylab("log(Culling counts)") +
            annotate("text", x = min(rast_df$median) + 0.2,
                     y = max(log(rast_df$count_sum)) - 0.5,
                     label = paste0("Pearson's R = ", round(pearson, 2),
                                    "\nR\u00B2 = ", round(r2, 2)), 
                     size = 3, hjust = 0); corr_plot
      
      # assign(paste0(gsub(" ", "_", species_choice), "_cull_plot"), corr_plot, envir = .GlobalEnv)
      
      # ------------------------------------------------------------------------------------------------------------------------
      # Model validation - Gab Bag Census data
      # ------------------------------------------------------------------------------------------------------------------------
      OS_grid_100km2 <- vect('spatial_data_raw/OSGB_Grid_1km.shp') %>% 
            terra::project(crs(km_projection)) %>% 
            mask(., ROI) %>% 
            rast(., resolution = 10) 
      
      GBC_rast <- terra::rasterize(validation_GBC, OS_grid_100km2, field = "NGC", background=0) %>% mask(., ROI) 
      
      preds_raster_100km2 <- resample(prediction_raster$median, GBC_rast, method = 'average') %>% mask(., ROI) 
      
      rast_stack_GBC <- c(preds_raster_100km2, GBC_rast)
      
      rast_df_GBC <- rast_stack_GBC %>% 
            drop_na() %>% 
            data.frame() %>% 
            filter(NGC_last > 0)
      
      r2_GBC <- rsq(rast_df_GBC$median, log(rast_df_GBC$NGC_last)); r2_GBC
      print(paste0(species_choice, "_GBC_N: ", NROW(rast_df_GBC)))
      # Fit linear regression model
      lm_fit_GBC <- lm(log(NGC_last) ~ median, data = rast_df_GBC)
      intercept_GBC <- as.numeric(coef(lm_fit_GBC)[1])
      slope_GBC <- as.numeric(coef(lm_fit_GBC)[2]); slope_GBC
      pearson_GBC <- cor(rast_df_GBC$median, log(rast_df_GBC$NGC_last)); pearson_GBC
      
      corr_plot_GBC <- ggplot(rast_df_GBC, aes(x = median, y = log(NGC_last))) +
            geom_point() +
            geom_abline(intercept = intercept_GBC, slope = slope_GBC, 
                        col = "blue", alpha = 0.5, linewidth = 1.2) +
            theme_bw() +
            ggtitle('Validation - Game Bag Census') +
            xlab("log(Estimated relative abundance)") +
            ylab("log(Game Bag Census Index)") +
            annotate("text", x = min(rast_df_GBC$median) + 0.1,
                     y = max(log(rast_df_GBC$NGC_last)) - 0.2,
                     label = paste0("Pearson's R = ", round(pearson_GBC, 2),
                                    "\nR\u00B2 = ", round(r2_GBC, 2)), 
                     size = 3, hjust = 0); corr_plot_GBC
      # assign(paste0(gsub(" ", "_", species_choice), "_GBC_plot"), corr_plot_GBC, envir = .GlobalEnv)
      
      validation_plots <- gridExtra::grid.arrange(corr_plot, corr_plot_GBC, nrow = 1, top = grid::textGrob(species_common, gp=grid::gpar(fontsize=15,font=8)))
      
      # assign(paste0(gsub(" ", "_", species_choice), "_val_plots"), validation_plots, envir = .GlobalEnv)
      ggsave(plot = validation_plots, filename = paste0("figures/", gsub(" ", "_", species_choice), "val_plots.png"), width = 25, height = 15, units = "cm")

      # ------------------------------------------------------------------------------------------------------------------------
      # Adjust abundance
      # ------------------------------------------------------------------------------------------------------------------------
      raster_values <- terra::extract(prediction_raster, deer_density_polygons, fun = mean, na.rm = TRUE)
      deer_density_polygons$raster_median <- raster_values$q0.5
      deer_density_polygons$raster_LCL <- raster_values$q0.025
      deer_density_polygons$raster_UCL <- raster_values$q0.975
      
      extraction_df <- as.data.frame(deer_density_polygons) %>% filter(!is.nan(density), !is.nan(raster_median))
      
      lm_fit_median <- lm(log(density) ~ raster_median, data = extraction_df)
      intercept_median <- as.numeric(coef(lm_fit_median)[1])
      slope_median <- as.numeric(coef(lm_fit_median)[2])
      
      r2_adjust <- round(rsq(extraction_df$raster_median, log(extraction_df$density)), digits = 2); r2_adjust
      pearson_adjust <- cor(extraction_df$raster_median, log(extraction_df$density)); pearson_adjust
      
      corr_plot_adjust <- ggplot(extraction_df, aes(x = raster_median, y = log(density))) +
            geom_point() +
            geom_abline(intercept = intercept_median, slope = slope_median, 
                        col = "blue", alpha = 0.5, linewidth = 1.2) +
            theme_bw() +
            ggtitle('Correlation: Estimate - observed') +
            xlab("log(Estimated relative abundance)") +
            ylab("log(ind/km2)") +
            annotate("text", x = min(extraction_df$raster_median) + 0.01,
                     y = max(log(extraction_df$density)) - 0.1,
                     label = paste0("Intercept = ", round(intercept_median, 2), 
                                    "\nSlope = ", round(slope_median, 2),
                                    "\nR\u00B2 = ", round(r2_adjust, 2)), 
                     size = 3, hjust = 0); corr_plot_adjust
      # assign(paste0(gsub(" ", "_", species_choice), "_adjust_plot"), corr_plot_adjust, envir = .GlobalEnv)
      
      prediction_raster$adjusted_abundance_median <- exp((prediction_raster$q0.5 * slope_median) + intercept_median)
      writeRaster(prediction_raster$adjusted_abundance_median, filename = paste0("model-output/", gsub(" ", "_", species_choice), "_adjusted.tif"), overwrite = T)
      
      # ------------------------------------------------------------------------------------------------------------------------
      # Plot model output - smoothed spatial estimates
      # ------------------------------------------------------------------------------------------------------------------------
      pred_plot_median <- ggplot() + 
            tidyterra::geom_spatvector(data = gb_ROI, fill = "white", colour = "black") +
            tidyterra::geom_spatvector(data = ROI, fill = "white", colour = "black") +
            geom_spatraster(data = med_rast) +
            tidyterra::geom_spatvector(data = ROI, fill = NA, colour = "black") +
            theme_bw() +
            xlab(element_blank()) +
            ylab(element_blank()) +
            ggtitle(species_common) + 
            scale_fill_viridis(name = "Median", na.value="transparent",
                               option = "H") +
            coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE); pred_plot_median
      
      assign(paste0(gsub(" ", "_", species_choice), "_pred_plot_median"), pred_plot_median, envir = .GlobalEnv)
      
      pop_size_median <- round(sum(values(prediction_raster$adjusted_abundance_median), na.rm = TRUE)*4); pop_size_median

      pred_plot_median_adjusted <- ggplot() +
            tidyterra::geom_spatvector(data = gb_ROI, fill = "white", colour = "black") +
            tidyterra::geom_spatvector(data = ROI, fill = "white", colour = "black") +
            geom_spatraster(data = prediction_raster$adjusted_abundance_median) +
            tidyterra::geom_spatvector(data = ROI, fill = NA, colour = "black") +
            theme_bw() +
            xlab(element_blank()) +
            ylab(element_blank()) +
            ggtitle("Adjusted abundance") +
            scale_fill_viridis(name = "ind/km\u00B2", na.value="transparent",
                               option = "H") +
            coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE); pred_plot_median_adjusted

      pred_plot_sd <- ggplot() +
            tidyterra::geom_spatvector(data = gb_ROI, fill = "white", colour = "black") +
            tidyterra::geom_spatvector(data = ROI, fill = "white", colour = "black") +
            geom_spatraster(data = prediction_raster$sd) +
            tidyterra::geom_spatvector(data = ROI, fill = NA, colour = "black") +
            theme_bw() +
            xlab(element_blank()) +
            ylab(element_blank()) +
            scale_fill_viridis(name = "Median", na.value="transparent",
                               option = "H") +
            coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE); pred_plot_sd
      
      assign(paste0(gsub(" ", "_", species_choice), "_pred_plot_sd"), pred_plot_sd, envir = .GlobalEnv)
      
      # pred_plot <- pred_plot_median | pred_plot_sd + plot_annotation(title = species_choice, tag_levels = 'A')
      # assign(paste0(gsub(" ", "_", species_choice), "_pred_plots"), pred_plot, envir = .GlobalEnv)
      adjusted_abundance_plots <- gridExtra::grid.arrange(corr_plot_adjust, pred_plot_median_adjusted, nrow = 1, top = grid::textGrob(species_common, gp=grid::gpar(fontsize=15,font=8)))
      assign(paste0(gsub(" ", "_", species_choice), "_adjusted_abundance"), adjusted_abundance_plots, envir = .GlobalEnv)
      
      # ------------------------------------------------------------------------------------------------------------------------
      # Plot model output - effect plots
      # ------------------------------------------------------------------------------------------------------------------------
      # Effect plots, for regression spline fits
      combined_data <- data.frame(values = c(effects_data$quant_05, 
                                             effects_data$quant_0025,
                                             effects_data$quant_0975))
      combined_data$rescaled_value <- range01(combined_data$values)
      n <- nrow(effects_data)
      effects_data$median_rescaled <- combined_data$rescaled_value[1:n]
      effects_data$quant_0025_rescaled <- combined_data$rescaled_value[(n + 1):(2 * n)]
      effects_data$quant_0975_rescaled <- combined_data$rescaled_value[(2 * n + 1):(3 * n)]
      
      effects_plot <- ggplot(effects_data) +
            geom_line(aes(x = orig_values, y = median_rescaled)) +
            geom_line(aes(x = orig_values, y = quant_0025_rescaled), lty = 2, alpha = .5) +
            geom_line(aes(x = orig_values, y = quant_0975_rescaled), lty = 2, alpha = .5) +
            theme_bw() +
            ggtitle(species_common) +
            facet_wrap( ~ covariate, scale = 'free', ncol = 4, labeller = labeller(
                  covariate = c(arable = "Arable cover [%]", 
                                builtup = "Built-up cover [%]", 
                                elevation = "Elevation [m]", 
                                tree_cover_density = "Tree cover density [%]")
            )) +
            xlab(element_blank()) +
            ylab("Relative abundance"); effects_plot
      assign(paste0(gsub(" ", "_", species_choice), "_effects_plot"), effects_plot, envir = .GlobalEnv)
      
      ggsave(plot = effects_plot, filename = paste0("figures/", gsub(" ", "_", species_choice), "_effects_plots.png"), width = 25, height = 8, units = "cm")
      
}

# ------------------------------------------------------------------------------------------------------------------------
# Plot model output - effect plots
# ------------------------------------------------------------------------------------------------------------------------
effects_data_combined <- bind_rows(effects_list)
combined_data <- data.frame(values = c(effects_data_combined$quant_05, 
                                       effects_data_combined$quant_0025,
                                       effects_data_combined$quant_0975))
combined_data$rescaled_value <- range01(combined_data$values)
n <- nrow(effects_data_combined)
effects_data_combined$median_rescaled <- combined_data$rescaled_value[1:n]
effects_data_combined$quant_0025_rescaled <- combined_data$rescaled_value[(n + 1):(2 * n)]
effects_data_combined$quant_0975_rescaled <- combined_data$rescaled_value[(2 * n + 1):(3 * n)]
effects_data_combined <- effects_data_combined %>% 
      rowwise() %>% 
      mutate(species_common = species_common_list[species_scientific == species])
      
      
# Effect plots, for regression spline fits
effects_plot_all <- ggplot(effects_data_combined) +
      geom_line(aes(x = orig_values, y = median_rescaled)) +
      geom_line(aes(x = orig_values, y = quant_0025_rescaled), lty = 2, alpha = .5) +
      geom_line(aes(x = orig_values, y = quant_0975_rescaled), lty = 2, alpha = .5) +
      theme_bw() +
      facet_grid(species_common ~ covariate, scale = 'free', labeller = labeller(
            covariate = c(arable = "Arable cover [%]", 
                          builtup = "Built-up area [%]", 
                          elevation = "Elevation [m]", 
                          tree_cover_density = "Tree cover density [%]")
      )) +
      xlab(element_blank()) +
      ylab("Relative abundance"); effects_plot_all

ggsave(plot = effects_plot_all, filename ="figures/effects_plot_all.png", width = 16, height = 11, units = "cm", dpi = 300)

# ------------------------------------------------------------------------------------------------------------------------
# Combined map plots
# ------------------------------------------------------------------------------------------------------------------------
comb_maps_med <- ggpubr::ggarrange(Capreolus_capreolus_pred_plot_median, Cervus_elaphus_pred_plot_median, Dama_dama_pred_plot_median, Muntiacus_reevesi_pred_plot_median,
                               common.legend = TRUE, legend = "right",
               ncol = 4,
               nrow = 1)

comb_maps_sd <- ggpubr::ggarrange(Capreolus_capreolus_pred_plot_sd, Cervus_elaphus_pred_plot_sd, Dama_dama_pred_plot_sd, Muntiacus_reevesi_pred_plot_sd,
                                   common.legend = TRUE, legend = "right",
                                   ncol = 4,
                                   nrow = 1)

comb_maps <- comb_maps_med / comb_maps_sd + plot_annotation(tag_levels = list(c("Relative abundance", "Standard deviation")))

png("figures/comb_maps.png", width = 23, height = 15, units = "cm", res = 200)
comb_maps
dev.off()
 
ggsave(plot = comb_maps_med, filename = paste0("figures/comb_maps_med.png"), width = 20, height = 20, units = "cm")


adjusted_abundance_comb <- gridExtra::grid.arrange(Capreolus_capreolus_adjusted_abundance, Muntiacus_reevesi_adjusted_abundance, nrow = 2)
ggsave(plot = adjusted_abundance_comb, filename = paste0("figures/adjusted_plots_combined.png"), width = 20, height = 21, units = "cm", dpi = 300)




