rm(list = ls())
library(tidyverse)
library(terra)
library(viridis)
library(egg)
library(tidyterra)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

# ------------------------------------------------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------------------------------------------------
load('model-output/Dama_dama_E13_R240.4_0.5_S9.6_0.1.RData')

validation_cull <- read.csv('data-processed/cull_final_2017_2020.csv') %>% 
      filter(species == species_choice) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE) 

deer_density <- read.csv('data-raw/deer_density.csv') %>% 
      filter(species == species_choice) %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE) 

ROI <- vect(ROI) %>% 
      terra::project(crs(km_projection))

# ------------------------------------------------------------------------------------------------------------------------
# Visualize predictions and validation data
# ------------------------------------------------------------------------------------------------------------------------
prediction_raster <- pred_data %>% 
      dplyr::select(x, y, median, sd) %>% 
      rast(., type = "xyz", crs = crs(km_projection))

cull_raster <- rasterize(validation_cull, prediction_raster, field = "count", fun = "sum") %>% 
      mask(., ROI) 

pred_plot_median <- ggplot() +
      tidyterra::geom_spatvector(data = ROI, fill = "white", colour = "black") +
      geom_spatraster(data = prediction_raster$median) +
      tidyterra::geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("log(Estimated relative abundance)") + 
      scale_fill_viridis(name = element_blank(), na.value="transparent",
                         option = "H"); pred_plot_median

cull_plot <- ggplot() +
      tidyterra::geom_spatvector(data = ROI, fill = "white", colour = "black") +
      geom_spatraster(data = cull_raster) +
      tidyterra::geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Culling counts") + 
      scale_fill_viridis(name = "Deer count", na.value="transparent",
                         option = "H"); cull_plot

# ------------------------------------------------------------------------------------------------------------------------
# Correlate estimates with culling data
# ------------------------------------------------------------------------------------------------------------------------
rast_stack <- c(prediction_raster$median, cull_raster)

rast_df <- rast_stack %>% 
      drop_na() %>% 
      data.frame()

rsq <- function (x, y) cor(x, y) ^ 2

r2 <- rsq(rast_df$median, log(rast_df$count_sum)); r2

# ------------------------------------------------------------------------------------------------------------------------
# Rescale estimates using culling data
# ------------------------------------------------------------------------------------------------------------------------
# Fit linear regression model
lm_fit <- lm(log(count_sum) ~ median, data = rast_df)
intercept <- as.numeric(coef(lm_fit)[1])
slope <- as.numeric(coef(lm_fit)[2])

corr_plot <- ggplot(rast_df, aes(x = median, y = log(count_sum))) +
      geom_point() +
      geom_abline(intercept = coef(lm_fit)[1], slope = coef(lm_fit)[2], 
                  col = "blue", alpha = 0.5, linewidth = 1.2) +
      theme_bw() +
      ggtitle('Correlation: Estimate - observed') +
      xlab("log(Estimated relative abundance)") +
      ylab("log(Culling counts)") +
      annotate("text", x = min(rast_df$median) + 0.2,
               y = max(log(rast_df$count_sum)) - 0.5,
               label = paste0("Intercept = ", round(intercept, 2), 
                              "\nSlope = ", round(slope, 2),
                              "\nR = ", round(pearson, 2)), 
               size = 3, hjust = 0); corr_plot
      
# Rescale estimated log relative abundance to true abundance
rast_stack$density_1km2 <- exp((rast_stack$median * slope) + intercept)/4 
pop_estimate_Eng <- round(sum(values(rast_stack$density_1km2), na.rm = T)*4, digits = 0)

species_scientific <- c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi")
species_common <- c("Red deer", "Roe deer", "Fallow deer", "Reeves's Muntjac")
spec_df <- data.frame(species_scientific, species_common)
spec_name <- spec_df$species_common[spec_df$species_scientific == species_choice]

estimated_density <- ggplot() +
      tidyterra::geom_spatvector(data = ROI, fill = "white", colour = "black") +
      geom_spatraster(data = rast_stack$density_1km2) +
      tidyterra::geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle(paste0("Estimated true density (", pop_estimate_Eng, ")")) + 
      scale_fill_viridis(name = "Deer density [ind/km\u00B2]", na.value="transparent",
                         option = "H"); estimated_density

final_plot <- egg::ggarrange(pred_plot_median, cull_plot, corr_plot, estimated_density, 
                       nrow = 2, labels = c("A", "B", "C", "D"))
final_plot_annotated <- ggpubr::annotate_figure(final_plot, top = ggpubr::text_grob(spec_name, face = "bold", size = 14))

png(paste0("figures/", gsub(" ", "_", species_choice), "_output_validation.png"), width = 25, height = 25, units = "cm", res = 300)
final_plot_annotated
dev.off()

png(paste0("figures/", gsub(" ", "_", species_choice), "_density.png"), width = 15, height = 15, units = "cm", res = 300)
ggpubr::annotate_figure(estimated_density, top = ggpubr::text_grob(spec_name, face = "bold", size = 14))
dev.off()

