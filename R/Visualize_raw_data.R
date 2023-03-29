rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

# ------------------------------------------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------------------------------------------
ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))

load('data-processed/deer_records.RData')

covariates <- raster::stack(readRDS('data-processed/raster_covariates.RDS'))
validation_cull <- read.csv('data-processed/cull_final_2017_2020.csv') %>% 
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE) 

species_list <- c("Red deer", "Roe deer", "Fallow deer", "Reeves's muntjac")

inat <- inaturalist_final_2004_2020 %>% 
      mutate(species = replace(species, species == "Capreolus capreolus", "Roe deer"),
             species = replace(species, species == "Cervus elaphus", "Red deer"),
             species = replace(species, species == "Dama dama", "Fallow deer"),
             species = replace(species, species == "Muntiacus reevesi", "Reeves's muntjac"),
             species = replace(species, species == "Cervus nippon", "Sika deer"),
             species = replace(species, species == "Hydropotes inermis", "Chinese water deer"),
             source = "iNaturalist",
             type = "Presence-only") %>% 
      dplyr::select(species, x, y, presence, source, type)

DVC <- DVC_final_2010_2020 %>% 
      mutate(species = replace(species, species == "Capreolus capreolus", "Roe deer"),
             species = replace(species, species == "Cervus elaphus", "Red deer"),
             species = replace(species, species == "Dama dama", "Fallow deer"),
             species = replace(species, species == "Muntiacus reevesi", "Reeves's muntjac"),
             species = replace(species, species == "Cervus nippon", "Sika deer"),
             species = replace(species, species == "Hydropotes inermis", "Chinese water deer"),
             source = "Deer Vehicle Collisions",
             type = "Presence-only",
             presence = 1) %>% 
      dplyr::select(species, x, y, presence, source, type)

BDS <- BDS_final_2005_2016 %>% 
      ungroup() %>% 
      mutate(species = replace(species, species == "Capreolus capreolus", "Roe deer"),
             species = replace(species, species == "Cervus elaphus", "Red deer"),
             species = replace(species, species == "Dama dama", "Fallow deer"),
             species = replace(species, species == "Muntiacus reevesi", "Reeves's muntjac"),
             species = replace(species, species == "Cervus nippon", "Sika deer"),
             species = replace(species, species == "Hydropotes inermis", "Chinese water deer"),
             source = "British Deer Survey",
             type = "Presence-absence") %>% 
      dplyr::select(species, x, y, presence, source, type)

BBS <- BBS_1995_2022_final %>% 
      mutate(species = replace(species, species == "Capreolus capreolus", "Roe deer"),
             species = replace(species, species == "Cervus elaphus", "Red deer"),
             species = replace(species, species == "Dama dama", "Fallow deer"),
             species = replace(species, species == "Muntiacus reevesi", "Reeves's muntjac"),
             source = "Breeding Bird Survey",
             type = "Presence-absence") %>% 
      dplyr::select(species, x, y, presence, source, type)

all_combined <- rbind(inat, DVC, BDS, BBS) %>% 
      filter(species %in% species_list)

# ------------------------------------------------------------------------------------------------------------------------
# Plot of raw observation data
# ------------------------------------------------------------------------------------------------------------------------
full_plot_panels <- ggplot(all_combined) +
      geom_point(aes(x = x, y = y, col = as.factor(presence)), pch = 15, cex = 0.7, alpha = 0.8)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      stat_density2d(aes(x = x, y = y), n = 100, alpha = .3, contour = T, col = "black") + 
      theme_bw() +
      theme(strip.text.x = element_text(size = 7.5),
            strip.text.y = element_text(size = 10)) +
      facet_grid(species ~ source) +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Raw deer observation data") +
      scale_color_manual(name = "Deer presence:", values = c("#D55E00", "#009E73")) +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme(legend.position = "bottom",
            legend.title = element_text(size=9),
            legend.text = element_text(size=9)) 

ggsave(plot = full_plot_panels, filename = "figures/raw_data_plot.png", width = 15, height = 18, units = "cm")      

for (i in 1:length(unique(all_combined$species))) {
      spec <- unique(all_combined$species)[i]
      full_plot_combined <- ggplot(all_combined %>% filter(species == spec)) +
            geom_point(aes(x = x, y = y, col = as.factor(presence)), cex = 0.4, alpha = 0.6)  +
            geom_spatvector(data = ROI, fill = NA, colour = "black") +
            stat_density2d(aes(x = x, y = y), n = 100, alpha = .3, contour = T, col = "black") + 
            theme_bw() +
            xlab(element_blank()) +
            ylab(element_blank()) +
            ggtitle("Raw observation data") +
            scale_color_manual(name = "Deer presence:", values = c("#D55E00", "#009E73")) +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.position = "bottom",
                  legend.title = element_text(size=14),
                  legend.text = element_text(size=13)) 
      ggsave(plot = full_plot_combined, filename = paste0("figures/", gsub(" ", "_", spec), "_raw_plot_combined.png"), width = 15, height = 15, units = "cm")      
      
}


red_deer_plot <- ggplot(all_combined %>% filter(species == "Red Deer")) +
      geom_point(aes(x = x, y = y, col = as.factor(presence), pch = type), cex = 0.9, alpha = 0.8)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      stat_density2d(aes(x = x, y = y), n = 100, alpha = .3, contour = T, col = "black") + 
      theme_bw() +
      facet_grid(. ~ source) +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Red deer observation data") +
      scale_color_manual(name = "Deer presence:", values = c("#D55E00", "#009E73")) +
      scale_shape_manual(name = "Data type:", values = c(18, 15)) +
      guides(shape = guide_legend(override.aes = list(size = 3)),
             color = guide_legend(override.aes = list(size = 3))) +
      theme(legend.position = "bottom",
            legend.title = element_text(size=12), 
            legend.text = element_text(size=11))

ggsave(plot = red_deer_plot, filename = "figures/red_deer_data_plot.png", width = 22, height = 9, units = "cm")    

# ------------------------------------------------------------------------------------------------------------------------
# Plot of covariate data
# ------------------------------------------------------------------------------------------------------------------------
pdf(file = "figures/covars_plot.pdf")
plot(covariates)
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# BDS effort is not great. Plot distribution of effort.
# ------------------------------------------------------------------------------------------------------------------------
ggplot(BDS_final_2005_2016) +
      geom_point(aes(x = x, y = y, col = as.factor(n_years_surveyed)), cex = 0.9, alpha = 0.8, pch = 15)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      facet_wrap(facets = vars(species), ncol = 4) +
      theme(legend.position = "bottom") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("British Deer Survey - Presence/Absence") +
      scale_color_manual(name = "Number of years surveyed:", values = c("#D55E00", "#009E73", "blue"))

# Mainly surveyed 3 times. Remove effort form model to simplify.


# ------------------------------------------------------------------------------------------------------------------------
# iNaturalist effort
# ------------------------------------------------------------------------------------------------------------------------
OS_grid_100km2 <- vect('spatial_data_raw/OSGB_Grid_1km.shp') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) %>% 
      rast(., resolution = 10) 

inat_vect <- inaturalist_final_2004_2020 %>% 
      mutate(distance_to_road = scale(log(distance_to_road+1)),
             n_users_100km2 = scale(log(n_users_100km2)),
             presence = NULL) %>%
      vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)
inat_rast <- rasterize(inat_vect, OS_grid_100km2, field = "n_users_100km2", background=0) %>% mask(., ROI) 

plot(inat_rast, main = "OS_grid_100km2")

# ------------------------------------------------------------------------------------------------------------------------
# Culling data
# ------------------------------------------------------------------------------------------------------------------------
ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))

OS_grid_15km <- vect('spatial_data_raw/OSGB_Grid_1km.shp') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) %>% 
      rast(., resolution = 15) 

cull_rast <- rasterize(validation_cull, OS_grid_15km, field = "count") %>% mask(., ROI) 
plot(cull_rast)














