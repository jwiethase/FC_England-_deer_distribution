rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(ggpubr)
library(viridis)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
options(scipen = 999)
# ------------------------------------------------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------------------------------------------------
ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))

traffic_data_2010_2020 <- read.csv("data-raw/England_traffic_data_2010_2020.csv")

load('data-processed/deer_records.RData')

inat <- inaturalist_final_2004_2020 %>% 
      filter(species == "Dama dama") 

DVC <- DVC_final_2010_2020 %>% 
      filter(species == "Dama dama") 

BDS_final_2005_2016 <- BDS_final_2005_2016 %>% 
      filter(species == "Dama dama") 

BBS_1995_2022_final <- BBS_1995_2022_final %>% 
      filter(species == "Dama dama") 


# ------------------------------------------------------------------------------------------------------------------------
# Plot of raw observation data
# ------------------------------------------------------------------------------------------------------------------------
BDS_plot <- ggplot(BDS_final_2005_2016) +
      geom_point(aes(x = x, y = y, col = as.factor(presence)), cex = 0.9, alpha = 0.8, pch = 15)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("British Deer Survey") +
      scale_color_manual(name = "Deer presence:", values = c("#D55E00", "#009E73"))

BBS_plot <- ggplot(BBS_1995_2022_final) +
      geom_point(aes(x = x, y = y, col = as.factor(presence)), cex = 0.9, alpha = 0.8, pch = 15)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("BTO Breeding Bird Survey") +
      scale_color_manual(name = "Deer presence:", values = c("#D55E00", "#009E73"))

inat_plot <- ggplot(inat) +
      geom_point(aes(x = x, y = y), col = "#E69F00", cex = 0.9, alpha = 0.8)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("iNaturalist")

DVC_plot <- ggplot(DVC) +
      geom_point(aes(x = x, y = y), col = "#56B4E9", cex = 0.9, alpha = 0.8)  +
      geom_spatvector(data = ROI, fill = NA, colour = "black") +
      theme_bw() +
      theme(legend.position = "bottom") +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Deer Vehicle Collisions")

final_plot <- ggpubr::ggarrange(DVC_plot, BDS_plot, BBS_plot, inat_plot, nrow = 1, common.legend = TRUE, legend="bottom")
ggsave(plot = final_plot, filename = "figures/red_deer_raw_data_plot.png", width = 25, height = 10, units = "cm")

# ------------------------------------------------------------------------------------------------------------------------
# Density of DVC records
# ------------------------------------------------------------------------------------------------------------------------
OS_grid_15km <- vect('spatial_data_raw/OSGB_Grid_1km.shp') %>% 
      terra::project(crs(km_projection)) %>% 
      mask(., ROI) %>% 
      rast(., resolution = 15) 

DVC_vect <- DVC %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)
DVC_rast <- rasterize(DVC_vect, OS_grid_15km, fun=function(i){length(i)}, background=0) %>% mask(., ROI) 

png("figures/DVC_density.png", width = 15, height = 15, units = "cm", res = 150)
plot(DVC_rast, main = "Red deer DVC records density")
plot(ROI, add = T)
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# Density and values of traffic data
# ------------------------------------------------------------------------------------------------------------------------
traffic_vect_2010_2020 <- traffic_data_2010_2020 %>% vect(geom = c("x", "y"), crs = crs(km_projection), keepgeom = FALSE)
traffic_point_rast <- rasterize(traffic_vect_2010_2020, OS_grid_15km, fun=function(i){length(i)}, background=0) %>% mask(., ROI) 
traffic_data_rast <-  rasterize(traffic_vect_2010_2020, OS_grid_15km, field = "All_motor_vehicles", fun=sum)  %>% mask(., ROI) 

png("figures/traffic_point_density.png", width = 15, height = 15, units = "cm", res = 150)
plot(traffic_point_rast, main = "Traffic data point density (2010-2020)")
plot(ROI, add = T)
dev.off()

png("figures/traffic_data.png", width = 15, height = 15, units = "cm", res = 150)
plot(traffic_data_rast, main = "Number of motor vehicles (2010-2020)")
plot(ROI, add = T)
dev.off()

png("figures/DVC_road_types.png", width = 20, height = 15, units = "cm", res = 150)
ggplot(DVC %>% mutate(Road_name = sub("^([[:alpha:]]*).*", "\\1", Road_name))) +
      geom_histogram(aes(x = Road_name), stat = "count") +
      theme_bw() +
      xlab("Road type") +
      ylab("Number of DVC records") +
      ggtitle("DVC Red deer records by road type")
dev.off()

png("figures/DVC_vehicle_density.png", width = 20, height = 15, units = "cm", res = 150)
plot(density(DVC$All_motor_vehicles), main = "Density of traffic counts associated with DVC Red deer records")
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# Estimated abundance
# ------------------------------------------------------------------------------------------------------------------------
load("model-output/HPC_Cervus_elaphus_E25_R150_0.05_S150_0.05.RData")
species_scientific <- c("Cervus elaphus", "Capreolus capreolus", "Dama dama", "Muntiacus reevesi")
species_common <- c("Red deer", "Roe deer", "Fallow deer", "Reeves's Muntjac")
spec_df <- data.frame(species_scientific, species_common)
spec_name <- spec_df$species_common[spec_df$species_scientific == species_choice]

png("figures/Red_deer_relative_abundance.png", width = 15, height = 15, units = "cm", res = 150)
ggplot(pred_data) +
      geom_tile(aes(x = x, y = y, fill = rescaled_median)) +
      tidyterra::geom_spatvector(data = vect(ROI), fill = NA, colour = "black") +
      theme_bw() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Relative abundance") + 
      scale_fill_viridis(name = element_blank())
dev.off()





