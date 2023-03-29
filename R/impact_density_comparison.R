library(tidyverse)
library(terra)
library(patchwork)
library(INLA)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

# ------------------------------------------------------------------------------------------------------------------------
# Global variables, load data
# ------------------------------------------------------------------------------------------------------------------------
ROI <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))
ROI_rgdal <- rgdal::readOGR('spatial_data/england_outline_27700_simple.shp') %>% 
      spTransform(., km_projection)

range01 <- function(x){(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))}

roe_deer_median <- rast("model-output/predRasterMedian_Capreolus_capreolus_E25_R256.7_0.5_S10.3_0.1.tif")
values(roe_deer_median) <- range01(values(roe_deer_median))
red_deer_median <- rast("model-output/predRasterMedian_Cervus_elaphus_E16_R122.4_0.5_S7.3_0.1.tif")
values(red_deer_median) <- range01(values(red_deer_median))
fallow_deer_median <- rast("model-output/predRasterMedian_Dama_dama_E14_R240.4_0.5_S9.6_0.1.tif")
values(fallow_deer_median) <- range01(values(fallow_deer_median))
muntjac_median <- rast("model-output/predRasterMedian_Muntiacus_reevesi_E24_R161.5_0.5_S4.9_0.1.tif")
values(muntjac_median) <- range01(values(muntjac_median))

impacts_vect <- vect('spatial_data/DIAdataall.shp') %>% 
      as.data.frame(.) %>% 
      # filter(Act_score >= 3) %>% 
      group_by(xcoord, ycoord) %>% 
      summarize(Impact_sco = mean(Impact_sco, na.rm = T)) %>% 
      ungroup() %>% 
      mutate(xcoord = as.numeric(xcoord),
             ycoord = as.numeric(ycoord)) %>% 
      select('xcoord', 'ycoord', 'Impact_sco') %>% 
      vect(., geom = c('xcoord', 'ycoord'), crs = crs("epsg:27700")) %>% 
      terra::project(crs(km_projection))  %>% 
      terra::crop(ROI)

# ------------------------------------------------------------------------------------------------------------------------
# Extract raster values at impact data points
# ------------------------------------------------------------------------------------------------------------------------
extr_roe_deer <- terra::extract(roe_deer_median, impacts_vect)
extr_red_deer <- terra::extract(red_deer_median, impacts_vect)
extr_fallow_deer <- terra::extract(fallow_deer_median, impacts_vect)
extr_muntjac <- terra::extract(muntjac_median, impacts_vect)

impacts_vect$roe_deer <- extr_roe_deer$median
impacts_vect$red_deer <- extr_red_deer$median
impacts_vect$fallow_deer <- extr_fallow_deer$median
impacts_vect$muntjac <- extr_muntjac$median

roe_df <- data.frame("Impact_sco" = impacts_vect$Impact_sco, "median_rel_abund" = impacts_vect$roe_deer,
                     "x" = crds(impacts_vect)[, 1], "y" = crds(impacts_vect)[, 2]) %>% drop_na()
red_df <- data.frame("Impact_sco" = impacts_vect$Impact_sco, "median_rel_abund" = impacts_vect$red_deer,
                     "x" = crds(impacts_vect)[, 1], "y" = crds(impacts_vect)[, 2]) %>% drop_na()
fallow_df <- data.frame("Impact_sco" = impacts_vect$Impact_sco, "median_rel_abund" = impacts_vect$fallow_deer,
                        "x" = crds(impacts_vect)[, 1], "y" = crds(impacts_vect)[, 2]) %>% drop_na()
muntjac_df <- data.frame("Impact_sco" = impacts_vect$Impact_sco, "median_rel_abund" = impacts_vect$muntjac,
                         "x" = crds(impacts_vect)[, 1], "y" = crds(impacts_vect)[, 2]) %>% drop_na()

# Simple INLA models
model_roe_deer_gauss <- inla(Impact_sco ~ median_rel_abund, data = roe_df,
                             control.compute = list(cpo = TRUE),
                             control.predictor=list(compute=TRUE))
model_red_deer_gauss <- inla(Impact_sco ~ median_rel_abund, data = red_df,
                             control.compute = list(cpo = TRUE),
                             control.predictor=list(compute=TRUE))
model_fallow_deer_gauss <- inla(Impact_sco ~ median_rel_abund, data = fallow_df, 
                                control.compute = list(cpo = TRUE),
                                control.predictor=list(compute=TRUE))
model_muntjac_gauss <- inla(Impact_sco ~ median_rel_abund, data = muntjac_df,
                            control.compute = list(cpo = TRUE),
                            control.predictor=list(compute=TRUE))

# Models with spatial random effect in inlabru
# Make the meshes and spde
max.edge = 25
region.bdry <- inla.sp2segment(ROI_rgdal)

mesh_roe <- inla.mesh.2d(loc = cbind(roe_df$x, roe_df$y),
                     boundary = region.bdry, 
                     cutoff = max.edge/2, 
                     max.edge = c(max.edge, max.edge*4), 
                     offset = c(max.edge, max.edge*5))
mesh_roe$crs <- km_projection
spde_roe <- inla.spde2.pcmatern(mesh_roe, prior.range = c(50, 0.5), prior.sigma = c(5, 0.1))

mesh_red <- inla.mesh.2d(loc = cbind(red_df$x, red_df$y),
                         boundary = region.bdry, 
                         cutoff = max.edge/2, 
                         max.edge = c(max.edge, max.edge*4), 
                         offset = c(max.edge, max.edge*5))
mesh_red$crs <- km_projection
spde_red <- inla.spde2.pcmatern(mesh_red, prior.range = c(50, 0.5), prior.sigma = c(5, 0.1))

mesh_fallow <- inla.mesh.2d(loc = cbind(fallow_df$x, fallow_df$y),
                         boundary = region.bdry, 
                         cutoff = max.edge/2, 
                         max.edge = c(max.edge, max.edge*4), 
                         offset = c(max.edge, max.edge*5))
mesh_fallow$crs <- km_projection
spde_fallow <- inla.spde2.pcmatern(mesh_fallow, prior.range = c(50, 0.5), prior.sigma = c(5, 0.1))

mesh_muntjac <- inla.mesh.2d(loc = cbind(muntjac_df$x, muntjac_df$y),
                         boundary = region.bdry, 
                         cutoff = max.edge/2, 
                         max.edge = c(max.edge, max.edge*4), 
                         offset = c(max.edge, max.edge*5))
mesh_muntjac$crs <- km_projection
spde_muntjac <- inla.spde2.pcmatern(mesh_muntjac, prior.range = c(50, 0.5), prior.sigma = c(5, 0.1))

c.c <- list(cpo = TRUE,
            dic = TRUE,
            waic = TRUE,
            config = TRUE)

roe_bru <- bru(Impact_sco ~ median_rel_abund + myfield(cbind(roe_df$x, roe_df$y), model = spde_roe),
            data = roe_df,
            family = "gaussian",
            options = list(control.compute = c.c))

red_bru <- bru(Impact_sco ~ median_rel_abund + myfield(cbind(red_df$x, red_df$y), model = spde_red),
               data = red_df,
               family = "gaussian",
               options = list(control.compute = c.c))

fallow_bru <- bru(Impact_sco ~ median_rel_abund + myfield(cbind(fallow_df$x, fallow_df$y), model = spde_fallow),
               data = fallow_df,
               family = "gaussian",
               options = list(control.compute = c.c))

muntjac_bru <- bru(Impact_sco ~ median_rel_abund + myfield(cbind(muntjac_df$x, muntjac_df$y), model = spde_muntjac),
               data = muntjac_df,
               family = "gaussian",
               options = list(control.compute = c.c))

roe_bru$summary.fixed
red_bru$summary.fixed
fallow_bru$summary.fixed
muntjac_bru$summary.fixed

# Model comparison
hist(roe_bru$cpo$pit)
hist(model_roe_deer_gauss$cpo$pit)

qqplot(qunif(ppoints(length(roe_bru$cpo$pit))),
       roe_bru$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(roe_bru$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
qqplot(qunif(ppoints(length(model_roe_deer_gauss$cpo$pit))),
       model_roe_deer_gauss$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_roe_deer_gauss$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
-sum(log(roe_bru$cpo$cpo))
-sum(log(model_roe_deer_gauss$cpo$cpo))

hist(red_bru$cpo$pit)
hist(model_red_deer_gauss$cpo$pit)

qqplot(qunif(ppoints(length(red_bru$cpo$pit))),
       red_bru$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(red_bru$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
qqplot(qunif(ppoints(length(model_red_deer_gauss$cpo$pit))),
       model_red_deer_gauss$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_red_deer_gauss$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
-sum(log(red_bru$cpo$cpo))
-sum(log(model_red_deer_gauss$cpo$cpo))

hist(fallow_bru$cpo$pit)
hist(model_fallow_deer_gauss$cpo$pit)

qqplot(qunif(ppoints(length(fallow_bru$cpo$pit))),
       fallow_bru$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(fallow_bru$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
qqplot(qunif(ppoints(length(model_fallow_deer_gauss$cpo$pit))),
       model_fallow_deer_gauss$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_fallow_deer_gauss$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
-sum(log(fallow_bru$cpo$cpo))
-sum(log(model_fallow_deer_gauss$cpo$cpo))

hist(muntjac_bru$cpo$pit)
hist(model_muntjac_gauss$cpo$pit)

qqplot(qunif(ppoints(length(muntjac_bru$cpo$pit))),
       muntjac_bru$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(muntjac_bru$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
qqplot(qunif(ppoints(length(model_muntjac_gauss$cpo$pit))),
       model_muntjac_gauss$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_muntjac_gauss$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
-sum(log(muntjac_bru$cpo$cpo))
-sum(log(model_muntjac_gauss$cpo$cpo))

roe_preds <- predict(roe_bru, 
                      data = data.frame(median_rel_abund = seq(0,1, length = 1000)), 
                      formula = Impact_sco ~ median_rel_abund + Intercept + myfield,
                      n.samples = 1000)
red_preds <- predict(red_bru, 
                     data = data.frame(median_rel_abund = seq(min(red_df$median_rel_abund, na.rm = T), max(red_df$median_rel_abund, na.rm = T), length = 100)), 
                     formula = Impact_sco ~ median_rel_abund + Intercept,
                     n.samples = 1000)
fallow_preds <- predict(fallow_bru, 
                     data = data.frame(median_rel_abund = seq(min(fallow_df$median_rel_abund, na.rm = T), max(fallow_df$median_rel_abund, na.rm = T), length = 100)), 
                     formula = Impact_sco ~ median_rel_abund + Intercept,
                     n.samples = 1000)
muntjac_preds <- predict(muntjac_bru, 
                        data = data.frame(median_rel_abund = seq(min(muntjac_df$median_rel_abund, na.rm = T), max(muntjac_df$median_rel_abund, na.rm = T), length = 100)), 
                        formula = Impact_sco ~ median_rel_abund + Intercept,
                        n.samples = 1000)

# Plot the outputs
roe_plot <- ggplot() +
      geom_point(data = roe_df, aes(x = median_rel_abund, y = Impact_sco), alpha = 0.5) +
      geom_line(data = roe_preds,
                aes(x = median_rel_abund, y = q0.5), lty = 1, alpha = 0.9, lwd = 1.2, col = "coral") +
      geom_ribbon(data = roe_preds,
                  aes(x = median_rel_abund, ymin = q0.025, ymax = q0.975), 
                  alpha = 0.1) +
      theme_linedraw() +
      ggtitle('Roe deer') +
      xlab(element_blank()) +
      ylab("Impact score")
roe_plot_rug <- ggExtra::ggMarginal(roe_plot, type = "densigram", fill = "lightblue")

red_plot <- ggplot() +
      geom_point(data = red_df, aes(x = median_rel_abund, y = Impact_sco), alpha = 0.5) +
      geom_line(data = red_preds,
                aes(x = median_rel_abund, y = q0.5), lty = 1, alpha = 0.9, lwd = 1.2, col = "coral") +
      geom_ribbon(data = red_preds,
                  aes(x = median_rel_abund, ymin = q0.025, ymax = q0.975), 
                  alpha = 0.1) +
      theme_linedraw() +
      ggtitle('Red deer') +
      xlab(element_blank()) +
      ylab(element_blank())
red_plot_rug <- ggExtra::ggMarginal(red_plot, type = "densigram", fill = "lightblue")

fallow_plot <- ggplot() +
      geom_point(data = fallow_df, aes(x = median_rel_abund, y = Impact_sco), alpha = 0.5) +
      geom_line(data = fallow_preds,
                aes(x = median_rel_abund, y = q0.5), lty = 1, alpha = 0.9, lwd = 1.2, col = "coral") +
      geom_ribbon(data = fallow_preds,
                  aes(x = median_rel_abund, ymin = q0.025, ymax = q0.975), 
                  alpha = 0.1) +
      theme_linedraw() +
      ggtitle('Fallow deer') +
      xlab("Relative abundance") +
      ylab("Impact score")
fallow_plot_rug <- ggExtra::ggMarginal(fallow_plot, type = "densigram", fill = "lightblue")

muntjac_plot <- ggplot() +
      geom_point(data = muntjac_df, aes(x = median_rel_abund, y = Impact_sco), alpha = 0.5) +
      geom_line(data = muntjac_preds,
                aes(x = median_rel_abund, y = q0.5), lty = 1, alpha = 0.9, lwd = 1.2, col = "coral") +
      geom_ribbon(data = muntjac_preds,
                  aes(x = median_rel_abund, ymin = q0.025, ymax = q0.975), 
                  alpha = 0.1) +
      theme_linedraw() +
      ggtitle('Muntjac') +
      xlab("Relative abundance") +
      ylab(element_blank()) 
muntjac_plot_rug <- ggExtra::ggMarginal(muntjac_plot, type = "densigram", fill = "lightblue")

combined_plots <- gridExtra::grid.arrange(roe_plot_rug, red_plot_rug, fallow_plot_rug, muntjac_plot_rug)

ggsave(combined_plots, filename = "figures/impacts_plot2.png", width = 20, height = 20, units = "cm")

# Function to extract predictions and confidence intervals
extract_model_data <- function(model, species_name) {
      pred <- model$summary.fixed[, "0.5quant"]
      lower <- model$summary.fixed[, "0.025quant"]
      upper <- model$summary.fixed[, "0.975quant"]
      data.frame(species = species_name, pred = pred, lower = lower, upper = upper)[1, ]
}

# Extract data for each species
roe_deer_data <- extract_model_data(roe_bru, "Roe deer")
red_deer_data <- extract_model_data(red_bru, "Red deer")
fallow_deer_data <- extract_model_data(fallow_bru, "Fallow deer")
muntjac_data <- extract_model_data(muntjac_bru, "Muntjac")
combined_data <- rbind(roe_deer_data, red_deer_data, fallow_deer_data, muntjac_data)

forest_plot <- ggplot(combined_data, aes(x = species, y = pred, ymin = lower, ymax = upper)) +
      geom_pointrange() +
      labs(title = "Model Results", x = "Species", y = "Predicted Impact slope") +
      theme_bw() +
      geom_hline(yintercept = 0, lty = 2); forest_plot
