rm(list = ls())
library(PointedSDMs)
library(inlabru)
library(tidyverse)
library(tidyterra)
library(INLA)
library(terra)
library(rgdal)
library(patchwork)
km_projection <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

# ------------------------------------------------------------------------------------------------------------------------
# Global variables, load data
# ------------------------------------------------------------------------------------------------------------------------
load('model-files/Cervus_elaphus_E23_R122.4_0.5_S4.9_0.1.RData')
load('data-processed/lincomb_files.RData')

res_predictions <- 2

ROI_vect <- vect('spatial_data/england_outline_27700_simple.shp') %>% 
      terra::project(crs(km_projection))

# ------------------------------------------------------------------------------------------------------------------------
# Run model
# ------------------------------------------------------------------------------------------------------------------------
# Subset lincombs variable to covariates used in model
lc_names <- unique(c(sapply(names(covariates_selection), function(x) gsub("_[0-9].*$", "", x))))
modified_lc_names <- c(outer(lc_names, 1:100, function(x, y) paste0(x, "_sortBase_lc", y)))

all_lc_sub <- all_lc[names(all_lc) %in% modified_lc_names]

names(covariates_selection)
start.time <- Sys.time()
print(start.time)
model <- PointedSDMs::fitISDM(model_setup, options = list(lincomb = all_lc_sub,
                                                          control.compute = list(cpo = T)))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

logCPO <- round(-sum(log(model$cpo$cpo[model$cpo$cpo != 0])), digits = 2)

# ------------------------------------------------------------------------------------------------------------------------
# Get model predictions
# ------------------------------------------------------------------------------------------------------------------------
required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / res_predictions)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / res_predictions)

data <- pixels(mesh = mesh, mask = ROI, nx = required_nx, ny = required_ny)

data <- cprod(data, data.frame(total_visit_duration = quantile(BBS_filtered$total_visit_duration, probs = 0.9)[[1]],
                               All_motor_vehicles = quantile(DVC_filtered$All_motor_vehicles, probs = 0.9)[[1]],
                               n_users_10km = quantile(inat_filtered$n_users_10km, probs = 0.9)[[1]]))

formula <- paste0("~ shared_spatial + n_users_100km2 + All_motor_vehicles +", paste(model[["names.fixed"]], collapse = " + "))

preds <- predict(model, data, fun = 'linear', n.samples = 1000, formula = as.formula(formula))

prediction_raster <- rast(preds[["predictions"]])

# ------------------------------------------------------------------------------------------------------------------------
# Get linear combination estimates
# ------------------------------------------------------------------------------------------------------------------------
original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
      mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
             sequence = as.numeric(gsub("\\D", "", rownames(.))))

effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                           sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                           quant_05 = model$summary.lincomb.derived$`0.5quant`,
                           quant_0025 = model$summary.lincomb.derived$`0.025quant`,
                           quant_0975 = model$summary.lincomb.derived$`0.975quant`)
effect_combs$covariate <- gsub('_sortBase', '', effect_combs$covariate)

effect_combs_m <- merge(original_values, effect_combs) %>%
      mutate(orig_values = if_else(str_detect(covariate, "log"), exp(orig_values) - 1, orig_values),
             species = species_choice)

# ------------------------------------------------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------------------------------------------------
writeRaster(prediction_raster, paste0("model-output/predRaster_", gsub("RData", "tif", model_name)), overwrite=TRUE)
save(model, species_choice, model_name, file = paste0("model-output/modelOut_", model_name))
write.csv(effect_combs_m, paste0("model-output/linCombs_", gsub("RData", "csv", model_name)))
