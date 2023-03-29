# Analysing deer distribution and abundance in England
R scripts to estimate the distribution of the four most common deer species in England, namely roe deer, Reeve's muntjac, fallow deer and red deer. This analysis is applying species distribution models with data integration and spline regression in INLA, using the R package 'pointedSDMs'.

## Table of Contents
| Script Name      | Description |
| ---------------- | ----------- |
| [Prepare_observation_data.R](R/Prepare_observation_data.R) | Preprocessing and cleaning of observation data |
| [Visualize_raw_data.R](R/Visualize_raw_data.R) | Plot the raw deer observation data |
| [Prepare_covariate_layers.R](R/Prepare_covariate_layers.R) | Preprocessing of environmental covariates |
| [Prepare_INLA_data.R](R/Prepare_INLA_data.R) | Get observation and covariate data ready for the INLA SDM |
| [Run_INLA_model.R](R/Run_INLA_model.R) | Apply the INLA model using the 'pointedSDMs' package, export results |
| [Get_INLA_results.R](R/Get_INLA_results.R) | Iterate through INLA model outputs and produce results plots |
| [impact_density_comparison.R](R/impact_density_comparison.R) | Analyse the deer density-impact relationship |
| [Folder: England_deer_shiny](England_deer_shiny) | Create a Shiny app to explore model outputs |

## Acknowledgements
This work was funded by the Forestry Commission (project FEE/1086). We would like to thank the organisations and individuals that kindly provided us with data for the development of the models. Deer vehicle collision data and helpful insights came from Dr Jochen Langbein, with additional data from National Highways. Deer count data from the Breeding Bird Survey were used by permission of the BTO. The BTO/JNCC/RSPB Breeding Bird Survey is a partnership jointly funded by the BTO, JNCC and RSPB, with fieldwork conducted by volunteers. Culling data were provided by Forestry England, and Woodland Impact Survey data came from the records of the Deer Initiative. Gamebag data came from the GWCT’s National Gamebag Census.
