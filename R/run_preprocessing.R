
# main script

# the purpose of this script is to proces the forcings and create the inputs for
# running the mHM model
library(reticulate)

# 1.read general configuration this configuration file should be the only input
# for the rest of the processing scripts
source("scripts/R/utils.R")
source("scripts/R/preprocess_forcings.R")
config <- "/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json"

# 2.process climate data
preprocess_climate_data(config, remove_temp = FALSE)
# 3.process lai data
preprocess_LAI_data(config, remove_temp = FALSE)
# 4.process dem data
preprocess_dem_data(config, remove_temp = FALSE)
# 5.process landcover data
preprocess_lc_data(config, remove_temp = FALSE)
# 6.process soil data
preprocess_soil_data(config, remove_temp = FALSE)
# 7.process geology data
preprocess_geo_data(config, remove_temp = FALSE)
# 8.create and appply common mask
create_roi_mask(config, remove_temp = FALSE, basin.mask = TRUE)
# 9.process streamflow data
preprocess_streamflow_data(config, remove_temp = FALSE)
# 10.create gauge raster
create_idgauges(config, remove_temp = FALSE)
# 11.create latlon file
use_condaenv("mhm_env", required = TRUE)

# # Run a whole Python script
py_run_file("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/run_create_latlon.py")
