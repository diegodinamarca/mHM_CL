
# main script

# the purpose of this script is to proces the forcings and create the inputs for
# running the mHM model
# List of required packages
required_packages <- c(
  "reticulate", "jsonlite", "terra", "sf", "magrittr", 
  "tidyverse", "whitebox", "lubridate"
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

#Load helper functions
source("scripts/R/utils.R")
source("scripts/R/preprocess_forcings.R")

# 1.set the path for the domain where all outputs will be written
domain_path = "/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/domain_zone_7"

# set up virtual environment with python packages.
# Recommended: use a virtual environment with mhm installation
  # instructions in: https://mhm-ufz.org/guides/install-unix/
# Additional packages needed
  # subprocess, os, json, re
# use_condaenv("mhm_env", required = TRUE)
# Specify the path to the Python interpreter
use_python("/Users/mhm/miniforge3/envs/mhm_env/bin/python", required = TRUE)

# The configuration file must be inside the domain_path
config_path <- file.path(domain_path, "preprocess_config.json")
config = read_json(config_path)

# Create folders inside the domain for writing outputs
dir.create(file.path(domain_path, config$out_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$header_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$latlon_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$morph_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$meteo_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$lai_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$lc_folder), showWarnings = FALSE)
dir.create(file.path(domain_path, config$gauge_folder), showWarnings = FALSE)


# Process climate data
preprocess_climate_data(domain_path, remove_temp = TRUE, large.raster = FALSE,
                        fix.negatives = TRUE, time_origin = config$start_date)
# Process lai data
preprocess_LAI_data(domain_path, remove_temp = TRUE, iter_num = 10, force.fill = TRUE)
# Process dem data
preprocess_dem_data(domain_path, remove_temp = TRUE)
# Process landcover data
preprocess_lc_data(domain_path, remove_temp = TRUE)
# Process soil data
preprocess_soil_data(domain_path, remove_temp = TRUE, iter_num = 10)
# Process geology data
preprocess_geo_data(domain_path, remove_temp = TRUE, source.file = "global")
# Create and appply common mask
create_roi_mask(domain_path, remove_temp = TRUE, basin.mask = FALSE)
# Process streamflow data
preprocess_streamflow_data(domain_path, remove_temp = TRUE)
# Create gauge raster
create_idgauges(domain_path, remove_temp = TRUE)

# Load python functions
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/run_create_latlon.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/create_geoparameter_block.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/modify_parameters_with_geoblock.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/update_nml_gauges.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/update_time_periods.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/calibrate_model_option.py")

# Call the function defined in the Python script
# create latlon.nc
run_create_latlon(domain_path)
# write geoparameter block for mhm_parameter.nml
write_geoparam_block(domain_path)
# update mhm_parameter.nml with the geoparameter block
update_geoparam_block(domain_path)
# update mhm.nml with the gauges inside the gauges folder
update_evaluation_gauges(domain_path)
# update simulation/calibration time period
update_time_periods(domain_path)
# add an option to calibrate the model modifying mhm.nml
update_calibrate_option(domain_path)
# PENDING - Shuffle initial parameter values
