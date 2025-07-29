# -------------------------------------------------------------
# Preprocess forcing data for mHM.
# Place `preprocess_config.yaml` inside the folder pointed to by
# `domain_path` and ensure `use_python()` references an mHM-enabled
# Python environment.
# Run with: source("R/run_preprocessing.R")
# -------------------------------------------------------------

# List of required packages
required_packages <- c(
  "reticulate", "yaml", "terra", "sf", "magrittr",
  "tidyverse", "whitebox", "lubridate"
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

#Load helper functions
source("R/utils.R")
source("R/preprocess/preprocess_forcings.R")

# 1.set the path for the domain where all outputs will be written
domain_path = "/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/domain_chile"

# set up virtual environment with python packages.
# Recommended: use a virtual environment with mhm installation
  # instructions in: https://mhm-ufz.org/guides/install-unix/
# Additional packages needed
  # subprocess, os, pyyaml, re
# use_condaenv("mhm_env", required = TRUE)
# Specify the path to the Python interpreter
use_python("/Users/mhm/miniforge3/envs/mhm_env/bin/python", required = TRUE)

# The configuration file must be inside the domain_path
config_path <- file.path(domain_path, "preprocess_config.yaml")
config = read_yaml(config_path)

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
preprocess_geo_data(domain_path, remove_temp = FALSE, source.file = "global")
# Create and appply common mask
create_roi_mask(domain_path, remove_temp = TRUE, basin.mask = TRUE)

# At this point we should check if the facc.asc coincides with the gauge points
#   in the est_fluv folder. Run the following functions when idgauges and facc
#   are well adjusted
# Process streamflow data
preprocess_streamflow_data(domain_path, remove_temp = TRUE)
# Create gauge raster
create_idgauges(domain_path, remove_temp = TRUE)

# Load python functions
source_python("python/run_create_latlon.py")
source_python("python/create_geoparameter_block.py")
source_python("python/modify_parameters_with_geoblock.py")
source_python("python/update_nml_gauges.py")
source_python("python/update_time_periods.py")
source_python("python/calibrate_model_option.py")

# Call the function defined in the Python script
# create latlon.nc
run_create_latlon(domain_path)
# system2("cp", args = c("~/Downloads/temp_latlon/latlon.nc", "/Volumes/KINGSTON/domain_7339001/latlon/"), stdout = TRUE, stderr = TRUE)
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