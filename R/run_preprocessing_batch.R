# -------------------------------------------------------------
# Batch preprocessing for mHM.
# Expects `preprocess_config.yaml` inside the domain folder supplied
# as a command-line argument and a Python environment referenced in
# `use_python()`.
# Usage: Rscript R/run_preprocessing_batch.R <domain_path>
# example: Rscript R/run_preprocessing_batch.r /Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/domain_11342001
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

# Load helper functions
source("./R/utils.R")
source("./R/preprocess/preprocess_forcings.R")

# Obtain domain path from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript R/run_preprocessing_batch.R <domain_path>")
}
domain_path <- args[1]

# set up virtual environment with python packages.
# use_condaenv("mhm_env", required = TRUE)
use_python("/Users/mhm/miniforge3/envs/mhm_env/bin/python", required = TRUE)

# The configuration file must be inside the domain_path
config_path <- file.path(domain_path, "preprocess_config.yaml")
config <- read_yaml(config_path)
# Create folders inside the domain for writing outputs
DirPaths <- c(config$out_folder, config$header_folder, config$latlon_folder,
              config$morph_folder, config$meteo_folder, config$lai_folder,
              config$lc_folder, config$gauge_folder)
for (d in DirPaths) dir.create(file.path(domain_path, d), showWarnings = FALSE)

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
create_roi_mask(domain_path, remove_temp = TRUE, basin.mask = TRUE)
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

