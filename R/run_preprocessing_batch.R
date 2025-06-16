# -------------------------------------------------------------
# Batch preprocessing for mHM.
# Expects `preprocess_config.json` inside the domain folder supplied
# as a command-line argument and a Python environment referenced in
# `use_python()`.
# Usage: Rscript R/run_preprocessing_batch.R <domain_path>
# -------------------------------------------------------------

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

# Load helper functions
source("scripts/R/utils.R")
source("scripts/R/preprocess_forcings.R")

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
config_path <- file.path(domain_path, "preprocess_config.json")
config <- read_json(config_path)
# Create folders inside the domain for writing outputs
DirPaths <- c(config$out_folder, config$header_folder, config$latlon_folder,
              config$morph_folder, config$meteo_folder, config$lai_folder,
              config$lc_folder, config$gauge_folder)
for (d in DirPaths) dir.create(file.path(domain_path, d), showWarnings = FALSE)

# Process data
preprocess_climate_data(domain_path, remove_temp = TRUE, large.raster = FALSE, fix.negatives = TRUE)
preprocess_LAI_data(domain_path, remove_temp = TRUE, iter_num = 10)
preprocess_dem_data(domain_path, remove_temp = TRUE)
preprocess_lc_data(domain_path, remove_temp = TRUE)
preprocess_soil_data(domain_path, remove_temp = TRUE, iter_num = 10)
preprocess_geo_data(domain_path, remove_temp = TRUE, source.file = "global")
create_roi_mask(domain_path, remove_temp = TRUE, basin.mask = FALSE)
preprocess_streamflow_data(domain_path, remove_temp = TRUE)
create_idgauges(domain_path, remove_temp = TRUE)

# Call Python scripts
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/run_create_latlon.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/create_geoparameter_block.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/modify_parameters_with_geoblock.py")
source_python("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/python/update_nml_gauges.py")
run_create_latlon(domain_path)
write_geoparam_block(domain_path)
update_geoparam_block(domain_path)
update_evaluation_gauges(domain_path)

