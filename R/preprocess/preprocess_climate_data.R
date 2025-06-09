library(sf)
library(terra)
library(jsonlite)
library(ncdf4)
library(abind)
source("scripts/R/utils.R")

# Load config
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")

roi_path <- config$roi_file
date1 <- config$start_date
date2 <- config$end_date
remove_temp <- config$remove_temp
variables <- config$variables_clim
cellsize_clim <- config$cellsize_clim
cellunit <- config$cellunit
clim_folder <- config$meteo_folder
header_folder <- config$header_folder

# Prepare folders
dir.create(clim_folder, showWarnings = FALSE, recursive = TRUE)
temp_folder <- file.path(clim_folder, "temp")
dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)

# Read ROI
roi <- st_read(roi_path, quiet = TRUE)
box = st_bbox(roi)
box.m = as.matrix(rbind(
  c(box[[1]], box[[2]]),
  c(box[[1]], box[[4]]),
  c(box[[3]], box[[4]]),
  c(box[[3]], box[[2]]),
  c(box[[1]], box[[2]])
))
box.pol = st_polygon(list(box.m)) %>% st_sfc(crs = 4326)
st_write(box.pol, file.path(temp_folder, "extent.geojson"), quiet = TRUE, append = FALSE)

# Buffered region
buffered_box <- st_buffer(box.pol, dist = 6000) # 6 km buffer
st_write(buffered_box, file.path(temp_folder, "buffered.geojson"), quiet = TRUE, append = FALSE)

bbox <- st_bbox(buffered_box)
X1 <- bbox["xmin"]; Y1 <- bbox["ymin"]; X2 <- bbox["xmax"]; Y2 <- bbox["ymax"]

mhm_varnames = c("pre","tmin","tmax","tavg","pet")
names(mhm_varnames) = names(variables)
# Loop over variables
for (VAR in names(variables)) {
  paths <- variables[[VAR]]
  input_dir <- paths$input_dir
  if (is.null(input_dir)) next
  
  input_files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE)
  cat("🔄 Processing variable:", VAR, "\n")
  cat("  📂 Found", length(input_files), "files in", input_dir, "\n")
  
  r = rast(input_files, subds = VAR)
  # cut to buffer
  r.cropped = crop(r, buffered_box)
  # select timeframe
  r.tf = subset(r.cropped, time(r.cropped) >= as.Date(date1) & time(r.cropped) <= as.Date(date2))
  varnames(r.tf) = "pre"
  ref = rast(r.tf)
  res(ref) = cellsize_clim
  r.rs = resample(r.tf, ref, method = "bilinear")
  
  resampled_path = file.path(temp_folder, paste0(mhm_varnames[[VAR]],"_resampled.nc"))
  writeCDF(r.rs, filename = resampled_path, overwrite = TRUE)
  
  # Trim to final bounding box
  buffered_3km <- st_buffer(box.pol, dist = 3000)
  r.rs.crop = crop(r.rs, buffered_3km)
  final_path = file.path(clim_folder, paste0(mhm_varnames[[VAR]],".nc"))
  writeCDF(r.rs.crop, filename = final_path, overwrite = TRUE, prec = "double", missval = -9999)

  # Header
  ncols <- ncol(r.rs.crop)
  nrows <- nrow(r.rs.crop)
  xllcorner <- terra::ext(r.rs.crop)[1]
  yllcorner <- terra::ext(r.rs.crop)[3]
  nodata_value <- -9999
  header_text <- sprintf(
    "ncols         %d\nnrows         %d\nxllcorner     %f\nyllcorner     %f\ncellsize      %f\nNODATA_value  %d\n",
    ncols, nrows, xllcorner, yllcorner, cellsize_clim, nodata_value
  )
  
  writeLines(header_text, file.path(clim_folder, "header.txt"))
  writeLines(header_text, file.path(header_folder, "header_meteo.txt"))
  cat("✅ Header saved.\n")
  
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}