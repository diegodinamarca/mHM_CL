
process_LAI = function(config, remove_temp = FALSE){
  library(terra)
  library(sf)
  library(jsonlite)
  library(tidyverse)
  # source("scripts/R/utils.R")
  
  # === Load config ===
  lai_path <- config$lai_file
  roi_path <- config$roi_file
  target_crs <- paste0("EPSG:", config$latlon_crs)
  target_res <- config$cellsize
  lai_folder <- config$lai_folder
  temp_folder <- file.path(lai_folder, "temp")
  dir.create(temp_folder, recursive = TRUE, showWarnings = FALSE)
  cellsize = config$cellsize
  
  # Output paths
  clipped_path <- file.path(temp_folder, "clipped_lai.tif")
  clipped_filled_path <- file.path(temp_folder, "clipped_filled_lai.tif")
  reprojected_path <- file.path(temp_folder, "reprojected_lai.tif")
  netcdf_path <- file.path(lai_folder, "lai.nc")
  ref_file = file.path(config$meteo_folder, "pre.nc")
  
  # === Load ROI ===
  roi <- st_read(roi_path, quiet = TRUE)
  box = get_extent(roi)
  bbox = st_buffer(box, dist = 5000)
  # === Load LAI raster ===
  lai <- rast(lai_path)
  
  # Clip LAI to ROI
  lai_clipped <- crop(lai, vect(bbox))
  # Save clipped original
  writeRaster(lai_clipped, clipped_path, overwrite = TRUE)
  
  # Fill no-data values
  lai_filled <- focal_repeat(lai_clipped*0.1, n = 10)
  
  # Save filled raster
  writeRaster(lai_filled, clipped_filled_path, overwrite = TRUE)
  
  # resample raster to ref grid
  ref = rast(ref_file)[[1]]
  res(ref) = cellsize
  lai_rs = resample(lai_filled, ref, method = "bilinear")
  names(lai_rs) = paste0("lai_",1:12)
  NAflag(lai_rs) = -9999
  writeCDF(lai_rs,
           filename = netcdf_path,
           varname = "lai",
           unit = "m2/m2",
           overwrite = TRUE,
           compression = 1,
           longname = "long term monthly mean LAI",
           history = "Created with terra::writeCDF()",
           prec = "double",
           missval = -9999
  )
  
  cat("Done! Final NetCDF written to:\n", netcdf_path, "\n")
}

