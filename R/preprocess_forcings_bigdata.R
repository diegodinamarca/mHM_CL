preprocess_climate_data <- function(domain_path, remove_temp = TRUE) {
  library(jsonlite)
  library(terra)
  library(sf)
  library(magrittr)
  write_large_raster_to_netcdf <- function(rast, filename, varname,
                                                   batch_size = 4000,
                                                   units = "",
                                                   longname = "",
                                                   time_origin = "1960-01-01",
                                                   missval = -9999,
                                                   prec = "double",
                                                   fix.negatives = FALSE) {
    library(terra)
    library(ncdf4)
    
    # === Extract metadata ===
    time_vec <- terra::time(rast)
    if (is.null(time_vec)) stop("Raster must have a time dimension.")
    
    time_days <- as.integer(difftime(time_vec, as.Date(time_origin), units = "days"))
    lon <- terra::xFromCol(rast, 1:ncol(rast))
    lat <- terra::yFromRow(rast, 1:nrow(rast))
    n_time <- nlyr(rast)
    
    # === Define dimensions ===
    dim_lon <- ncdim_def("lon", "degrees_east", lon)
    dim_lat <- ncdim_def("lat", "degrees_north", lat)
    dim_time <- ncdim_def("time", paste0("days since ", time_origin), vals = time_days, unlim = TRUE)
    
    # === Define variable ===
    var_def <- ncvar_def(name = varname,
                         units = units,
                         dim = list(dim_lon, dim_lat, dim_time),
                         missval = missval,
                         longname = longname,
                         prec = prec)
    
    # === Create NetCDF ===
    nc <- nc_create(filename, vars = list(var_def))
    
    # === Write in batches ===
    batch_indices <- split(1:n_time, ceiling(seq_along(1:n_time) / batch_size))
    
    for (batch in batch_indices) {
      cat(sprintf("Writing layers %d to %d\n", batch[1], tail(batch, 1)))
      r_batch <- rast[[batch]]
      
      # === Round and fix negatives if requested ===
      r_batch <- round(r_batch, 2)
      if (fix.negatives) {
        r_batch <- clamp(r_batch, lower = 0)
      }
      
      a <- as.array(r_batch)
      a[is.na(a)] <- missval
      a <- aperm(a, c(2, 1, 3))  # [lon, lat, time]
      
      ncvar_put(nc, varid = varname, vals = a,
                start = c(1, 1, batch[1]), count = c(-1, -1, length(batch)))
    }
    
    # === Close file ===
    nc_close(nc)
    cat("✅ NetCDF successfully written to:", filename, "\n")
  }
  # Load config
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  
  roi_path <- config$roi_file
  date1 <- config$start_date
  date2 <- config$end_date
  variables <- config$variables_clim
  cellsize_clim <- config$cellsize_clim
  cellunit <- config$cellunit
  
  clim_folder <- file.path(domain_path, config$meteo_folder)
  header_folder <- file.path(domain_path, config$header_folder)
  clim_folder = "meteo"
  header_folder = "headers"
  # Prepare folders
  temp_folder <- file.path(clim_folder, "temp")
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Read ROI and build buffered polygon
  roi <- st_read(roi_path, quiet = TRUE)
  box <- st_bbox(roi)
  box.m <- as.matrix(rbind(
    c(box[[1]], box[[2]]),
    c(box[[1]], box[[4]]),
    c(box[[3]], box[[4]]),
    c(box[[3]], box[[2]]),
    c(box[[1]], box[[2]])
  ))
  box.pol <- st_polygon(list(box.m)) %>% st_sfc(crs = 4326)
  st_write(box.pol, file.path(temp_folder, "extent.geojson"), quiet = TRUE, append = FALSE)
  
  # Buffered region
  buffered_box <- st_buffer(box.pol, dist = 6000) # 6 km buffer
  st_write(buffered_box, file.path(temp_folder, "buffered.geojson"), quiet = TRUE, append = FALSE)
  
  bbox <- st_bbox(buffered_box)
  X1 <- bbox["xmin"]; Y1 <- bbox["ymin"]; X2 <- bbox["xmax"]; Y2 <- bbox["ymax"]
  
  mhm_varnames <- c("pre", "tmin", "tmax", "tavg", "pet")
  names(mhm_varnames) <- names(variables)
  
  for (VAR in names(variables)[2:5]) {
    # VAR = "pr"
    paths <- variables[[VAR]]
    input_dir <- paths$input_dir
    if (is.null(input_dir)) next
    
    input_files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE)
    cat("🔄 Processing variable:", VAR, "\n")
    cat("  📂 Found", length(input_files), "files in", input_dir, "\n")
    
    r <- rast(input_files, subds = VAR)
    # r <- crop(r, buffered_box)
    r.tf <- subset(r, time(r) >= as.Date(date1) & time(r) <= as.Date(date2))
    varnames(r.tf) <- mhm_varnames[[VAR]]
    ref <- rast(r.tf)
    res(ref) <- cellsize_clim
    r.rs <- resample(r.tf, ref, method = "bilinear")
    # r.rs[r.rs < 0] = 0
    # resampled_path <- file.path(temp_folder, paste0(mhm_varnames[[VAR]], "_resampled.nc"))
    # writeCDF(r.rs, filename = resampled_path, overwrite = TRUE)
    
    # buffered_3km <- st_buffer(box.pol, dist = 3000)
    # r.rs.crop <- crop(r.rs, buffered_3km)
    final_path <- file.path(clim_folder, paste0(mhm_varnames[[VAR]], ".nc"))
    if (VAR %in% c("pr","et0")){
      write_large_raster_to_netcdf(r.rs, final_path, varname = mhm_varnames[[VAR]], fix.negatives = TRUE)
      
    }else{
      write_large_raster_to_netcdf(r.rs, final_path, varname = mhm_varnames[[VAR]], fix.negatives = FALSE)
      
    }
    # writeCDF(r.rs, filename = final_path, overwrite = TRUE, prec = "double", missval = -9999)
    
    ncols <- ncol(r.rs)
    nrows <- nrow(r.rs)
    xllcorner <- terra::ext(r.rs)[1]
    yllcorner <- terra::ext(r.rs)[3]
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
}

preprocess_LAI_data = function(domain_path, remove_temp = FALSE, iter_num = 10){
  library(terra)
  library(sf)
  library(jsonlite)
  library(tidyverse)
  # source("scripts/R/utils.R")
  
  # === Load config ===
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  
  lai_path <- config$lai_file
  roi_path <- config$roi_file
  
  lai_folder <- file.path(domain_path, config$lai_folder)
  lai_folder= "lai"
  dir.create(lai_folder)
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
           
           prec = "double",
           missval = -9999
  )
  
  cat("Done! Final NetCDF written to:\n", netcdf_path, "\n")
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

preprocess_dem_data<- function(domain_path, remove_temp = FALSE) {
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  library(whitebox)
  
  # === Load config ===
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- jsonlite::fromJSON(config_path)
  
  ref_file <- file.path(domain_path, config$meteo_folder, "pre.nc")
  dem_file <- config$dem_file
  roi_file <- config$roi_file
  hydronet_file <- config$hydro_network_file
  cellsize <- config$cellsize
  
  header_folder = file.path(domain_path, config$header_folder)
  morph_folder <- file.path(domain_path, config$morph_folder)
  temp_folder <- file.path(morph_folder, "temp")
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
  # === Read inputs ===
  dem <- rast(dem_file)
  roi <- st_read(roi_file, quiet = TRUE)
  hydronet <- st_read(hydronet_file, quiet = TRUE) %>% st_transform(4326)
  
  box = get_extent(roi)
  roi_buff <- st_buffer(box, dist = 3000)
  # hydronet <- st_crop(hydronet, roi_buff)
  # hydronet_crop_file <- file.path(temp_folder, "hydronet.shp")
  
  # Cortar DEM al buffer
  dem <- crop(dem, vect(roi_buff))
  dem_file = file.path(temp_folder, "dem_clipped.tif")
  writeRaster(dem, filename = dem_file, overwrite = TRUE)
  
  # resample DEM to L0 resolution
  # load reference file
  ref <- rast(ref_file)[[1]]
  res(ref) = cellsize
  ref
  
  dem = resample(dem, ref, method = "bilinear")
  dem_file = file.path(temp_folder, "dem_resampled.tif")
  writeRaster(dem, filename = dem_file, overwrite = TRUE)
  
  wbt_init()
  smooth_dem = file.path(temp_folder, "smoothed.tif")
  # wbt_feature_preserving_smoothing(dem_file, output = smooth_dem, filter = 7)
  
  filled_dem = file.path(temp_folder, "filled_dem.tif")
  # wbt_breach_depressions(dem = smooth_dem, output = filled_dem)
  
  facc = file.path(temp_folder, "facc.tif")
  # wbt_d_inf_flow_accumulation(input = filled_dem, output = facc)
  
  fdir = file.path(temp_folder, "fdir.tif")
  wbt_flow_accumulation_full_workflow(dem = dem_file, 
                                      out_dem = filled_dem, 
                                      out_pntr = fdir, 
                                      out_accum = facc, 
                                      out_type = "cells",
                                      esri_pntr = TRUE, verbose_mode = TRUE)
  
  slope = file.path(temp_folder, "slope.tif")
  wbt_slope(dem = filled_dem, output = slope)
  aspect = file.path(temp_folder, "aspect.tif")
  wbt_aspect(dem = filled_dem, output = aspect)
  
  
  # Guardar versiones remuestreadas
  writeRaster(rast(filled_dem), filename = file.path(morph_folder, "dem.asc"), overwrite = TRUE,
              datatype = "FLT8S", NAflag = -9999)
  writeRaster(rast(slope), filename = file.path(morph_folder, "slope.asc"), overwrite = TRUE,
              datatype = "FLT8S", NAflag = -9999)
  writeRaster(rast(aspect), filename = file.path(morph_folder, "aspect.asc"), overwrite = TRUE,
              datatype = "FLT8S", NAflag = -9999)
  writeRaster(rast(fdir), filename = file.path(morph_folder, "fdir.asc"), overwrite = TRUE,
              datatype = "INT4S", NAflag = -9999)
  writeRaster(rast(facc), filename = file.path(morph_folder, "facc.asc"), overwrite = TRUE,
              datatype = "INT4S", NAflag = -9999)
  
  cat("PROCESAMIENTO COMPLETO. Archivos guardados en:", morph_folder, "\n")
  
  dem_path = file.path(morph_folder, "dem.asc")
  extract_asc_header(dem_path, header_folder)
  
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

preprocess_lc_data = function(domain_path, remove_temp = FALSE){
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  # leer archivo de configuracion json con rutas de archivos
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  
  # Paths from config
  raster_path <- config$land_cover_file
  roi_path <- config$roi_file
  
  lc_folder <- file.path(domain_path, config$lc_folder)
  ref_file <- file.path(domain_path, config$meteo_folder, "pre.nc")
  cellsize <- config$cellsize
  
  temp_folder <- file.path(lc_folder, "temp")
  if (!dir.exists(temp_folder)) dir.create(temp_folder, recursive = TRUE)
  
  # Reclassification dictionary
  reclass_dict = as.matrix(rbind(
    c(100,150,2),
    c(200,252,1),
    c(300,410,2),
    c(411,429,1),
    c(430,640,2),
    c(790,810,3),
    c(890,910,2),
    c(909,911,3),
    c(912,929,2),
    c(930,1200,3)
  ))
  reclass_dict
  # Load ROI
  roi <- st_read(roi_path, quiet = TRUE)
  box = get_extent(roi)
  
  # Crear buffer de 4km (4000 metros)
  roi_buff <- st_buffer(box, dist = 4000)
  
  # Step 1: Load LC raster and clip to ROI
  lc <- rast(raster_path)
  lc_clipped <- crop(lc, vect(st_transform(roi_buff, 32719)))
  
  # Save clipped raster
  clipped_path <- file.path(temp_folder, "landcover_clipped.tif")
  writeRaster(lc_clipped, clipped_path, overwrite = TRUE)
  
  # Step 2: Reclassify using dict
  lc_reclass = classify(lc_clipped, rcl = reclass_dict, include.lowest = T, right = T)
  lc_reclass[lc_reclass == 0] = NA
  # plot(lc_reclass)
  
  # Save reclassified raster
  reclass_path <- file.path(temp_folder, "landcover_reclassified.tif")
  writeRaster(lc_reclass, reclass_path, overwrite = TRUE)
  
  # Step 3: Reproject to EPSG:4326 (should already be)
  lc_reproj <- project(lc_reclass, "EPSG:4326", method = "near")
  
  # Step 4: Read DEM and align LC to its grid
  ref <- rast(ref_file)[[1]]
  res(ref) = cellsize
  # roi_buff <- st_buffer(box.pol, dist = 3000)
  # lc_3km= crop(lc_reproj, roi_buff)
  lc_rs <- resample(lc_reproj, ref, method = "near")
  
  # Step 5: Export final LC raster as ASCII Grid
  # Crear buffer de 5km (3000 metros)
  final_lc_path <- file.path(lc_folder, "landcover.asc")
  writeRaster(lc_rs, final_lc_path, overwrite = TRUE, NAflag = -9999,
              datatype = "INT4S")
  cat("Final landcover.asc saved to:", final_lc_path, "\n")
  
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}


preprocess_geo_data <- function(domain_path, remove_temp = FALSE, source.file = "local") {
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  # leer archivo de configuracion json con rutas de archivos
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  
  geo_file <- config$geo_file
  geo_file_global <- config$geo_file_global
  
  roi_file <- config$roi_file
  morph_folder <- file.path(domain_path, config$morph_folder)
  ref_file = file.path(domain_path, config$meteo_folder,"pre.nc")
  cellsize = config$cellsize
  
  temp_folder <- file.path(morph_folder, "temp")
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
  
  if (source.file == "local"){
    # Leer ROI
    roi <- st_read(roi_file, quiet = TRUE)
    box = get_extent(roi)
    bbox = st_buffer(box, dist = 3000)
    
    # Leer shapefile geológico
    geo <- st_read(geo_file, quiet = TRUE)
    
    # Recortar geología a la extensión de la ROI
    sf_use_s2(FALSE)
    geo_clipped <- st_crop(geo, bbox)
    sf_use_s2(TRUE)
    
    # Crear LUT para 'cd_geol' (string → número)
    unique_vals <- sort(unique(na.omit(geo_clipped$cd_geol)))
    val_map <- setNames(seq_along(unique_vals), unique_vals)
    geo_clipped$raster_val <- val_map[geo_clipped$cd_geol]
    
    # Grilla de referencia para rasterización
    ref = rast(ref_file)
    res(ref) = cellsize
    
    # Rasterizar
    geo_vect <- vect(geo_clipped)
    # writeVector(geo_vect, file.path(temp_dir, "roi_mask.shp"), overwrite=TRUE)
    rasterized <- rasterize(geo_vect, ref, field = "raster_val", touches = TRUE)
    # plot(rasterized)
    # values(rasterized)[is.na(values(rasterized))] <- -9999
    
    # Guardar raster como ASCII Grid
    raster_path <- file.path(morph_folder, "geology_class.asc")
    writeRaster(rasterized, raster_path, overwrite = TRUE, NAflag = -9999,
                datatype = "INT4S")
    cat("Process completed: geology raster generated from local shapefile.\n")
    
    
  }else if (source.file == "global"){
    roi <- st_read(roi_file, quiet = TRUE)
    r <- rast(geo_file_global)
    
    # Buffer 10 km in UTM
    roi_buff = st_buffer(roi, dist = 50000)
    
    # Clip to buffered ROI
    clipped <- crop(r, roi_buff)
    buffered_path <- file.path(temp_dir, "geo_clipped.tif")
    writeRaster(clipped, buffered_path, overwrite = TRUE)
    
    # Resample to target resolution
    ref = rast(ref_file)
    res(ref) = cellsize
    
    resampled <- resample(clipped, ref, method = "near")
    resampled_path <- file.path(temp_dir, "geo_resampled.tif")
    writeRaster(resampled, resampled_path, overwrite = TRUE)
    
    # Reclassification
    values <- values(resampled)
    values[is.na(values)] <- NA
    unique_vals <- sort(na.omit(unique(as.vector(values))))
    val_map <- setNames(seq_along(unique_vals), unique_vals)
    
    reclassified <- classify(resampled, rbind(unique_vals, val_map) |> t(), include.lowest = TRUE)
    # reclassified[is.na(reclassified)] <- -9999
    geology_class_path <- file.path(morph_folder, "geology_class.asc")
    writeRaster(reclassified, geology_class_path, overwrite = TRUE, NAflag = -9999,
                datatype = "INT4S")
    cat("Process completed: geology raster generated from global raster (GLiM).\n")
    
  }
  # Guardar LUT
  lut_df <- data.frame(
    GeoClass = names(val_map),
    Value = as.integer(val_map)
  )
  write_csv(lut_df, file.path(temp_folder, "LUT_geo_chile.csv"))
  
  # Crear archivo geology_classdefinition.txt
  n_units <- length(val_map)
  if (n_units > 25){
    cat("mHM does not support more than 25 geology classes. Current number of geo classes found: ",n_units,"\n")
  }
  param_lines <- c(
    sprintf("nGeo_Formations  %d", n_units),
    "",
    "GeoParam(i)   ClassUnit     Karstic      Description"
  )
  for (i in seq_along(unique_vals)) {
    param_lines <- c(param_lines, sprintf("%10d %13d %11d      GeoUnit-%d", i, i, 0, i))
  }
  
  param_lines <- c(
    param_lines,
    "!<-END",
    "!***********************************",
    "! NOTES",
    "!***********************************",
    "1 = Karstic",
    "0 = Non-karstic",
    "",
    "IMPORTANT ::",
    "   Ordering has to be according to the ordering in mhm_parameter.nml",
    "   (namelist: geoparameter)"
  )
  
  writeLines(param_lines, file.path(morph_folder, "geology_classdefinition.txt"))
  cat("class LUT written in:", file.path(morph_folder, "geology_classdefinition.txt"),"\n")
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

preprocess_soil_data = function(domain_path, remove_temp = FALSE){
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  
  config_path <- file.path(domain_path, "preprocess_config.json")
  config_path = "preprocess_config_windows.json"
  config <- fromJSON(config_path)
  
  roi_file = config$roi_file
  soil_files = config$soil_folder
  cellsize = config$cellsize
  
  ref_file = file.path(domain_path, config$meteo_folder,"pre.nc")
  ref_file = "meteo/pre.nc"
  morph_folder = "morph"
  morph_folder = file.path(domain_path, config$morph_folder)
  temp_folder = file.path(morph_folder, "temp");temp_folder
  dir.create(temp_folder, showWarnings = FALSE)
  
  # read ROI
  roi = read_sf(roi_file)
  box = get_extent(roi)
  # Crear buffer de 1km (3000 metros)
  roi_buff <- st_buffer(box, dist = 3000)
  # read reference grid
  ref = rast(ref_file)
  res(ref) = cellsize
  # refernce morph data
  reproject_to_mhm_morph <- function(r) {
    # r = crop(r, roi_buff)
    # r = crop(r, ref)
    r = resample(r, ref, method = "bilinear")
  }
  
  files = list.files(soil_files, full.names = TRUE)
  bulkd_files = files[grep("Bulkd", files)]
  clay_files = files[grep("Clay", files)]
  sand_files = files[grep("Sand", files)]
  
  bulkd_r = reproject_to_mhm_morph(rast(bulkd_files))
  clay_r = reproject_to_mhm_morph(rast(clay_files))
  sand_r = reproject_to_mhm_morph(rast(sand_files))
  
  
  focal_repeat = function(r, n){
    if (n != 0){
      r = terra::focal(x = r, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
      return(focal_repeat(r, n-1))
    }else{
      return(r)
    }
  }
  
  bulkd_filled = focal_repeat(bulkd_r, 10)
  # plot(bulkd_r)
  # plot(bulkd_filled)
  # writeRaster(bulkd_filled, filename = file.path(temp_folder, "bulkd_filled.tif"), overwrite = TRUE)
  
  clay_filled = focal_repeat(clay_r, 10)
  # plot(clay_r)
  # plot(clay_filled)
  # writeRaster(clay_filled, filename = file.path(temp_folder, "clay_filled.tif"), overwrite = TRUE)
  
  sand_filled = focal_repeat(sand_r, 10)
  # plot(sand_r)
  # plot(sand_filled)
  # writeRaster(sand_filled, filename = file.path(temp_folder, "sand_filled.tif"), overwrite = TRUE)
  
  depths = c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm","100-200cm")
  # df to hold unique combinations from all the layers
  data.depths <- data.frame()
  
  # imagenes por profundidad
  get_combinations <- function(r, var_names = c("Bulkd","Clay","Sand")) {
    lapply(depths, function(d){
      cat("Cargando mapas de profundidad ", d,"\n")
      index = grep(d, names(r))
      imgs = r[[index]]
      cat("Extrayendo valores\n")
      imgs.df <- values(imgs)
      colnames(imgs.df) <- var_names
      imgs.df %>% as_tibble()
    }) %>% 
      bind_rows()
  }
  
  all_soil = c(
    round(bulkd_filled,digits=2),
    round(sand_filled,digits=0),
    round(clay_filled,digits=0)
  )
  
  soil.mask = sum(!is.na(all_soil))
  soil.mask[soil.mask != 18] = NA
  all_soil = mask(all_soil, soil.mask)
  
  data.depths = get_combinations(all_soil)
  glimpse(data.depths)
  
  # contar combinaciones presentes
  comb.data <- data.depths %>% group_by_all() %>% count() %>% drop_na()
  
  # definir ID de tipo de suelo
  comb.data$comb_id <- seq(1, nrow(comb.data));comb.data
  
  # beepr::beep(sound = 5)
  
  for (i in 1:6) {
    # i=1
    d <- depths[i]
    # cargar imagenes
    cat("Cargando mapas de profundidad ", d,"\n")
    index = grep(d, names(all_soil))
    imgs = all_soil[[index]]
    # extraer valores
    cat("Extrayendo valores...\n")
    imgs.df <- values(imgs)
    colnames(imgs.df) <- c("Bulkd","Clay","Sand")
    
    # redondear texturas a enteros y Da a 1 decimal
    imgs.df = imgs.df %>% 
      as_tibble()
    #asignar el ID de suelo a cada pixel
    cat("Asignando clases de suelo...\n")
    imgs.class <- left_join(imgs.df, comb.data, by = c("Clay","Sand","Bulkd"), keep = F)
    imgs.class
    #asignar valores al raster
    v <- rast(all_soil[[1]])
    values(v) <- imgs.class$comb_id
    names(v) = paste0("soilclass_0",i)
    
    #exportar resultado
    cat("Exportando resultados...\n")
    names(v) = paste0("soil_class_horizon_0",i,".tif")
    v %>% writeRaster(str_c(temp_folder,"/soil_class_horizon_0",i,".tif"), overwrite = T)
  }
  
  # ruta de imagenes
  imgs <- list.files(temp_folder, full.names = T, pattern = "soil_class_horizon_[0-9]+.tif$") %>% rast;imgs
  
  EX.TABLE = comb.data %>% 
    mutate(ID = comb_id) %>% 
    mutate(Clay = as.integer(Clay),
           Sand = as.integer(Sand)) %>% 
    select(
      ID,
      `CLAY[%]` = Clay,
      `SAND[%]` = Sand,
      `Bd_mu[gcm-3]` = Bulkd,
      nSamples = n
    )
  
  
  nSoiltypes = nrow(comb.data);nSoiltypes
  write_lines(str_c("nSoil_Types\t" ,as.integer(nSoiltypes)) ,str_c(morph_folder, "/soil_classdefinition_iFlag_soilDB_1.txt"))
  write_delim(EX.TABLE, str_c(morph_folder, "/soil_classdefinition_iFlag_soilDB_1.txt"), 
              append = TRUE, 
              col_names = TRUE, delim = "\t")
  
  # transform to ASCII
  rr <- rast(imgs)
  names(imgs)
  for (i in 1:6) {
    # i = 1
    print(i)
    img <- imgs[[i]]
    rr <- c(rr, img)
  }
  is.na(rr) %>% plot
  
  # beepr::beep(sound = 5)
  names(rr) = paste0("soil_class_horizon_0",c(1,2,3,4,5,6),".tif")
  writeRaster(rr,
              paste0(morph_folder,"/soil_class_horizon_0",c(1,2,3,4,5,6),".asc"),
              overwrite = T, datatype = "INT4S", NAflag = -9999)
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

create_roi_mask = function(config_path, remove_temp = FALSE, basin.mask = FALSE){
  library(terra)
  library(jsonlite)
  
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  morph_folder = file.path(domain_path, config$morph_folder)
  lc_folder = file.path(domain_path, config$lc_folder)
  lai_folder = file.path(domain_path, config$lai_folder)
  
  roi_file = config$roi_file
  dem_file = file.path(morph_folder, "dem.asc")
  geo_file = file.path(morph_folder, "geology_class.asc")
  lc_file = file.path(lc_folder, "landcover.asc")
  lai_file = file.path(lai_folder, "lai.nc")
  
  dem = rast(dem_file)
  geo = rast(geo_file)
  lc = rast(lc_file)
  lai = rast(lai_file)
  
  dem.m = !is.na(dem)
  geo.m = !is.na(geo)
  lc.m = !is.na(lc)
  lai.m = !is.na(lai[[1]])
  
  roi_mask = dem.m + geo.m + lc.m + lai.m
  max.val = max(values(roi_mask), na.rm = TRUE)
  roi_mask[roi_mask != max.val] = NA
  # plot(roi_mask)
  writeRaster(roi_mask, file.path(morph_folder, "roi_mask.tif"), overwrite = TRUE)
  
  # apply mask to dem
  roi_mask = rast(file.path(morph_folder, "roi_mask.tif"))

  aspect_file = file.path(morph_folder, "aspect.asc")
  slope_file = file.path(morph_folder, "slope.asc")
  fdir_file = file.path(morph_folder, "fdir.asc")
  facc_file = file.path(morph_folder, "facc.asc")
  
  dem = dem_file %>% rast %>% mask(roi_mask)
  aspect = aspect_file %>% rast %>% mask(roi_mask)
  slope = slope_file %>% rast %>% mask(roi_mask)
  fdir = fdir_file %>% rast %>% mask(roi_mask)
  facc = facc_file %>% rast %>% mask(roi_mask)
  
  if (basin.mask){
    roi_mask = vect(roi_file)
    dem = mask(dem, roi_mask)
    aspect = mask(aspect, roi_mask)
    slope = mask(slope, roi_mask)
    fdir = mask(fdir, roi_mask)
    facc = mask(facc, roi_mask)
  }
  writeRaster(dem, file.path(morph_folder, "dem.asc"), overwrite = TRUE, NAflag = -9999)
  writeRaster(aspect, file.path(morph_folder, "aspect.asc"), overwrite = TRUE, NAflag = -9999)
  writeRaster(slope, file.path(morph_folder, "slope.asc"), overwrite = TRUE, NAflag = -9999)
  writeRaster(fdir, file.path(morph_folder, "fdir.asc"), overwrite = TRUE, NAflag = -9999, datatype = "INT4S")
  writeRaster(facc, file.path(morph_folder, "facc.asc"), overwrite = TRUE, NAflag = -9999, datatype = "INT4S")
  
  cat("Mask applied to dem, aspect, slope, fdir and facc. Files overwritten")
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
  
}


preprocess_streamflow_data = function(domain_path, remove_temp = FALSE){
  library(jsonlite)
  library(lubridate)
  library(tidyverse)
  
  # === Load configuration ===
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  
  gauges_file = config$fluv_station_file
  roi_file = config$roi_file
  gauge_folder <- file.path(domain_path, config$gauge_folder)
  dir.create(gauge_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Find gauges inside roi
  roi = read_sf(roi_file)
  gauges = read_sf(gauges_file) %>% st_transform(crs = st_crs(roi))
  
  gauge_list = as.character(st_intersection(roi, gauges)$ID)
  cat(length(gauge_list),"Gauge stations found inside ROI: ", paste(gauge_list, collapse = ", "))
  
  dir.create(gauge_folder, showWarnings = FALSE)
  
  # === Read CSV assuming the first column is date ===
  df <- read_csv(streamflow_data_file)
  colnames(df)[1] <- "Date"
  
  for (i in seq_along(gauge_list)) {
    gauge_id <- gauge_list[[i]]
    if (!(gauge_id %in% colnames(df))) {
      cat(sprintf("Gauge ID %s not found in data.\n", gauge_id))
      next
    }
    
    series_valid <- df %>% select(Date, !!gauge_id) %>% drop_na()
    series_valid
    
    if (nrow(series_valid) == 0) {
      cat(sprintf("No valid data for gauge %s. Skipping.\n", gauge_id))
      next
    }
    
    start_date <- min(series_valid$Date)
    end_date <- max(series_valid$Date)
    
    
    full_range <- data.frame(Date = seq(start_date, end_date, by = "1 day"))
    filled_series <- full_range %>%
      left_join(series_valid, by = "Date") %>%
      mutate(Value = ifelse(is.na(.data[[gauge_id]]), -9999, .data[[gauge_id]])) %>%
      mutate(Hour = 0, Minute = 0)
    
    # Header lines
    output_lines <- c(
      sprintf("%s Gauge %d (daily discharge)", gauge_id, i),
      "nodata   -9999",
      "n       1       measurements per day [1, 1440]",
      sprintf("start  %s", format(start_date, "%Y %m %d %H %M")),
      sprintf("end    %s", format(end_date, "%Y %m %d %H %M"))
    )
    
    # Data lines
    data_lines <- filled_series %>%
      mutate(Line = sprintf("%4d  %02d  %02d  %02d  %02d   %8.3f",
                            year(Date), month(Date), day(Date), Hour, Minute, Value)) %>%
      pull(Line)
    
    output_lines <- c(output_lines, data_lines)
    
    # Write output
    output_file <- file.path(gauge_folder, paste0(gauge_id, ".day"))
    writeLines(output_lines, output_file)
    
    cat(sprintf("Wrote file: %s\n", output_file))
  }
}

create_idgauges = function(domain_path, remove_temp = FALSE){
  library(sf)
  library(terra)
  library(readr)
  library(jsonlite)
  
  # === Load config ===
  config_path <- file.path(domain_path, "preprocess_config.json")
  config <- fromJSON(config_path)
  
  gauge_file <- config$fluv_station_file
  gauge_folder = file.path(domain_path, config$gauge_folder)
  morph_folder <- file.path(domain_path, config$morph_folder)
  meteo_folder = file.path(domain_path, config$meteo_folder)
  nodata_value <- -9999
  cellsize = config$cellsize
  
  gauge_list = list.files(gauge_folder, pattern = ".day") %>% 
    str_sub(end = -5) %>% 
    as.numeric()
  
  gauges = read_sf(gauge_file) %>% 
    filter(ID %in% gauge_list) %>% 
    st_transform(4326)
  
  # === Load reference raster ===
  ref_path <- file.path(meteo_folder, "pre.nc")
  ref <- rast(ref_path)
  res(ref) = cellsize
  
  # === Rasterize station ID ===
  gdf_vect <- vect(gauges)
  id_raster <- rasterize(gdf_vect, ref, field = "ID")
  
  # === Write raster to ASCII ===
  output_path <- file.path(morph_folder, "idgauges.asc")
  writeRaster(id_raster, output_path, overwrite = TRUE, NAflag = nodata_value, datatype = "INT4S")
  
  cat("idgauges.asc saved at",output_path,"\n")
}