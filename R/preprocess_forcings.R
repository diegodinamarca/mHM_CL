preprocess_climate_data <- function(config_path, remove_temp = TRUE) {
  library(jsonlite)
  library(terra)
  library(sf)
  library(magrittr)
  
  # Load config
  config <- fromJSON(config_path)
  
  roi_path <- config$roi_file
  date1 <- config$start_date
  date2 <- config$end_date
  variables <- config$variables_clim
  cellsize_clim <- config$cellsize_clim
  cellunit <- config$cellunit
  clim_folder <- config$meteo_folder
  header_folder <- config$header_folder
  
  # Prepare folders
  dir.create(clim_folder, showWarnings = FALSE, recursive = TRUE)
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
  
  for (VAR in names(variables)) {
    paths <- variables[[VAR]]
    input_dir <- paths$input_dir
    if (is.null(input_dir)) next
    
    input_files <- list.files(input_dir, pattern = "\\.nc$", full.names = TRUE)
    cat("🔄 Processing variable:", VAR, "\n")
    cat("  📂 Found", length(input_files), "files in", input_dir, "\n")
    
    r <- rast(input_files, subds = VAR)
    r.cropped <- crop(r, buffered_box)
    r.tf <- subset(r.cropped, time(r.cropped) >= as.Date(date1) & time(r.cropped) <= as.Date(date2))
    varnames(r.tf) <- "pre"
    ref <- rast(r.tf)
    res(ref) <- cellsize_clim
    r.rs <- resample(r.tf, ref, method = "bilinear")
    
    resampled_path <- file.path(temp_folder, paste0(mhm_varnames[[VAR]], "_resampled.nc"))
    writeCDF(r.rs, filename = resampled_path, overwrite = TRUE)
    
    buffered_3km <- st_buffer(box.pol, dist = 3000)
    r.rs.crop <- crop(r.rs, buffered_3km)
    final_path <- file.path(clim_folder, paste0(mhm_varnames[[VAR]], ".nc"))
    writeCDF(r.rs.crop, filename = final_path, overwrite = TRUE, prec = "double", missval = -9999)
    
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
}

preprocess_LAI_data = function(config_path, remove_temp = FALSE){
  library(terra)
  library(sf)
  library(jsonlite)
  library(tidyverse)
  # source("scripts/R/utils.R")
  
  # === Load config ===
  config <- fromJSON(config_path)
  
  lai_path <- config$lai_file
  roi_path <- config$roi_file
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
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}



preprocess_dem_data<- function(config_path, remove_temp = FALSE) {
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  library(rdwplus)

  # === Load config ===
  config <- jsonlite::fromJSON(config_path)
  outdir <- config$morph_folder
  temp_folder <- file.path(outdir, "temp")
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
  
  ref_file <- file.path(config$meteo_folder, "pre.nc")
  dem_file <- config$dem_file
  roi_file <- config$roi_file
  hydronet_file <- config$hydro_network_file
  cellsize <- config$cellsize
  header_folder = config$header_folder
  
  # === Read inputs ===
  dem <- rast(dem_file)
  roi <- st_read(roi_file, quiet = TRUE)
  hydronet <- st_read(hydronet_file, quiet = TRUE) %>% st_transform(4326)
  
  box = get_extent(roi)
  roi_buff <- st_buffer(box, dist = 3000)
  hydronet <- st_crop(hydronet, roi_buff)
  hydronet_crop_file <- file.path(temp_folder, "hydronet.shp")
  
  dem_crop <- crop(dem, vect(roi_buff))
  dem_crop_file <- file.path(temp_folder, "dem_clipped.tif")
  writeRaster(dem_crop, filename = dem_crop_file, overwrite = TRUE)
  
  ref <- rast(ref_file)[[1]]
  res(ref) <- cellsize
  dem_rs <- resample(dem_crop, ref, method = "bilinear")
  dem_rs_file <- file.path(temp_folder, "dem_resampled.tif")
  writeRaster(dem_rs, filename = dem_rs_file, overwrite = TRUE, datatype = "INT4S")
  
  # === GRASS setup ===
  my_grass <- "/Applications/GRASS-8.4.app/Contents/Resources"
  Sys.setenv(PATH = paste(file.path(my_grass, "bin"),
                          file.path(my_grass, "scripts"),
                          Sys.getenv("PATH"), sep = ":"))
  initGRASS(my_grass, mapset = "PERMANENT", override = TRUE)
  check_running()
  set_envir(dem_rs_file)
  
  vector_to_mapset(hydronet_crop_file, overwrite = TRUE)
  vibe_check()
  
  rasterise_stream("hydronet", file.path(temp_folder,"hydronet.tif"), TRUE)
  reclassify_streams("hydronet.tif", "hydronet01.tif", overwrite = TRUE)
  
  burn_in("dem_resampled.tif", "hydronet01.tif", "burndem.tif", overwrite = TRUE)
  
  nsinks = function(sink_file){
    sinks = rast(sink_file)
    cells = ncell(sink)
    # plot(sinks)
    nsink = sinks[sinks != -1] %>% nrow
  }
  
  fill_sinks("burndem.tif", "filldem.tif", "fd1.tif","sinks.tif", overwrite = TRUE)
  sink_file = file.path(temp_folder, "sinks.tif")
  retrieve_raster("sinks.tif", sink_file, overwrite = TRUE)
  # print(nsinks(sink_file))
  for (i in 1:1) {
    # print(i)
    fill_sinks(
      dem = "filldem.tif", 
      out_dem = "filldem.tif", 
      out_fd = "fd1.tif", # flow-direction grid that cannot be used in subsequent steps but is required by GRASS
      overwrite = T,
      out_sinks = "sinks.tif"
    )
    sink_file = file.path(temp_folder, "sinks.tif")
    retrieve_raster("sinks.tif", sink_file, overwrite = TRUE)
    nsinks_out = nsinks(sink_file)
    nsinks_cell = ncell(rast(sink_file))
    
    if (nsinks_out == nsinks_cell){
      cat("\n All sinks filled")
    }
  }
  derive_flow("filldem.tif", "flowdir.tif", "flowacc.tif", overwrite = TRUE)
  
  # plot_GRASS("filldem.tif", colours = topo.colors(15))
  # plot_GRASS("flowdir.tif", colours = topo.colors(6))
  # plot_GRASS("flowacc.tif", colours = topo.colors(6))
  
  retrieve_raster("filldem.tif", file.path(temp_folder, "filldem.tif"), overwrite = TRUE)
  retrieve_raster("flowacc.tif", file.path(temp_folder, "flowacc.tif"), overwrite = TRUE)
  retrieve_raster("flowdir.tif", file.path(temp_folder, "flowdir.tif"), overwrite = TRUE)
  
  # === Flowdir reclassification ===
  flowdir <- abs(rast(file.path(temp_folder, "flowdir.tif")))
  mat.reclass <- rbind(c(8,1), c(7,2), c(6,4), c(5,8),
                       c(4,16), c(3,32), c(2,64), c(1,128))
  flowdir_class <- terra::classify(flowdir, rcl = as.matrix(mat.reclass))
  writeRaster(flowdir_class, file.path(temp_folder, "flowdir_reclass.tif"), overwrite = TRUE)
  
  # === Terrain variables ===
  dem_filled <- rast(file.path(temp_folder, "filldem.tif"))
  slope <- terrain(dem_filled, v = "slope", unit = "degrees")
  aspect <- terrain(dem_filled, v = "aspect", unit = "degrees")
  writeRaster(slope, file.path(temp_folder, "slope.tif"), overwrite = TRUE)
  writeRaster(aspect, file.path(temp_folder, "aspect.tif"), overwrite = TRUE)
  
  # === Resample and export to mHM ===
  facc <- rast(file.path(temp_folder, "flowacc.tif"))
  
  dem_rs <- resample(dem_filled, ref, method = "bilinear")
  slope_rs <- resample(slope, ref, method = "bilinear")
  aspect_rs <- resample(aspect, ref, method = "bilinear")
  fdir_rs <- resample(flowdir_class, ref, method = "near")
  facc_rs <- resample(facc, ref, method = "bilinear")
  
  writeRaster(dem_rs, file.path(outdir, "dem.asc"), overwrite = TRUE, datatype = "FLT8S", NAflag = -9999)
  writeRaster(slope_rs, file.path(outdir, "slope.asc"), overwrite = TRUE, datatype = "FLT8S", NAflag = -9999)
  writeRaster(aspect_rs, file.path(outdir, "aspect.asc"), overwrite = TRUE, datatype = "FLT8S", NAflag = -9999)
  writeRaster(fdir_rs, file.path(outdir, "fdir.asc"), overwrite = TRUE, datatype = "INT4S", NAflag = -9999)
  writeRaster(facc_rs, file.path(outdir, "facc.asc"), overwrite = TRUE, datatype = "FLT8S", NAflag = -9999)
  
  cat("PROCESAMIENTO COMPLETO. Archivos guardados en:", outdir, "\n")
  
  dem_path = file.path(config$morph_folder, "dem.asc")
  extract_asc_header(dem_path, header_folder)
  
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

preprocess_lc_data = function(config_path, remove_temp = FALSE){
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  # leer archivo de configuracion json con rutas de archivos
  config <- fromJSON(config_path)
  
  # Paths from config
  raster_path <- config$land_cover_file
  roi_path <- config$roi_file
  morph_folder <- config$morph_folder
  lc_folder <- config$lc_folder
  ref_file <- file.path(config$meteo_folder, "pre.nc")
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


preprocess_geo_data <- function(config_path, remove_temp = FALSE) {
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  # leer archivo de configuracion json con rutas de archivos
  config <- fromJSON(config_path)
  
  geo_file <- config$geo_file
  roi_file <- config$roi_file
  morph_folder <- config$morph_folder
  ref_file = file.path(config$meteo_folder,"pre.nc")
  cellsize = config$cellsize
  
  dir.create(morph_folder, showWarnings = FALSE, recursive = TRUE)
  temp_folder <- file.path(morph_folder, "temp")
  dir.create(temp_folder, showWarnings = FALSE, recursive = TRUE)
  
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
  
  # Guardar LUT
  lut_df <- data.frame(
    GeoClass = names(val_map),
    Value = as.integer(val_map)
  )
  write_csv(lut_df, file.path(temp_folder, "LUT_geo_chile.csv"))
  
  # Crear archivo geology_classdefinition.txt
  n_units <- length(val_map)
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
  
  cat("Process completed: geology raster generated from shapefile.\n")
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

preprocess_soil_data = function(config_path, remove_temp = FALSE){
  library(tidyverse)
  library(terra)
  library(sf)
  library(jsonlite)
  config <- fromJSON(config_path)
  
  roi_file = config$roi_file
  soil_files = config$soil_folder
  ref_file = file.path(config$meteo_folder,"pre.nc")
  cellsize = config$cellsize
  outdir = config$morph_folder
  temp_folder = str_c(outdir, "/temp");temp_folder
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
    r = crop(r, roi_buff)
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
  writeRaster(bulkd_filled, filename = paste0(temp_folder, "/bulkd_filled.tif"), overwrite = TRUE)
  
  clay_filled = focal_repeat(clay_r, 10)
  # plot(clay_r)
  # plot(clay_filled)
  writeRaster(clay_filled, filename = paste0(temp_folder, "/clay_filled.tif"), overwrite = TRUE)
  
  sand_filled = focal_repeat(sand_r, 10)
  # plot(sand_r)
  # plot(sand_filled)
  writeRaster(sand_filled, filename = paste0(temp_folder, "/sand_filled.tif"), overwrite = TRUE)
  
  depths = c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm","100-200cm")
  # df to hold unique combinations from all the layers
  data.depths <- data.frame()
  
  # imagenes por profundidad
  get_unique_combinations <- function(r, var_names = c("Bulkd","Clay","Sand")) {
    lapply(depths, function(d){
      cat("Cargando mapas de profundidad ", d,"\n")
      index = grep(d, names(r))
      imgs = r[[index]]
      cat("Extrayendo valores\n")
      imgs.df <- values(imgs)
      colnames(imgs.df) <- var_names
      imgs.df <- imgs.df %>% as_tibble %>% 
        mutate(Clay = round(Clay, digits = 0),
               Sand = round(Sand, digits = 0),
               Bulkd = round(Bulkd, digits = 2))
    }) %>% 
      bind_rows() %>% 
      group_by_all() %>%
      slice(1) %>%
      ungroup()
    
  }
  
  all_soil = c(
    bulkd_filled,
    sand_filled,
    clay_filled
  )

  data.depths = get_unique_combinations(all_soil)
  glimpse(data.depths)
  
  # contar combinaciones presentes
  comb.data <- data.depths %>% group_by_all() %>% count
  
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
      as_tibble() %>% 
      mutate(Clay = as.integer(round(Clay, digits = 0)),
             Sand = as.integer(round(Sand, digits = 0)),
             Bulkd = round(Bulkd, digits = 2))
    #asignar el ID de suelo a cada pixel
    cat("Asignando clases de suelo...\n")
    imgs.class <- left_join(imgs.df, comb.data, by = c("Clay","Sand","Bulkd"), keep = F)
    #asignar valores al raster
    v <- rast(all_soil[[1]])
    values(v) <- imgs.class$comb_id
    names(v) = paste0("soilclass_0",i)
    
    #exportar resultado
    cat("Exportando resultados...\n")
    names(v) = paste0("soil_class_horizon_0",i,".tif")
    v %>% writeRaster(str_c(temp_folder,"/soil_class_horizon_0",i,".tif"), overwrite = T)
    gc()
  }
  
  # exportar tabla de clasificaciones
  comb.data %>% write_csv(paste0(temp_folder,"/soilclass_withNA.csv"))
  
  # ruta de imagenes
  imgs <- list.files(temp_folder, full.names = T, pattern = "soil_class_horizon_[0-9]+.tif$") %>% rast;imgs
  files <- list.files(temp_folder, full.names = T, pattern = "soil_class_horizon_[0-9]+.tif$");files
  
  # Remover clases con valores NA
  rm.class <- comb.data %>% ungroup() %>% filter(is.na(Sand) | is.na(Clay) | is.na(Bulkd)) %>%
    select(comb_id) %>% unique
  rm.class <- rm.class$comb_id
  rm.class
  
  # exportar clases sin valores NA
  comb.data = comb.data %>% 
    filter(!(comb_id %in% rm.class)) %>%
    write_csv(paste0(temp_folder,"/soilclasses.csv"))
  
  EX.TABLE = comb.data %>% 
    mutate(ID = cur_group_id()) %>% 
    mutate(Clay = as.integer(Clay),
           Sand = as.integer(Sand)) %>% 
    select(
      ID,
      `CLAY[%]` = Clay,
      `SAND[%]` = Sand,
      `Bd_mu[gcm-3]` = Bulkd,
      nSamples = comb_id
    )
  
  
  nSoiltypes = nrow(comb.data);nSoiltypes
  write_lines(str_c("nSoil_Types\t" ,as.integer(nSoiltypes)) ,str_c(outdir, "/soil_classdefinition_iFlag_soilDB_1.txt"))
  write_delim(EX.TABLE, str_c(outdir, "/soil_classdefinition_iFlag_soilDB_1.txt"), 
              append = TRUE, 
              col_names = TRUE, delim = "\t")
  
  # transform to ASCII
  rr <- rast(imgs)
  for (i in 1:6) {
    # i = 1
    print(i)
    img <- imgs[[i]]
    img[img %in% rm.class] <- NA
    rr <- c(rr, img)
  }

  names(rr) = paste0("soil_class_horizon_0",c(1,2,3,4,5,6),".tif")
  writeRaster(rr, 
              paste0(temp_folder,"/soil_class_horizon_0",c(1,2,3,4,5,6),"_NA_removed.tif"),
              overwrite = T)
  writeRaster(rr,
              paste0(outdir,"/soil_class_horizon_0",c(1,2,3,4,5,6),".asc"),
              overwrite = T, datatype = "INT4S", NAflag = -9999)
  if (remove_temp) {
    unlink(temp_folder, recursive = TRUE)
    cat("🗑️ Removed temp folder:", temp_folder, "\n")
  }
}

create_roi_mask = function(config_path, remove_temp = FALSE, basin.mask = FALSE){
  config <- fromJSON(config_path)
  morph_folder = config$morph_folder
  lc_folder = config$lc_folder
  lai_folder = config$lai_folder
  
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


preprocess_streamflow_data = function(config_path, remove_temp = FALSE){
  library(jsonlite)
  library(lubridate)
  library(tidyverse)
  
  # === Load configuration ===
  config <- fromJSON(config_path)
  
  streamflow_data_file <- config$streamflow_data_file
  gauge_list <- config$gauge_list
  gauge_folder <- config$gauge_folder
  
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
    output_file <- file.path(gauge_folder, paste0(gauge_id, ".txt"))
    writeLines(output_lines, output_file)
    
    cat(sprintf("Wrote file: %s\n", output_file))
  }
}

create_idgauges = function(config_path, remove_temp = FALSE){
  library(sf)
  library(terra)
  library(readr)
  library(jsonlite)
  
  # === Load config ===
  config <- fromJSON(config_path)
  
  csv_path <- config$fluv_station_file
  morph_folder <- config$morph_folder
  meteo_folder = config$meteo_folder
  nodata_value <- -9999
  cellsize = config$cellsize
  
  # === Load reference raster ===
  ref_path <- file.path(meteo_folder, "pre.nc")
  ref <- rast(ref_path)
  res(ref) = cellsize
  
  # === Read CSV and convert to sf ===
  gdf = read_csv(csv_path) %>% st_as_sf(coords = c("LON","LAT"), crs = 4326)
  
  # === Rasterize station ID ===
  gdf_vect <- vect(gdf)
  id_raster <- rasterize(gdf_vect, ref, field = "ID")
  
  # === Write raster to ASCII ===
  output_path <- file.path(morph_folder, "idgauges.asc")
  writeRaster(id_raster, output_path, overwrite = TRUE, NAflag = nodata_value, datatype = "INT4S")
  
  cat("idgauges.asc saved at",output_path,"\n")
}