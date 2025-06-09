# procesar DEM para mHM
{
  # HEADER --------------------------------------------
  #
  # Author: Diego Dinamarca
  # Email:  ddinamarcamuller@gmail.com
  # 
  # Date:
  #
  # Script Name: 
  #
  # Script Description: preparar datos de soilgrids para modelacion
  #
  # Notes:
  #
  #
  # INSTALL PACKAGES & LOAD LIBRARIES -----------------
  cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
  packages <- c("tidyverse",
                "terra", 
                "sf",
                "here",
                "jsonlite",
                "rdwplus") # list of packages to load
  n_packages <- length(packages) # count how many packages are required
  
  new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
  
  # install missing packages
  if(length(new.pkg)){
    install.packages(new.pkg)
  }
  
  # load all requried libraries
  for(n in 1:n_packages){
    cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
    lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
    eval(parse(text = lib_load)) # evaluate the string to load the library
  }
  # SET WORKING DIRECTORY -----------------------------
  cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
  wd <- here::here()
  setwd(wd)
  cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
  
  # SET OPTIONS ---------------------------------------
  cat("SETTING OPTIONS... \n\n", sep = "")
  options(scipen = 999) # turns off scientific notation
  options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
  
  
  # CONFLICTS ---------------------------------------------------------------
  conflicted::conflict_prefer("select", "dplyr")
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("extract", "terra")
  conflicted::conflict_prefer("last", "dplyr")
  # LOAD FUNCTIONS ------------------------------------
}

# leer archivo de configuracion json con rutas de archivos
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")

# output directory
outdir = config$morph_folder
# temporal files
temp_folder = file.path(outdir, "temp");temp_folder
dir.create(temp_folder)
# leer raster de referencia pre.nc
ref_file = file.path(config$meteo_folder, "pre.nc")
# ref_file = "/Users/mhm/Desktop/FONDECYT_CAMILA/7339001/morph/dem.asc"

# leer dem
dem_file = config$dem_file
# leer shp
roi_file = config$roi_file
# leer red hidrografica
hydronet_file = config$hydro_network_file
# target cellsize
cellsize = config$cellsize

# === LEER DATOS ================================================
dem <- rast(dem_file)
roi <- st_read(roi_file)

hydronet = st_read(hydronet_file) %>% st_transform(4326)
box = st_bbox(roi)
box.m = as.matrix(rbind(
  c(box[[1]], box[[2]]),
  c(box[[1]], box[[4]]),
  c(box[[3]], box[[4]]),
  c(box[[3]], box[[2]]),
  c(box[[1]], box[[2]])
))
box.pol = st_polygon(list(box.m)) %>% st_sfc(crs = 4326)
# box.pol %>% write_sf(file.path(out_temp, "bounding_box.geojson"))

# === PROCESAMIENTO =============================================
# Reproyectar ROI a CRS del DEM si es necesario
if (st_crs(roi)$epsg != crs(dem, proj=TRUE)) {
  roi <- st_transform(roi, crs(dem))
}

# Crear buffer de 1km (1000 metros)
roi_buff <- st_buffer(box.pol, dist = 3000)
# write_sf(roi_buff, file.path(out_temp, "roi_buff_3km.geojson"))

# cortar red hidrografica
# hydronet = st_crop(hydronet, roi_buff)
# hydronet_crop_file = file.path(temp_folder, "hydronet.shp")
# write_sf(hydronet, hydronet_crop_file)

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
writeRaster(rast(filled_dem), filename = file.path(outdir, "dem.asc"), overwrite = TRUE,
            datatype = "FLT8S", NAflag = -9999)
writeRaster(rast(slope), filename = file.path(outdir, "slope.asc"), overwrite = TRUE,
            datatype = "FLT8S", NAflag = -9999)
writeRaster(rast(aspect), filename = file.path(outdir, "aspect.asc"), overwrite = TRUE,
            datatype = "FLT8S", NAflag = -9999)
writeRaster(rast(fdir), filename = file.path(outdir, "fdir.asc"), overwrite = TRUE,
            datatype = "INT4S", NAflag = -9999)
writeRaster(rast(facc), filename = file.path(outdir, "facc.asc"), overwrite = TRUE,
            datatype = "INT4S", NAflag = -9999)

cat("\n PROCESAMIENTO COMPLETO. Archivos guardados en:", outdir, "\n")
