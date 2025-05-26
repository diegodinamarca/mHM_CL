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
                "jsonlite") # list of packages to load
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
  source("scripts/R/utils.R")
  
}

# leer archivo de configuracion json con rutas de archivos
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")

# Paths from config
raster_path <- config$land_cover_file
roi_path <- config$roi_file
morph_folder <- config$morph_folder
lc_folder <- config$lc_folder
ref_file <- file.path(config$meteo_folder, "pre.nc")
cellsize <- config$cellsize

temp_dir <- file.path(lc_folder, "temp")
if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

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
clipped_path <- file.path(temp_dir, "landcover_clipped.tif")
writeRaster(lc_clipped, clipped_path, overwrite = TRUE)

# Step 2: Reclassify using dict
lc_reclass = classify(lc_clipped, rcl = reclass_dict, include.lowest = T, right = T)
lc_reclass[lc_reclass == 0] = NA
plot(lc_reclass)

# Save reclassified raster
reclass_path <- file.path(temp_dir, "landcover_reclassified.tif")
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
writeRaster(lc_rs, final_lc_path, overwrite = TRUE, NAflag = -9999, datatype = "INT4S")

cat("✅ Final landcover.asc saved to:", final_lc_path, "\n")
