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
ref_file = "/Users/mhm/Desktop/FONDECYT_CAMILA/7339001/morph/dem.asc"

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
hydronet = st_crop(hydronet, roi_buff)
hydronet_crop_file = file.path(temp_folder, "hydronet.shp")
# write_sf(hydronet, hydronet_crop_file)

# Cortar DEM al buffer
dem_crop <- crop(dem, vect(roi_buff))
dem_crop_file = file.path(temp_folder, "dem_clipped.tif")
writeRaster(dem_crop, filename = dem_crop_file, overwrite = TRUE)

# resample DEM to L0 resolution
# load reference file
ref <- rast(ref_file)[[1]]
res(ref) = cellsize
ref

dem_rs = resample(dem_crop, ref, method = "bilinear")
dem_rs
# dem_rs = dem_crop
dem_rs_file = file.path(temp_folder, "dem_resampled.tif")
writeRaster(dem_rs, filename = dem_rs_file, overwrite = TRUE, datatype = "INT4S")

# Rellenar deplot()# Rellenar depresiones (fill sinks)
my_grass = "/Applications/GRASS-8.4.app/Contents/Resources"

# Add bin and scripts to PATH
Sys.setenv(PATH = paste(
  file.path(my_grass, "bin"),
  file.path(my_grass, "scripts"),
  Sys.getenv("PATH"),
  sep = ":"
))

initGRASS(my_grass, mapset = "PERMANENT", override = TRUE)
check_running()
set_envir(dem_rs_file)

# cargar hydronet a grass
vector_to_mapset(
  hydronet_crop_file, # vector of inputs (give file paths)
  overwrite = TRUE
)
# ver variable cargadas
vibe_check()

rasterise_stream(
  streams = "hydronet", # input of vector stream lines
  out = file.path(temp_folder,"hydronet.tif"), # name of output raster
  TRUE
)
reclassify_streams(
  "hydronet.tif", # input raster
  "hydronet01.tif", # output raster
  overwrite = TRUE
)

# burn streams
burn_in(
  dem = "dem_resampled.tif", 
  stream = "hydronet01.tif", # should be in 0-1 format, 0 for non-stream cells, 1 for stream cells,
  out = "burndem.tif", 
  overwrite = TRUE
)

nsinks = function(sink_file){
  sinks = rast(sink_file)
  cells = ncell(sink)
  plot(sinks)
  nsink = sinks[sinks != -1] %>% nrow
}

fill_sinks(
  dem = "burndem.tif", 
  out_dem = "filldem.tif", 
  out_fd = "fd1.tif", # flow-direction grid that cannot be used in subsequent steps but is required by GRASS
  overwrite = T,
  out_sinks = "sinks.tif"
)
sink_file = file.path(temp_folder, "sinks.tif")
retrieve_raster("sinks.tif", sink_file, overwrite = TRUE)
print(nsinks(sink_file))
for (i in 1:5) {
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


# Derive flow direction and accumulation rasters
derive_flow(
  dem = "filldem.tif", 
  flow_dir = "flowdir.tif", # uses d8 by default
  flow_acc = "flowacc.tif", 
  overwrite = TRUE
)

# Filled DEM 
plot_GRASS("filldem.tif", colours = topo.colors(15))
# Flow Direction
plot_GRASS("flowdir.tif", colours = topo.colors(6)) # grid has fewer than 8 directions
# Flow Accumulation
plot_GRASS("flowacc.tif", topo.colors(6))

retrieve_raster("filldem.tif", file.path(temp_folder, "filldem.tif"), overwrite = TRUE)
retrieve_raster("flowacc.tif", file.path(temp_folder, "flowacc.tif"), overwrite = TRUE)
retrieve_raster("flowdir.tif", file.path(temp_folder, "flowdir.tif"), overwrite = TRUE)

# reclassify flow direction
flowdir_file = file.path(temp_folder, "flowdir.tif")
flowdir = abs(rast(flowdir_file))
mat.reclass = rbind(
  c(8,1),
  c(7,2),
  c(6,4),
  c(5,8),
  c(4,16),
  c(3,32),
  c(2,64),
  c(1,128)
)
flowdir_class = terra::classify(flowdir, rcl = as.matrix(mat.reclass))
writeRaster(flowdir_class, file.path(temp_folder, "flowdir_reclass.tif"), overwrite = TRUE)

# Calcular slope, aspect, flow direction, flow accumulation
dem_filled_file = file.path(temp_folder, "filldem.tif")
dem_filled = rast(dem_filled_file)
slope <- terrain(dem_filled, v = "slope", unit = "degrees")
aspect <- terrain(dem_filled, v = "aspect", unit = "degrees")

# flowdir = terrain(dem_filled, v = "flowdir")
# facc = flowAccumulation(flowdir)
# plot(flowdir)
# plot(facc)

# Escribir las capas
writeRaster(slope, filename = file.path(temp_folder, "slope.tif"), overwrite = TRUE)
writeRaster(aspect, filename = file.path(temp_folder, "aspect.tif"), overwrite = TRUE)

# read flow accumulation
facc = rast(file.path(temp_folder, "flowacc.tif"))

# === OPCIONAL: RESAMPLE TO REFERENCE ===========================

# resample
dem_rs = resample(dem_filled, ref, method = "bilinear")
slope_rs <- resample(slope, ref, method = "bilinear")
aspect_rs <- resample(aspect, ref, method = "bilinear")
fdir_rs <- resample(flowdir_class, ref, method = "near")
facc_rs <- resample(facc, ref, method = "bilinear")


r.mask = rast(ref_file)
dem_rs = mask(dem_rs, r.mask)
slope_rs = mask(slope_rs, r.mask)
aspect_rs = mask(aspect_rs, r.mask)
fdir_rs = mask(fdir_rs, r.mask)
facc_rs = mask(facc_rs, r.mask)

terrain(x = dem_rs, "flowdir")
flowAccumulation()
# dem_rs = ddem_rs# dem_rs = dem_filled
# slope_rs <- slope
# aspect_rs <- aspect
# fdir_rs <- flowdir_class
# facc_rs <- facc

# Guardar versiones remuestreadas
writeRaster(dem_rs, filename = file.path(outdir, "dem.asc"), overwrite = TRUE,
            datatype = "FLT8S", NAflag = -9999)
writeRaster(slope_rs, filename = file.path(outdir, "slope.asc"), overwrite = TRUE,
            datatype = "FLT8S", NAflag = -9999)
writeRaster(aspect_rs, filename = file.path(outdir, "aspect.asc"), overwrite = TRUE,
            datatype = "FLT8S", NAflag = -9999)
writeRaster(fdir_rs, filename = file.path(outdir, "fdir.asc"), overwrite = TRUE,
            datatype = "INT4S", NAflag = -9999)
writeRaster(facc_rs, filename = file.path(outdir, "facc.asc"), overwrite = TRUE,
            datatype = "INT4S", NAflag = -9999)

cat("\n PROCESAMIENTO COMPLETO. Archivos guardados en:", outdir, "\n")
