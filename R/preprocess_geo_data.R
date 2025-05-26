library(sf)
library(terra)
library(dplyr)
library(readr)
library(jsonlite)
source("scripts/R/utils.R")

process_vector_to_raster <- function(geo_file, roi_file, target_res, output_folder, 
                                     ref_grid, cellsize, clean_temp = TRUE) {
  
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  temp_dir <- file.path(output_folder, "temp")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Leer ROI
  roi <- st_read(roi_file, quiet = TRUE)
  box = get_extent(roi)
  bbox = st_buffer(box, dist = 3000)
  
  # Leer shapefile geológico
  geo <- st_read(geo_file, quiet = TRUE)

  # Recortar geología a la extensión de la ROI
  sf_use_s2(FALSE)
  geo_clipped <- st_crop(geo, bbox)

  # Crear LUT para 'cd_geol' (string → número)
  unique_vals <- sort(unique(na.omit(geo_clipped$cd_geol)))
  val_map <- setNames(seq_along(unique_vals), unique_vals)
  geo_clipped$raster_val <- val_map[geo_clipped$cd_geol]
  
  # Grilla de referencia para rasterización
  ref = rast(ref_grid)
  res(ref) = cellsize
  
  # Rasterizar
  geo_vect <- vect(geo_clipped)
  # writeVector(geo_vect, file.path(temp_dir, "roi_mask.shp"), overwrite=TRUE)
  rasterized <- rasterize(geo_vect, ref, field = "raster_val", touches = TRUE)
  plot(rasterized)
  # values(rasterized)[is.na(values(rasterized))] <- -9999
  
  # Guardar raster como ASCII Grid
  raster_path <- file.path(output_folder, "geology_class.asc")
  writeRaster(rasterized, raster_path, overwrite = TRUE, NAflag = -9999,
              datatype = "INT4S")
  
  # Guardar LUT
  lut_df <- data.frame(
    GeoClass = names(val_map),
    Value = as.integer(val_map)
  )
  write_csv(lut_df, file.path(temp_dir, "LUT_geo_chile.csv"))
  
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
  
  writeLines(param_lines, file.path(output_folder, "geology_classdefinition.txt"))
  
  if (clean_temp && dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
    cat("🧹 Carpeta 'temp' eliminada.\n")
  }
  
  cat("✅ Proceso completado: raster generado desde shapefile.\n")
}

# === Cargar config ===
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")

geo_file <- config$geo_file
roi_file <- config$roi_file
target_res <- config$cellsize
morph_folder <- config$morph_folder
ref_file = file.path(config$meteo_folder,"pre.nc")
cellsize = config$cellsize

# === Ejecutar función principal ===
process_vector_to_raster(
  geo_file = geo_file,
  roi_file = roi_file,
  target_res = 0.01,  # en grados
  output_folder = morph_folder,
  clean_temp = FALSE,
  ref_grid = ref_file,
  cellsize = cellsize
)
