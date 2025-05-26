library(sf)
library(terra)
library(readr)
library(jsonlite)

# === Load config ===
config <- fromJSON("/Users/mhm/Desktop/FONDECYT_CAMILA/mhm_snow/SCRIPTS/preprocess_config.json")

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

