library(terra)
library(sf)
library(jsonlite)
library(dplyr)
library(readr)

# === Load config ===
config_path <- file.path(domain_path, "preprocess_config.json")
config <- fromJSON(config_path)

geo_file_global <- config$geo_file_global
roi_file <- config$roi_file
morph_folder <- file.path(domain_path, config$morph_folder)
ref_file = file.path(domain_path, config$meteo_folder,"pre.nc")
cellsize = config$cellsize
remove_temp = config$remove_temp

# === Run ===
temp_dir <- file.path(morph_folder, "temp")
dir.create(temp_dir, showWarnings = FALSE)

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

# LUT CSV
lut_df <- data.frame(
  OriginalValue = unique_vals,
  ReclassifiedValue = unname(val_map)
)
write_csv(lut_df, file.path(temp_dir, "lookup_table.csv"))

# Create geology_classdefinition.txt
n_units <- length(val_map)
lines <- c(
  sprintf("nGeo_Formations  %d", n_units),
  "",
  "GeoParam(i)   ClassUnit     Karstic      Description"
)

for (i in 1:n_units) {
  lines <- c(lines, sprintf("%10d %13d %11d      GeoUnit-%d", i, i, 0, i))
}

lines <- c(
  lines,
  "!<-END\n",
  "!***********************************",
  "! NOTES",
  "!***********************************",
  "1 = Karstic",
  "0 = Non-karstic\n",
  "IMPORTANT ::",
  "   Ordering has to be according to the ordering in mhm_parameter.nml",
  "   (namelist: geoparameter)"
)

writeLines(lines, file.path(morph_folder, "geology_classdefinition.txt"))

if (clean_temp && dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  cat("🧹 Carpeta 'temp' eliminada.\n")
}

cat("✅ Proceso completado.\n")