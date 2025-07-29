#' Get bounding box polygon of an ROI
#'
#' Constructs a polygon from the bounding box of the provided
#' region of interest (ROI).
#'
#' @param roi An [`sf`] object representing the ROI.
#'
#' @return An [`sf`] polygon representing the ROI extent in EPSG:4326.
#' @keywords internal
get_extent <- function(roi) {
  box = st_bbox(roi)
  box.m = as.matrix(rbind(
    c(box[[1]], box[[2]]),
    c(box[[1]], box[[4]]),
    c(box[[3]], box[[4]]),
    c(box[[3]], box[[2]]),
    c(box[[1]], box[[2]])
  ))
  box.pol = st_polygon(list(box.m)) %>% st_sfc(crs = 4326)
}

#' Iterative focal mean filter
#'
#' Applies a 3x3 mean filter recursively `n` times.
#'
#' @param r A [`terra::rast`] object to smooth.
#' @param n Integer. Number of iterations to apply.
#'
#' @return Smoothed [`terra::rast`] object.
#' @keywords internal
focal_repeat = function(r, n){
  if (n != 0){
    r = terra::focal(x = r, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
    return(focal_repeat(r, n-1))
  }else{
    return(r)
  }
}

#' Extract header from an ESRI ASCII file
#'
#' Reads the first six lines of a `.asc` raster file and writes them
#' to `header_morph.txt` in the specified output folder.
#'
#' @param asc_file Path to the input `.asc` file.
#' @param output_folder Directory where the header file will be created.
#'
#' @return Invisibly returns `NULL`. The header file is written to disk.
#' @keywords internal
extract_asc_header <- function(asc_file, output_folder) {
  # Check if input file exists
  if (!file.exists(asc_file)) {
    stop("The input .asc file does not exist.")
  }
  
  # Check if output folder exists, if not, create it
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Read the first 6 lines of the .asc file (header)
  header_lines <- readLines(asc_file, n = 6)
  
  # Define the path to the output header file
  header_file <- file.path(output_folder, "header_morph.txt")
  
  # Write the header to the file
  writeLines(header_lines, con = header_file)
  
  cat("Header successfully written to", header_file, "\n")
}

#' Compute annual mean or sum from a SpatRaster with time dimension or provided dates
#'
#' @param r SpatRaster object.
#' @param fun Aggregation function: "mean" or "sum".
#' @param dates_vector Optional vector of dates (POSIXct or Date) to use when time(r) is NULL.
#' @return SpatRaster of the mean across years.
annual_mean <- function(r, fun = c("mean", "sum"), dates_vector = NULL) {
  # browser()
  fun <- match.arg(fun)
  
  tvec <- terra::time(r)
  
  if (is.null(tvec) || all(is.na(tvec))) {
    if (is.null(dates_vector)) {
      stop("Raster has no time information and no dates_vector provided.")
    }
    tvec <- dates_vector
  }
  
  years <- format(tvec, "%Y")
  idx <- split(seq_along(years), years)
  
  yearly <- lapply(idx, function(i) {
    if (fun == "sum") {
      sum(r[[i]], na.rm = TRUE)
    } else {
      mean(r[[i]], na.rm = TRUE)
    }
  })
  
  yearly <- rast(yearly)
  mean(yearly, na.rm = TRUE)
}


#' Aggregate daily data to monthly values
#'
#' Converts a SpatRaster with daily layers into monthly sums or means depending
#' on `fun`.
#'
#' @param r A [`terra::SpatRaster`] with a daily time dimension.
#' @param fun Character. Either "mean" or "sum" specifying how to aggregate the
#'   days of each month.
#' @return A SpatRaster with one layer per month.
#' @export
daily_to_monthly <- function(r, fun = "mean") {
  # fun <- match.arg(fun)
  tvec <- terra::time(r)
  ym <- format(tvec, "%Y-%m")
  r_ann <- tapp(r, ym, fun = fun, na.rm = TRUE)
  return(r_ann)
}

#' Aggregate monthly data to yearly values
#'
#' Converts a SpatRaster with monthly layers into yearly sums or means depending
#' on `fun`.
#'
#' @param r A [`terra::SpatRaster`] with a monthly time dimension.
#' @param fun Character. Either "mean" or "sum" specifying how to aggregate the
#'   months of each year.
#' @return A SpatRaster with one layer per year.
#' @export
monthly_to_yearly <- function(r, fun = "mean") {
  # fun <- match.arg(fun)
  tvec <- terra::time(r)
  ym<- format(tvec, "%Y")
  r_ann <- tapp(r, y, fun = fun, na.rm = TRUE)
  return(r_ann)
}

#' Extract ROI mean time series
#'
#' Computes the spatial mean of a variable inside a region of interest (ROI)
#' for every time step of `mHM_Fluxes_States.nc`.
#'
#' @param nc_path Path to `mHM_Fluxes_States.nc`.
#' @param var_name Name of the variable to extract.
#' @param roi_file Either a path to a vector file or a [`terra::SpatVector`]
#'   representing the ROI. If polygons are provided, the mean value inside each
#'   polygon is extracted. For points, the value at each point is returned.
#' @param out_file Optional path to write the resulting CSV. If `NULL`, the
#'   data frame is returned invisibly.
#'
#' @return A data frame with columns `time` and `mean`.
#' @examples
#' extract_roi_timeseries("domain/OUT/mHM_Fluxes_States.nc",
#'                        var_name = "Q",
#'                        roi_file = "domain/roi.shp")
#' @export
extract_roi_timeseries <- function(nc_path, var_name, roi_file, out_file = NULL) {
  library(terra)

  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }

  if (inherits(roi_file, "SpatVector")) {
    roi <- roi_file
  } else if (is.character(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  } else {
    stop("`roi_file` must be a path or a SpatVector object")
  }

  r <- rast(nc_path, subds = var_name)
  geom <- unique(geomtype(roi))
  if (length(geom) != 1) {
    stop("ROI must contain a single geometry type")
  }

  if (geom == "points") {
    vals <- extract(r, roi, ID = FALSE)
  } else if (geom %in% c("polygons", "lines")) {
    vals <- extract(r, roi, fun = mean, na.rm = TRUE, ID = FALSE)
  } else {
    stop("Unsupported geometry type: ", geom)
  }

  vals <- t(as.matrix(vals))
  tvec <- terra::time(r)
  df <- data.frame(time = tvec, vals)
  colnames(df)[-1] <- paste0("roi", seq_len(ncol(vals)))
  
  if (!is.null(out_file)) {
    write.csv(df, out_file, row.names = FALSE)
  }
  
  invisible(df)
}

#' Extract ROI mean time series for all variables
#'
#' Computes the spatial mean inside a region of interest (ROI) for
#' every variable available in `mHM_Fluxes_States.nc`.
#'
#' @param nc_path Path to `mHM_Fluxes_States.nc`.
#' @param roi_file Either a path to a vector file or a [`terra::SpatVector`]
#'   representing the ROI. If polygons are provided, the mean value inside each
#'   polygon is extracted. For points, the value at each point is returned.
#' @param out_file Optional path to write the resulting CSV. If `NULL`, the
#'   data frame is returned invisibly.
#'
#' @return A data frame with a `time` column and one column per variable.
#' @examples
#' extract_roi_timeseries_all("domain/OUT/mHM_Fluxes_States.nc",
#'                            roi_file = "domain/roi.shp")
#' @export
extract_roi_timeseries_all <- function(nc_path, roi_file, out_file = NULL) {
  library(terra)
  library(ncdf4)

  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }

  if (inherits(roi_file, "SpatVector")) {
    roi <- roi_file
  } else if (is.character(roi_file)) {
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
  } else {
    stop("`roi_file` must be a path or a SpatVector object")
  }

  geom <- unique(geomtype(roi))
  if (length(geom) != 1) {
    stop("ROI must contain a single geometry type")
  }

  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  vars <- names(nc$var)

  r <- rast(nc_path)
  time_vec <- as.Date(terra::time(r))
  res_list <- list()

  for (v in vars) {
    r <- rast(nc_path, subds = v)
    if (geom == "points") {
      vals <- extract(r, roi, ID = FALSE)
    } else if (geom %in% c("polygons", "lines")) {
      vals <- extract(r, roi, fun = mean, na.rm = TRUE, ID = FALSE)
    } else {
      stop("Unsupported geometry type: ", geom)
    }
    vals <- t(as.matrix(vals))
    colnames(vals) <- paste0(v, "_roi", seq_len(ncol(vals)))
    res_list[[v]] <- vals
  }

  res_mat <- do.call(cbind, res_list)
  df <- data.frame(time = time_vec, res_mat)
  
  if (!is.null(out_file)) {
    write.csv(df, out_file, row.names = FALSE)
  }
  
  invisible(df)
}

#' Mosaic mHM outputs from multiple domains
#'
#' Reads `mHM_Fluxes_States.nc` from each domain, applies the ROI mask defined
#' in the respective `preprocess_config.yaml` and merges all domains into
#' separate NetCDF files, one per variable.
#'
#' @param domains Character vector with paths to the domain folders. If `NULL`,
#'   all folders in the current working directory matching "domain" are used.
#' @param out_dir Folder where the mosaicked NetCDFs will be created.
#' @param vars Character vector of variable names to mosaic. By default the
#'   variables used in `visualize_annual_outputs()` except precipitation and
#'   potential evapotranspiration are processed.
#' @examples
#' mosaic_outputs(domains = c("domain_1", "domain_2"),
#'                out_dir = "domain_chile/OUT")
#' @export
mosaic_outputs <- function(domains = NULL,
                           out_dir = "domain_chile/OUT",
                           vars = c("snowpack","SM_Lall","satSTW","aET","Q",
                                    "SM_L01","SM_L02","SM_L03","SM_L04","SM_L05","SM_L06"),
                           config_name = "preprocess_config.yaml") {
  library(terra)
  library(yaml)
  
  if (is.null(domains)) {
    domains <- dir(pattern = "domain", full.names = TRUE)
  }
  if (length(domains) == 0) {
    stop("No domain folders found")
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  written <- vector("character", length(vars))
  for (v in vars) {
    message("Mosaicking variable: ", v)
    rasters <- list()
    for (d in domains) {
      cfg_path <- file.path(d, config_name)
      if (!file.exists(cfg_path)) {
        stop("Configuration file not found in ", d)
      }
      cfg <- read_yaml(cfg_path)
      nc_file <- file.path(d, cfg$out_folder, "mHM_Fluxes_States.nc")
      if (!file.exists(nc_file)) {
        stop("NetCDF file not found: ", nc_file)
      }
      r <- rast(nc_file, subds = v)
      if (!is.null(cfg$roi_file) && file.exists(cfg$roi_file)) {
        roi <- vect(cfg$roi_file)
        r <- mask(r, roi)
      }
      rasters[[length(rasters) + 1]] <- r
    }
    
    mosaic_r <- do.call(mosaic, rasters)
    out_file <- file.path(out_dir, paste0(v, ".nc"))
    writeCDF(mosaic_r, filename = out_file, varname = v, overwrite = TRUE)
    written[which(vars == v)] <- out_file
  }
  
  invisible(written)
}

#' Mosaic meteorological forcings from multiple domains
#'
#' Reads daily NetCDF files located in the `meteo` folder of each domain and
#' merges them into a single file per variable. The ROI mask defined in each
#' `preprocess_config.yaml` is applied prior to mosaicking.
#'
#' @param domains Character vector with paths to the domain folders. If `NULL`,
#'   all folders in the current working directory matching "domain" are used.
#' @param out_dir Folder where the mosaicked NetCDFs will be created.
#' @param vars Character vector of forcing names to mosaic. Defaults to
#'   `c("pre", "pet", "tmin", "tmax", "tavg")` which correspond to
#'   precipitation, potential evapotranspiration and temperature metrics.
#' @examples
#' mosaic_meteo(domains = c("domain_1", "domain_2"),
#'              out_dir = "domain_chile/OUT")
#' @export
mosaic_meteo <- function(domains = NULL,
                         out_dir = "domain_chile/OUT",
                         vars = c("pre", "pet", "tmin", "tmax", "tavg"),
                         config_name = "preprocess_config.yaml") {
  library(terra)
  library(yaml)
  
  if (is.null(domains)) {
    domains <- dir(pattern = "domain", full.names = TRUE)
  }
  if (length(domains) == 0) {
    stop("No domain folders found")
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  written <- vector("character", length(vars))
  for (v in vars) {
    message("Mosaicking variable: ", v)
    rasters <- list()
    for (d in domains) {
      cfg_path <- file.path(d, config_name)
      if (!file.exists(cfg_path)) {
        stop("Configuration file not found in ", d)
      }
      cfg <- read_yaml(cfg_path)
      meteo_file <- file.path(d, "meteo", paste0(v, ".nc"))
      if (!file.exists(meteo_file)) {
        stop("NetCDF file not found: ", meteo_file)
      }
      r <- rast(meteo_file, subds = v)
      if (!is.null(cfg$roi_file) && file.exists(cfg$roi_file)) {
        roi <- vect(cfg$roi_file)
        r <- mask(r, roi)
      }
      rasters[[length(rasters) + 1]] <- r
    }
    
    mosaic_r <- do.call(mosaic, rasters)
    out_file <- file.path(out_dir, paste0(v, ".nc"))
    writeCDF(mosaic_r, filename = out_file, varname = v, overwrite = TRUE)
    written[which(vars == v)] <- out_file
  }
  
  invisible(written)
}


#' Export mHM outputs by time step
#'
#' Reads `mHM_Fluxes_States.nc` located in the `OUT` folder of a domain and
#' writes the selected variable either as monthly or yearly layers. Results are
#' saved inside `OUT/<var_name>` in the chosen format.
#'
#' @param domain_path Path to the model domain containing the `OUT` folder.
#' @param var_name Name of the variable inside `mHM_Fluxes_States.nc` to export.
#' @param ts Time step of the output. Either "month" or "year".
#' @param roi_mask Logical. Apply the ROI mask to the output.
#' @param roi_file Optional path to a ROI file. Used when `roi_mask` is `TRUE`.
#'   If `NULL`, the `roi_file` defined in `preprocess_config.yaml` is used.
#' @param out.format Output format: "tif" for GeoTIFF or "nc" for NetCDF.
#' @param out.opt Optional path to the output folder. If provided, overrides
#'   the `out_folder` defined in `preprocess_config.yaml`.
#'
#' @return Invisibly returns the paths of the written files.
#' @examples
#' write_output("/path/to/domain", "Q", "year")
#' @export
write_output <- function(domain_path, var_name, ts,
                         roi_mask = TRUE, roi_file = NULL,
                         out.format = c("tif", "nc"), out.opt = NULL,
                         config_name = "preprocess_config.yaml") {
  library(terra)
  library(yaml)
  
  out.format <- match.arg(out.format)
  ts <- match.arg(tolower(ts), c("month", "year"))
  
  config_path <- file.path(domain_path, config_name)
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  config <- read_yaml(config_path)

  out_folder <- if (!is.null(out.opt)) out.opt else config$out_folder
  nc_path <- file.path(domain_path, out_folder, "mHM_Fluxes_States.nc")
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }
  
  r <- rast(nc_path, subds = var_name)
  
  if (roi_mask) {
    if (is.null(roi_file)) {
      roi_file <- config$roi_file
    }
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }
  
  if (ts == "month") {
    if (var_name %in% c("aET","Q")){
      r <- daily_to_monthly(r, fun = "sum")
    }else{
      r <- daily_to_monthly(r, fun = "mean")
    }
  } else if (ts == "year") {
    r <- monthly_to_yearly(r, fun = "mean")
  }
  
  out_dir <- file.path(domain_path, out_folder, var_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (out.format == "tif") {
    fname <- file.path(out_dir, paste0(var_name,"_",ts,".tif"))
    writeRaster(r, fname, overwrite = TRUE)
    paths <- fname
    paths <- vector("character", nlyr(r))
  } else {
    fname <- file.path(out_dir, paste0(var_name,"_",ts,".nc"))
    writeCDF(r, filename = fname, varname = var_name, overwrite = TRUE)
    paths <- fname
  }
  
  invisible(paths)
}

#' Export meteorological forcings by time step
#'
#' Reads daily forcing NetCDF files from the `meteo` folder and aggregates
#' them to monthly or yearly values. The result is written either as GeoTIFF or
#' NetCDF inside `OUT/<var_name>`.
#'
#' @param domain_path Path to the model domain.
#' @param var_name Name of the forcing variable to export.
#' @param ts Time step of the output: "month" or "year".
#' @param roi_mask Logical. Apply the ROI mask to the output.
#' @param roi_file Optional path to a ROI file. Used when `roi_mask` is `TRUE`.
#' @param out.format Output format: "tif" for GeoTIFF or "nc" for NetCDF.
#' @param config_name Name of the configuration file.
#'
#' @return Invisibly returns the paths of the written files.
#' @export
write_clim = function(domain_path, var_name, ts,
                      roi_mask = TRUE, roi_file = NULL,
                      out.format = c("tif", "nc"),
                      config_name = "preprocess_config.yaml"){ 
  out.format <- match.arg(out.format)
  ts <- match.arg(tolower(ts), c("month", "year"))
  
  config_path <- file.path(domain_path, config_name)
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  config <- read_yaml(config_path)
  
  out_folder <- config$out_folder
  meteo_folder = config$meteo_folder
  
  
  nc_path <- file.path(domain_path, meteo_folder, paste0(var_name, ".nc"))
  if (!file.exists(nc_path)) {
    stop("NetCDF not found: ", nc_path)
  }
  
  r <- rast(nc_path, subds = var_name)
  
  if (roi_mask) {
    if (is.null(roi_file)) {
      roi_file <- config$roi_file
    }
    if (!file.exists(roi_file)) {
      stop("ROI file not found: ", roi_file)
    }
    roi <- vect(roi_file)
    r <- mask(r, roi)
  }
  
  if (ts == "month") {
    if (var_name %in% c("pre","pet")){
      r <- daily_to_monthly(r, fun = "sum")
    }else{
      r <- daily_to_monthly(r, fun = "mean")
    }
  } else if (ts == "year") {
    r <- monthly_to_yearly(r, fun = "mean")
  }
  
  out_dir <- file.path(domain_path, out_folder, var_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (out.format == "tif") {
    fname <- file.path(out_dir, paste0(var_name,"_",ts,".tif"))
    writeRaster(r, fname, overwrite = TRUE)
    paths <- fname
    paths <- vector("character", nlyr(r))
  } else {
    fname <- file.path(out_dir, paste0(var_name,"_",ts,".nc"))
    writeCDF(r, filename = fname, varname = var_name, overwrite = TRUE)
    paths <- fname
  }
  
  invisible(paths)
}

#' Get plotting limits for a raster
#'
#' Computes the 10\% and 90\% quantiles of the raster values which are
#' typically used to define the limits of a colour scale.
#'
#' @param r A [`terra::SpatRaster`] object.
#' @param round.digit Integer indicating the number of digits to round the
#'   limits.
#' @param lower.lim Numeric lower quantile used when computing limits.
#' @param upper.lim Numeric upper quantile used when computing limits.
#'
#' @return Numeric vector of length two containing the lower and upper limits.
#' @keywords internal
get_scale_lim = function(r, round.digit = 0, lower.lim = 0.1, upper.lim = 0.9){
  values(r) %>%
    quantile(., c(lower.lim, upper.lim), na.rm = TRUE) %>%
    round(., digits = round.digit)
}

#' Plot raster facets
#'
#' Displays multiple layers of a raster in a faceted layout with a common
#' colour scale.
#'
#' @param r A [`terra::SpatRaster`] object.
#' @param var Character. Title for the legend.
#' @param limits Numeric vector providing the limits of the colour scale.
#' @param facet_col Attribute column used to create the facets.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @keywords internal
plot_raster_facet <- function(r, var, limits, facet_col = "attributes") {
  library(stars)
  library(cowplot)
  color.scale = scale_fill_distiller(type = "seq", palette = 16, direction = 1,
                                     na.value = "transparent",
                                     name = var,
                                     limits = limits,
                                     oob = scales::squish)
  stars_r = stars::st_as_stars(r)
  ggplot() +
    geom_stars(data = stars_r) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    facet_wrap(.~get(facet_col))+
    theme(legend.title = element_text(size = 12))
  
}

#' Plot a single raster
#'
#' Creates a map of one raster layer using a sequential palette.
#'
#' @param r A [`terra::SpatRaster`] to plot.
#' @param var Character string used as legend title.
#' @param limits Numeric vector with limits for the colour scale.
#' @param palette Integer palette number passed to `scale_fill_distiller`.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @keywords internal
plot_raster <- function(r, var, limits, palette = 16) {
  color.scale = scale_fill_distiller(type = "seq", palette = palette, direction = 1,
                                     na.value = "transparent",
                                     name = var,
                                     limits = limits,
                                     oob = scales::squish)
  ggplot() +
    geom_stars(data = stars::st_as_stars(r)) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    theme(legend.title = element_text(size = 12))
}


#' Get streamflow table in millimetres. Extracts runoff values from gauge stations.
#'
#' Reads simulated runoff from `mHM_Fluxes_States.nc` together with CAMELS
#' observations and returns a table.
#'
#' @param domain_path Path to the model domain directory.
#' @param crop_to_roi Logical. If `TRUE`, only gauges inside the ROI defined in
#'   `preprocess_config.yaml` are considered.
#'
#' @return A tibble with columns `ID`, `date`, `Q_sim` and `Q_obs`.
#' @export
get_qmm_table <- function(domain_path, crop_to_roi = TRUE,
                          config_name = "preprocess_config.yaml") {
  library(terra)
  library(sf)
  library(yaml)
  library(tidyverse)
  # browser()
  # === Load config and paths ===
  config_path <- file.path(domain_path, config_name)
  config <- read_yaml(config_path)
  
  streamflow_data_file <- config$streamflow_data_file_mm
  gauges_file <- config$fluv_station_file
  roi_file <- config$roi_file
  gauge_folder <- file.path(domain_path, config$gauge_folder)
  nc_path <- file.path(domain_path, config$out_folder, "mHM_Fluxes_States.nc")
  
  # === Read spatial data ===
  roi <- read_sf(roi_file)
  sf_use_s2(FALSE)
  gauges <- read_sf(gauges_file)
  
  # === Filter gauges by ROI if needed ===
  if (crop_to_roi) {
    gauge_list <- as.character(st_intersection(roi, gauges)$ID)
    message(length(gauge_list), " gauge stations found in ROI: ", paste(gauge_list, collapse = ", "))
  } else {
    gauge_list <- as.character(gauges$ID)
    message(length(gauge_list), " gauge stations used: ", paste(gauge_list, collapse = ", "))
  }
  
  gauges_roi <- gauges %>% filter(ID %in% gauge_list)
  
  # === Read Q simulation from NetCDF ===
  nc <- rast(nc_path, subds = "Q")
  names(nc) <- as.character(time(nc))  # Assign date names to layers
  
  basins = "/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/camels_cl_boundaries/camels_cl_boundaries.shp"
  basins = read_sf(basins)
  basins = basins %>% filter(gauge_id %in% gauge_list)
  
  basins_area = basins %>% mutate(area_m2 = area_km2*(1000^2)) %>% pull(area_km2)
  
  val <- terra::extract(nc, basins)
  # val$ID <- basins$gauge_id
  df_sim <- val %>%
    as_tibble() %>%
    group_by(ID) %>% 
    summarise_all(mean, na.rm = TRUE) %>%
    select(-ID) %>% 
    bind_cols(ID = basins$gauge_id, 
              area = basins$area_km2) %>%
    pivot_longer(cols = !c(starts_with("ID"), starts_with("area")), names_to = "date", values_to = "Q_sim") %>%
    mutate(date = as.Date(date))
  
  # === Read CAMELS observation ===
  df_obs <- read_csv(streamflow_data_file) %>%
    select(date, all_of(gauge_list)) %>%
    pivot_longer(cols = -date, names_to = "ID", values_to = "Q_obs") %>%
    mutate(ID = as.numeric(ID), date = as.Date(date))
  
  # === Join simulated and observed ===
  df <- full_join(df_sim, df_obs, by = c("ID", "date"))
  sf_use_s2(TRUE)
  
  return(df)
}

#' Get streamflow table in m\eqn{^3}/s. Reads discharge simulated with the routing
#' option of m
#'
#' Extracts simulated and observed discharge from a NetCDF file (discharge.nc)
#' generated by mHM
#'
#' @param nc_path Path to the NetCDF file containing variables `Qsim_<id>` and
#'   `Qobs_<id>`.
#'
#' @return A tibble with columns `ID`, `date`, `Q_sim` and `Q_obs`.
#' @export
get_qm3s_table <- function(nc_path) {
  library(tidyverse)
  library(ncdf4)
  
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))
  
  # Extract and convert time
  time <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  origin <- str_extract(time_units, "\\d{4}-\\d{2}-\\d{2}")
  dates <- as.Date(time / 24, origin = origin)  # assuming time is in hours
  
  # Get variable names
  var_names <- names(nc$var)
  sim_vars <- var_names[str_detect(var_names, "^Qsim_")]
  obs_vars <- var_names[str_detect(var_names, "^Qobs_")]
  
  # Extract gauge IDs
  sim_ids <- str_remove(sim_vars, "^Qsim_")
  obs_ids <- str_remove(obs_vars, "^Qobs_")
  common_ids <- intersect(sim_ids, obs_ids)
  
  # Function to extract time series
  extract_ts <- function(var) {
    vals <- ncvar_get(nc, var)
    vals[vals == -9999] <- NA
    vals
  }
  
  # Combine all into tibble
  result <- purrr::map_dfr(common_ids, function(id) {
    qsim <- extract_ts(paste0("Qsim_", id))
    qobs <- extract_ts(paste0("Qobs_", id))
    tibble(
      ID = as.numeric(id),
      date = dates,
      Q_sim = as.numeric(qsim),
      Q_obs = as.numeric(qobs)
    )
  })
  
  return(result)
}

#' Calculate hydrological error metrics
#'
#' Computes a range of performance metrics comparing observed and simulated
#' streamflow.
#'
#' @param obs Numeric vector of observations.
#' @param sim Numeric vector of simulations.
#' @param norm Character string passed to `hydroGOF::nrmse` indicating the
#'   normalisation method.
#'
#' @return A tibble with one row containing the calculated metrics.
#' @keywords internal
calculate_error_metrics <- function(obs, sim, norm = "sd") {
  # browser()
  N <- length(obs)
  bias <- mean(sim - obs)
  pbias <- hydroGOF::pbias(sim, obs)
  r2 <- hydroGOF::R2(sim, obs)
  rmse <- hydroGOF::rmse(sim, obs)
  ubrmse <- hydroGOF::ubRMSE(sim, obs)
  nrmse <- hydroGOF::nrmse(sim, obs, norm = norm)
  nse <-hydroGOF::NSE(sim, obs)
  kge_components <- hydroGOF::KGE(sim, obs, method="2009", out.type = "full")
  kge <- kge_components$KGE.value
  alpha <- kge_components$KGE.elements[["Alpha"]]
  beta <- kge_components$KGE.elements[["Beta"]]
  r <- kge_components$KGE.elements[["r"]]
  mae <- hydroGOF::mae(sim, obs)
  hbe <- (sum(sim)-sum(obs))/sum(obs)
  return(tibble(
    kge = kge,
    r2 = r2,
    bias = bias,
    pbias = pbias,
    rmse = rmse,
    ubrmse = ubrmse,
    nrmse = nrmse,
    nse = nse,
    alpha = alpha,
    beta = beta,
    hbe = hbe,
    r = r,
    mae = mae,
    N = N
  ))
}



#' Evaluate streamflow metrics at gauge stations
#'
#' Calculates performance metrics for each gauge station in the supplied data
#' frame and optionally writes the results to disk.
#'
#' @param df Data frame containing at least the columns `ID`, `date`,
#'   `Q_obs_mm`, `Q_sim_mm`, `Q_obs_m3s` and `Q_sim_m3s`.
#' @param config_path Path to `preprocess_config.yaml` with station metadata.
#' @param start_date Start date of the evaluation period.
#' @param end_date End date of the evaluation period.
#' @param out_folder Optional output folder where the table of metrics is
#'   written.
#' @param av_threshold Minimum percentage of available observations required to
#'   include a gauge.
#'
#' @return A tibble with the computed metrics and station coordinates.
#' @export
evaluate_station_metrics <- function(df,
                                     config_path = "preprocess_config.yaml",
                                     start_date = as.Date("1980-01-01"),
                                     end_date = as.Date("2020-12-31"),
                                     out_folder = NULL,
                                     av_threshold = 80) {
  # Read config file
  config <- yaml::read_yaml(config_path)
  # browser()
  calculate_non_na_percentage <- function(df, start_date = as.Date("1980-01-01"), end_date = as.Date("2020-12-31")) {
    df %>%
      filter(date >= start_date & date <= end_date) %>%
      group_by(ID) %>%
      summarise(
        total = n(),
        non_na_Q_obs_mm = sum(!is.na(Q_obs_mm)),
        non_na_Q_sim_mm = sum(!is.na(Q_sim_mm)),
        non_na_Q_obs_m3s = sum(!is.na(Q_obs_m3s)),
        non_na_Q_sim_m3s = sum(!is.na(Q_sim_m3s)),
        .groups = "drop"
      ) %>%
      mutate(
        perc_Q_obs_mm = 100 * non_na_Q_obs_mm / total,
        perc_Q_sim_mm = 100 * non_na_Q_sim_mm / total,
        perc_Q_obs_m3s = 100 * non_na_Q_obs_m3s / total,
        perc_Q_sim_m3s = 100 * non_na_Q_sim_m3s / total
      ) %>%
      select(ID, starts_with("perc_"))
  }
  # Calculate non-NA percentages (assumes the function is defined elsewhere)
  perc <- calculate_non_na_percentage(df)
  
  selected_gauges <- perc %>%
    filter(perc_Q_obs_m3s > av_threshold) %>%
    pull(ID)
  if (length(selected_gauges) == 0){
    selected_gauges = perc$ID
    message("All guages have less streamflow data than %", av_threshold, ", using all data available")
  }
  
  # Filter and compute metrics (assumes calculate_error_metrics() is defined)
  metrics <- df %>%
    filter(!is.na(Q_obs_m3s)) %>%
    filter(date >= start_date, date <= end_date) %>%
    filter(ID %in% selected_gauges) %>%
    group_by(ID) %>%
    do(calculate_error_metrics(.$Q_obs_m3s, .$Q_sim_m3s))
  print(metrics)
  # Read station coordinates from config
  message("Reading config file. Fluvio stations coordinates read at: ", config_path)
  fluv_station <- read_sf(config$fluv_station_file)
  
  # Join coordinates
  metrics <- metrics %>%
    left_join(fluv_station, by = "ID")
  
  # Optionally write output
  if (!is.null(out_folder)) {
    dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
    write_csv(metrics, file.path(out_folder, "eval_metrics_stations.csv"))
  }
  
  return(metrics)
}

#' Plot raster background with station metrics
#'
#' Generates a map showing a raster layer together with station performance
#' metrics coloured by the chosen variable.
#'
#' @param rast A [`terra::SpatRaster`] used as background.
#' @param metrics_sf An [`sf`] object containing the metric values.
#' @param metric Name of the metric column to colour the points.
#' @param metric_name Optional title for the point legend.
#' @param raster_name Optional title for the raster legend.
#' @param raster_palette Viridis palette used for the raster fill.
#' @param point_palette Viridis palette used for the point colours.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @keywords internal
plot_metric_map <- function(rast, metrics_sf, metric = "kge", metric_name = NULL,
                            raster_name = NULL,
                            raster_palette = "viridis", point_palette = "plasma") {
  # scale_fill_viridis_c() and scale_color_viridis_c()
  # support various palettes: "viridis", "magma", "plasma", "inferno", "cividis".
  # Convert raster to stars object
  stars_rast <- st_as_stars(rast)
  if (is.null(metric_name)){
    metric_name = metric
  }
  if (is.null(raster_name)){
    raster_name = varnames(rast)
  }
  # Base plot
  p <- ggplot() +
    geom_stars(data = stars_rast, na.rm = TRUE) +
    scale_fill_viridis_c(option = raster_palette, na.value = NA, name = raster_name) +
    geom_sf(data = metrics_sf, aes_string(color = metric), size = 2, na.rm = TRUE) +
    scale_color_viridis_c(option = point_palette, name = metric_name, na.value = NA) +
    # theme_minimal() +
    theme(legend.position = "right") +
    labs(x = "Longitude", y = "Latitude")
  
  return(p)
}

#' Summarise streamflow record length
#'
#' Calculates the number of observations and time span of a gauge time series.
#'
#' @param df Optional data frame of daily streamflow with columns `date` and
#'   gauge IDs. If `NULL`, a default CAMELS-CL dataset is read.
#' @param id Gauge identifier to analyse. If `NULL` the function returns `NULL`.
#'
#' @return A tibble with columns `n`, `start`, `end` and `years`.
#' @keywords internal
streamflow_data_range <- function(df = NULL, id = NULL) {
  if (is.null(df)){
    df = read_csv("/Volumes/KINGSTON/FONDECYT_CAMILA/mHM_CL/DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/q_m3s_day.csv")
  }
  if (!is.null(id))
    df %>% dplyr::select(date, all_of(id)) %>% 
    drop_na() %>% 
    summarise(n = n(),
              start = min(date),
              end = max(date),
              years = n/365)
}

#' Calculate and store annual mean rasters for key hydrological variables
#'
#' Reads either monthly `.tif` files or extracts from mHM_Fluxes_States.nc
#' and computes annual means, saving them into `annual_mean` subfolder
#' inside `out_folder` as defined in `config_name`.
#'
#' @param domain_path Path to the model domain containing config and data.
#' @param mask_roi Logical. If `TRUE`, apply ROI mask.
#' @param roi_file Optional path to an ROI file.
#' @param config_name Configuration filename (default: preprocess_config.yaml).
#' @export
calculate_annual_mean <- function(domain_path, mask_roi = FALSE, roi_file = NULL,
                                  config_name = "preprocess_config.yaml", process.clim = TRUE,
                                  variables = c("pre", "pet", "snowpack", "SM_Lall", "satSTW", "aET", "Q",
                                                "SWC_L01", "SWC_L02", "SWC_L03", "SWC_L04", "SWC_L05", "SWC_L06")
) {
  library(terra)
  library(yaml)
  library(lubridate)
  
  config_path <- file.path(domain_path, config_name)
  config <- yaml::read_yaml(config_path)
  out_folder <- config$out_folder
  res_file <- file.path(domain_path, out_folder, "mHM_Fluxes_States.nc")
  meteo_folder <- config$meteo_folder
  clim_file <- file.path(domain_path, meteo_folder, "pre.nc")
  
  roi <- NULL
  if (!is.null(roi_file)) {
    roi <- vect(roi_file)
  } else if (mask_roi) {
    roi_path <- config$roi_file
    if (!grepl("^(/|[A-Za-z]:)", roi_path)) {
      roi_path <- file.path(domain_path, roi_path)
    }
    roi <- vect(roi_path)
  }
  
  out_dir <- file.path(domain_path, out_folder, "annual_mean")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (v in variables) {
    if (!process.clim && v %in% c("pre", "pet")) next
    
    tif_path <- file.path(domain_path, out_folder, v, paste0(v, "_month.tif"))
    r <- NULL
    dates_vector <- NULL
    
    if (file.exists(tif_path)) {
      message("Raster ", basename(tif_path), " exists. Creating annual mean from monthly values")
      r <- rast(tif_path)
      tvec <- time(r)
      
      if (is.null(tvec) || all(is.na(tvec))) {
        message("Raster ", basename(tif_path), " has empty or NA time. Recovering from NC...")
        src_file <- if (v %in% c("pre", "pet")) clim_file else res_file
        r_nc <- try(rast(src_file, subds = v), silent = TRUE)
        if (inherits(r_nc, "try-error")) {
          warning("Variable ", v, " not found in NC for fallback.")
          next
        }
        dates_vector <- unique(floor_date(time(r_nc), "month"))
      }
    } else {
      message("Raster ", v, "_month.tif does not exist. Reading from NetCDF...")
      src_file <- if (v %in% c("pre", "pet")) file.path(domain_path, meteo_folder, paste0(v, ".nc")) else res_file
      r <- try(rast(src_file, subds = v), silent = TRUE)
      if (inherits(r, "try-error")) {
        warning("Variable ", v, " not found in sources.")
        next
      }
    }
    
    fun <- if (v %in% c("snowpack", "SM_Lall", "satSTW", 
                        "SWC_L01", "SWC_L02", "SWC_L03", "SWC_L04", "SWC_L05", "SWC_L06")) "mean" else "sum"
    
    r_ann <- annual_mean(r, fun = fun, dates_vector = dates_vector)
    
    if (!is.null(roi)) {
      r_ann <- crop(r_ann, roi) |> mask(roi)
    }
    
    names(r_ann) <- v
    out_file <- file.path(out_dir, paste0(v, "_annual_mean.tif"))
    writeRaster(r_ann, out_file, overwrite = TRUE)
    message("Saved annual mean for ", v, " at ", out_file)
  }
}



#' Visualize annual mean rasters as 2x5 panel figure
#'
#' Reads annual mean rasters from `annual_mean` folder inside `out_folder`
#' defined in `config_name` and generates a PNG figure showing selected variables.
#'
#' @param domain_path Path to the model domain.
#' @param mask_roi Logical. If `TRUE`, apply ROI mask.
#' @param roi_file Optional path to ROI file.
#' @param png_filename Path to save the PNG figure.
#' @param config_name Configuration filename (default: preprocess_config.yaml).
#' @export
visualize_annual_mean <- function(domain_path, mask_roi = FALSE, roi_file = NULL,
                                  png_filename = NULL, config_name = "preprocess_config.yaml",
                                  variables = c("pre", "pet", "snowpack", "SM_Lall", "satSTW", "aET", "Q",
                                                "SWC_L01", "SWC_L02", "SWC_L03", "SWC_L04", "SWC_L05", "SWC_L06")
) {
  library(terra)
  library(yaml)
  # browser()
  config_path <- file.path(domain_path, config_name)
  config <- yaml::read_yaml(config_path)
  out_folder <- config$out_folder
  ann_dir <- file.path(domain_path, out_folder, "annual_mean")
  
  roi <- NULL
  if (!is.null(roi_file)) {
    roi <- vect(roi_file)
  } else if (mask_roi) {
    roi_path <- config$roi_file
    if (!grepl("^(/|[A-Za-z]:)", roi_path)) {
      roi_path <- file.path(domain_path, roi_path)
    }
    roi <- vect(roi_path)
  }
  
  ann_rasters <- list()
  
  for (v in variables) {
    tif_path <- file.path(ann_dir, paste0(v, "_annual_mean.tif"))
    if (file.exists(tif_path)) {
      r <- rast(tif_path)
      if (!is.null(roi)) {
        r <- crop(r, roi) |> mask(roi)
      }
      ann_rasters[[v]] <- r
    } else {
      warning("Annual mean for ", v, " not found.")
    }
  }
  
  # Derived variables
  disp_ann <- ann_rasters[["pre"]] - ann_rasters[["aET"]]
  # ann_rasters[["SWC_total"]] <- Reduce(`+`, ann_rasters[grep("^SWC_L", names(ann_rasters))])
  ann_rasters[["SWC_total"]] <-  ann_rasters[["SWC_L01"]] + ann_rasters[["SWC_L02"]] + ann_rasters[["SWC_L03"]] +
    ann_rasters[["SWC_L04"]] + ann_rasters[["SWC_L05"]] + ann_rasters[["SWC_L06"]]
  flux_range <- range(c(values(ann_rasters[["pre"]]), values(ann_rasters[["pet"]]),
                        values(ann_rasters[["aET"]]), values(ann_rasters[["Q"]]),
                        values(disp_ann)), na.rm = TRUE)
  
  state_range <- range(c(values(ann_rasters[["SWC_total"]]),
                         values(ann_rasters[["snowpack"]]),
                         values(ann_rasters[["satSTW"]])), na.rm = TRUE)
  
  # Plot
  png(png_filename, width = 2000, height = 1000, res = 150)
  par(mfrow = c(2,5), mar = c(4,4,2,5))
  
  plot(ann_rasters[["pre"]], main = "Precipitation", zlim = flux_range)
  plot(ann_rasters[["pet"]], main = "Potential ET", zlim = flux_range)
  plot(ann_rasters[["aET"]], main = "Actual ET", zlim = flux_range)
  plot(ann_rasters[["Q"]], main = "Runoff", zlim = flux_range)
  plot(disp_ann, main = "Pr - ET", zlim = flux_range)
  
  plot(ann_rasters[["SWC_total"]], main = "Soil Moisture", zlim = state_range)
  plot(ann_rasters[["snowpack"]], main = "Snowpack", zlim = state_range)
  plot(ann_rasters[["satSTW"]], main = "GW Level", zlim = state_range)
  plot.new()
  text(0.5, 0.5, "Summary Maps", cex = 1.2)
  
  dev.off()
  message("Figure saved to: ", png_filename)
}

#' Read model domain area
#'
#' Parses `Configfile.log` generated by mHM and extracts the total domain area
#' reported in square kilometres.
#'
#' @param domain_path Path to the model domain directory.
#' @param config_name Name of the configuration file (e.g. preprocess_config.yaml).
#'
#' @return The domain area in km^2 is printed and returned invisibly.
#' @keywords internal
read_domain_area = function(domain_path, config_name){
  # === Load config and paths ===
  config_path <- file.path(domain_path, config_name)
  config <- read_yaml(config_path)
  
  # Read all lines of the file
  file = file.path(domain_path, config$out_folder, "Configfile.log")
  lines <- readLines(file)
  
  # Find the line that contains 'Total[km2]'
  total_area_line <- grep("Total\\[km2\\]", lines, value = TRUE)
  
  # Extract the number using regular expressions
  area <- as.numeric(sub(".*Total\\[km2\\]\\s+", "", total_area_line))
  
  # Print the result
  print(area)

}

#' Merge streamflow tables from one or multiple domains
#'
#' Reads `streamflow_data.csv` produced during post-processing and merges
#' all domains if `use_multiple_domains` is `TRUE`. Erroneous values in the
#' simulated series are filtered out.
#'
#' @param config_name Name of the configuration file used to locate outputs.
#' @return A list with elements `df` (data frame) and `config` (parsed yaml).
#' @export
read_and_merge_streamflow <- function(config_name = "preprocess_config.yaml") {
  df <- tibble()
  if (use_multiple_domains) {
    dir.list <- dir("..", pattern = "domain_zone", full.names = TRUE)
    for (domain in dir.list) {
      config_path <- file.path(domain, config_name)
      config <- read_yaml(config_path)
      x <- read_csv(file.path(domain, config$out_folder, "streamflow", "streamflow_data.csv"))
      df <- bind_rows(df, x)
    }
  } else {
    config_path <- file.path(domain_path, config_name)
    config <- read_yaml(config_path)
    df <- read_csv(file.path(domain_path, config$out_folder, "streamflow", "streamflow_data.csv"))
  }

  fault.ID <- df %>% filter(Q_sim_m3s > 1e6 | Q_sim_m3s < 0) %>% pull(ID) %>% unique()
  df <- df %>% filter(!(ID %in% fault.ID))

  list(df = df, config = config)
}

#' Prepare folders and paths for streamflow graphics
#'
#' Creates the metrics output folder and the `FIGS` directory, returning
#' their paths together with the location of the routing raster used for maps.
#'
#' @param config Configuration list read from YAML.
#' @return List with `out_figs`, `out_folder` and `map_file`.
#' @export
setup_output_paths <- function(config) {
  if (use_multiple_domains) {
    out_figs <- "../FIGS"
    out_folder <- "../domain_Chile/out/metrics"
    map_file <- "../domain_Chile/out/Q.nc"
  } else {
    out_figs <- file.path(domain_path, "FIGS")
    out_folder <- file.path(domain_path, config$out_folder, "metrics")
    map_file <- file.path(domain_path, config$out_folder, "Q", "Q_month.tif")
  }
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)
  list(out_figs = out_figs, out_folder = out_folder, map_file = map_file)
}

#' Compute streamflow performance metrics
#'
#' Wrapper around `evaluate_station_metrics()` using the global workflow
#' configuration. Results are written inside the metrics folder.
#'
#' @param df Data frame returned by `read_and_merge_streamflow()`.
#' @param config Configuration list read from YAML.
#' @param out_folder Output folder for the metrics table.
#' @return Tibble of evaluation metrics.
#' @export
calculate_metrics <- function(df, config, out_folder) {
  evaluate_station_metrics(
    df = df,
    config_path = "preprocess_config.yaml",
    start_date = start_date,
    end_date = end_date,
    out_folder = out_folder,
    av_threshold = av_threshold
  )
}

#' Plot maps of R2 and KGE metrics
#'
#' Generates side-by-side maps showing the spatial distribution of the
#' performance metrics computed for each gauge.
#'
#' @param metrics Data frame with columns `LON`, `LAT`, `r2` and `kge`.
#' @param map_file Raster file used as background.
#' @param out_figs Folder where the figure is written.
#' @param plot_name Filename of the output PNG.
#' @param r2_threshold Optional numeric filter for R2 values.
#' @param kge_threshold Optional numeric filter for KGE values.
#' @export
plot_metrics_maps <- function(metrics, map_file, out_figs,
                              plot_name = "gauge_metrics.png",
                              r2_threshold = NULL,
                              kge_threshold = NULL) {
  dates_raw <- str_sub(names(rast(map_file)), start = 2)
  dates_vec <- as.Date(paste0(dates_raw, ".01"), format = "%Y.%m.%d")
  metrics_sp <- metrics %>% st_as_sf(coords = c("LON", "LAT"), crs = 4326)
  r <- rast(map_file) %>% annual_mean(fun = "sum", dates_vector = dates_vec)

  metrics_r2 <- metrics_sp
  metrics_kge <- metrics_sp
  if (!is.null(r2_threshold)) {
    metrics_r2 <- metrics_r2 %>% filter(r2 >= r2_threshold)
  }
  if (!is.null(kge_threshold)) {
    metrics_kge <- metrics_kge %>% filter(kge >= kge_threshold)
  }

  p1 <- plot_metric_map(r, metrics_r2, metric = "r2", raster_name = "Q [mm]", metric_name = bquote(R^2)) +
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- plot_metric_map(r, metrics_kge, metric = "kge", raster_name = "Q [mm]", metric_name = "KGE") +
    theme(plot.title = element_text(hjust = 0.5))

  p.mat <- ggarrange(p1, p2, nrow = 1, labels = c("a", "b"))
  ggsave(filename = file.path(out_figs, plot_name), plot = p.mat, width = 10, height = 5)
}

#' Boxplot of streamflow metrics
#'
#' Creates a simple boxplot of the selected metric.
#'
#' @param metrics Data frame returned by `calculate_metrics()`.
#' @param metric_name Name of the metric column to plot.
#' @param out_figs Folder to write the plot.
#' @param plot_name Optional filename of the PNG.
#' @export
plot_metrics_boxplot <- function(metrics, metric_name, out_figs, plot_name = NULL) {
  if (is.null(plot_name)) {
    plot_name <- paste0("boxplot_", metric_name, ".png")
  }

  p <- ggplot(metrics, aes(y = .data[[metric_name]])) +
    geom_boxplot(fill = "lightblue", outlier.shape = 21, outlier.fill = "red") +
    labs(y = metric_name, x = NULL, title = paste("Boxplot of", metric_name))

  ggsave(file.path(out_figs, plot_name), plot = p, width = 4, height = 6)
}

#' Cumulative distribution function of a metric
#'
#' Plots the percentage of gauges exceeding a given metric value.
#'
#' @inheritParams plot_metrics_boxplot
#' @export
plot_cdf <- function(metrics, metric_name, out_figs, plot_name = NULL) {
  plot_data <- metrics %>%
    ungroup() %>%
    select(ID, !!sym(metric_name)) %>%
    arrange(desc(!!sym(metric_name))) %>%
    mutate(rank = row_number(), pct_above = 100 * (rank - 1) / n())

  p.cdf <- ggplot(plot_data, aes(x = !!sym(metric_name), y = 100 - pct_above)) +
    geom_line(color = "blue") +
    labs(x = metric_name, y = "% of IDs with metric ≥ x", title = paste("Cumulative distribution of", metric_name)) +
    theme_bw()

  if (is.null(plot_name)) {
    plot_name <- paste0("cummulative_distribution_", metric_name, "_gauges.png")
  }

  ggsave(file.path(out_figs, plot_name), plot = p.cdf, width = 4, height = 4)
}

#' Compare simulated and observed streamflow
#'
#' Generates scatter plots of observed vs simulated runoff with and without
#' routing at daily and monthly scales.
#'
#' @param df Data frame returned by `read_and_merge_streamflow()`.
#' @param out_figs Folder to write the figure.
#' @param plot_name Name of the PNG file.
#' @export
plot_streamflow_comparison <- function(df, out_figs, plot_name = "streamflow_with_and_without_routing.png") {
  labels <- c("m3s" = "Runoff with routing (m3/s)", "mm" = "Runoff w/o routing (mm)")

  p1 <- df %>%
    filter(!is.na(Q_obs_m3s)) %>%
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>%
    mutate(var = str_c(var, "_", source)) %>%
    select(ID, date, var, unit, value) %>%
    pivot_wider(names_from = var, values_from = value) %>%
    ggplot(aes(x = Q_obs, y = Q_sim)) +
    facet_wrap(. ~ unit, ncol = 2, scales = "free", labeller = as_labeller(labels)) +
    geom_point() +
    ggpubr::stat_cor() +
    ggtitle("Daily streamflow")

  df.month <- df %>%
    filter(!is.na(Q_obs_m3s)) %>%
    mutate(date = floor_date(date, "month")) %>%
    group_by(ID, date) %>%
    summarise(count = n(),
              Q_obs_mm = sum(Q_obs_mm),
              Q_sim_mm = sum(Q_sim_mm),
              Q_obs_m3s = mean(Q_obs_m3s),
              Q_sim_m3s = mean(Q_sim_m3s),
              .groups = "drop")

  p2 <- df.month %>%
    filter(count > 28) %>%
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>%
    mutate(var = str_c(var, "_", source)) %>%
    select(ID, date, var, unit, value) %>%
    pivot_wider(names_from = var, values_from = value) %>%
    ggplot(aes(x = Q_obs, y = Q_sim)) +
    facet_wrap(. ~ unit, ncol = 2, scales = "free", labeller = as_labeller(labels)) +
    geom_point() +
    ggpubr::stat_cor() +
    ggtitle("Monthly streamflow")

  p.streamflow <- ggpubr::ggarrange(p1, p2, ncol = 1, labels = c("a", "b"))
  ggsave(filename = file.path(out_figs, plot_name), plot = p.streamflow, width = 8, height = 8)
}

#' Compare routing effects by basin
#'
#' Creates scatter plots for a set of CAMELS-CL basins showing the relation
#' between routed and non-routed runoff volumes.
#'
#' @inheritParams plot_streamflow_comparison
#' @param config Configuration list with the gauge shapefile path.
#' @param basins Vector of CAMELS gauge IDs to include.
#' @param single_gauge Logical. If `TRUE` each basin is compared only against its
#'   own gauge.
#' @export
plot_basin_streamflow_comparison <- function(df, config, out_figs, basins, single_gauge = TRUE,
                                             plot_name = "simulated_streamflow_routing_by_basin.png") {
  camels <- read_sf("../DATA/SHP/Cuencas_CAMELS/CAMELS_CL_v202201/camels_cl_boundaries/camels_cl_boundaries.shp")
  gauges_file <- config$fluv_station_file
  gauges <- read_sf(gauges_file)

  roi <- camels %>% filter(gauge_id %in% basins)
  sf::sf_use_s2(FALSE)
  inter <- st_intersection(roi, gauges) %>%
    as_tibble() %>%
    select(gauge_id, gauge_name, ID)

  df.basins <- df %>%
    left_join(inter, by = "ID") %>%
    filter(!is.na(gauge_id), !is.na(Q_obs_m3s))

  if (single_gauge) {
    df.basins <- df.basins %>% filter(gauge_id == ID)
  }
  area <- read_domain_area(domain_path, config_name)
  df.basins <- df.basins %>% mutate(Q_sim_m3s = 60*60*24*Q_sim_m3s/(area*1000))

  p1 <- df.basins %>%
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>%
    mutate(var = str_c(var, "_", unit)) %>%
    select(gauge_id, gauge_name, ID, date, var, source, value) %>%
    pivot_wider(names_from = var, values_from = value) %>%
    filter(source == "sim") %>%
    ggplot(aes(x = Q_mm, y = Q_m3s)) +
    facet_wrap(. ~ gauge_name, ncol = 3, scales = "free") +
    geom_point() +
    ggpubr::stat_cor() +
    labs(y = "Q with routing [mm]", x = "Q without routing [mm]") +
    ggtitle("Daily streamflow")

  df.basins.month <- df.basins %>%
    mutate(date = floor_date(date, "month")) %>%
    group_by(gauge_id, gauge_name, ID, date) %>%
    summarise(count = n(),
              Q_obs_mm = sum(Q_obs_mm),
              Q_sim_mm = sum(Q_sim_mm),
              Q_obs_m3s = mean(Q_obs_m3s),
              Q_sim_m3s = sum(Q_sim_m3s),
              .groups = "drop")

  p2 <- df.basins.month %>%
    filter(count > 28) %>%
    pivot_longer(cols = starts_with("Q"), names_sep = "_", names_to = c("var", "source", "unit")) %>%
    mutate(var = str_c(var, "_", unit)) %>%
    select(gauge_id, gauge_name, ID, date, var, source, value) %>%
    pivot_wider(names_from = var, values_from = value) %>%
    filter(source == "sim") %>%
    ggplot(aes(x = Q_mm, y = Q_m3s)) +
    facet_wrap(. ~ gauge_name, ncol = 3, scales = "free") +
    geom_point() +
    ggpubr::stat_cor() +
    labs(y = "Q with routing [mm]", x = "Q without routing [mm]") +
    ggtitle("Monthly streamflow")

  p.streamflow <- ggpubr::ggarrange(p1, p2, ncol = 1, labels = c("a", "b"))
  ggsave(filename = file.path(out_figs, plot_name), plot = p.streamflow, width = 12, height = 8)
}
